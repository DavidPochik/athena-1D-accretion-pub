//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_fluxes.cpp
//! \brief Calculate hydro/MHD fluxes

// C headers

// C++ headers
#include <algorithm>   // min,max

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"   // reapply floors to face-centered reconstructed states
#include "../field/field.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#include "../gravity/gravity.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../scalars/scalars.hpp"
#include "hydro.hpp"
#include "hydro_diffusion/hydro_diffusion.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// empty array to pass in place of scalars when NSCALARS=0
static AthenaArray<Real> empty;

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes
//! \brief Calculate Hydrodynamic Fluxes using the Riemann solver

void Hydro::CalculateFluxes(AthenaArray<Real> &w, FaceField &b,
                            AthenaArray<Real> &bcc, AthenaArray<Real> &r,
                            const int order) {
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;

  // scalars
  PassiveScalars *ps = pmb->pscalars;
  AthenaArray<Real> &rl_ = (NSCALARS) ? ps->rl_ : empty;
  AthenaArray<Real> &rr_ = (NSCALARS) ? ps->rr_ : empty;
  AthenaArray<Real> &rlb_ = (NSCALARS) ? ps->rlb_ : empty;
  AthenaArray<Real> &rl3d_ = (NSCALARS) ? ps->rl3d_ : empty;
  AthenaArray<Real> &rr3d_ = (NSCALARS) ? ps->rr3d_ : empty;
  AthenaArray<Real> &s_flux_fc = (NSCALARS) ? ps->scr1_nkji_ : empty;
  AthenaArray<Real> &s_laplacian_all_fc = (NSCALARS) ? ps->scr2_nkji_ : empty;

  // b,bcc are passed as fn parameters becausse clients may want to pass different bcc1,
  // b1, b2, etc., but the remaining members of the Field class are accessed directly via
  // pointers because they are unique. NOTE: b, bcc are nullptrs if no MHD.
#if MAGNETIC_FIELDS_ENABLED
  // used only to pass to (up-to) 2x RiemannSolver() calls per dimension:
  // x1:
  AthenaArray<Real> &b1 = b.x1f, &w_x1f = pmb->pfield->wght.x1f,
                  &e3x1 = pmb->pfield->e3_x1f, &e2x1 = pmb->pfield->e2_x1f;
  // x2:
  AthenaArray<Real> &b2 = b.x2f, &w_x2f = pmb->pfield->wght.x2f,
                  &e1x2 = pmb->pfield->e1_x2f, &e3x2 = pmb->pfield->e3_x2f;
  // x3:
  AthenaArray<Real> &b3 = b.x3f, &w_x3f = pmb->pfield->wght.x3f,
                  &e1x3 = pmb->pfield->e1_x3f, &e2x3 = pmb->pfield->e2_x3f;
#endif
  AthenaArray<Real> &flux_fc = scr1_nkji_;
  AthenaArray<Real> &laplacian_all_fc = scr2_nkji_;

  //--------------------------------------------------------------------------------------
  // i-direction

  AthenaArray<Real> &x1flux = flux[X1DIR];
  AthenaArray<Real> &s_x1flux = ps->s_flux[X1DIR];

  // set the loop limits
  jl = js, ju = je, kl = ks, ku = ke;
  if (MAGNETIC_FIELDS_ENABLED || order == 4) {
    if (pmb->block_size.nx2 > 1) {
      if (pmb->block_size.nx3 == 1) // 2D
        jl = js-1, ju = je+1, kl = ks, ku = ke;
      else // 3D
        jl = js-1, ju = je+1, kl = ks-1, ku = ke+1;
    }
  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // reconstruct L/R states
      if (order == 1) {
        pmb->precon->DonorCellX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
        if (NSCALARS) pmb->precon->DonorCellX1(k, j, is-1, ie+1, r, rl_, rr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, w, bcc, wl_, wr_, r, rl_, rr_);
      } else {
        pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, w, bcc, wl_, wr_, r, rl_,rr_);
      }

      pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
      RiemannSolver(k, j, is, ie+1, IVX, wl_, wr_, x1flux, dxw_, rl_, rr_, s_x1flux);
#else  // MHD:
      // x1flux(IBY) = (v1*b2 - v2*b1) = -EMFZ
      // x1flux(IBZ) = (v1*b3 - v3*b1) =  EMFY
      RiemannSolver(k, j, is, ie+1, IVX, b1, wl_, wr_, x1flux, e3x1, e2x1, w_x1f, dxw_,
                    rl_, rr_, s_x1flux);
#endif

      if (order == 4) {
        for (int n=0; n<NWAVE; n++) {
          for (int i=is; i<=ie+1; i++) {
            wl3d_(n,k,j,i) = wl_(n,i);
            wr3d_(n,k,j,i) = wr_(n,i);
          }
        }
        for (int n=0; n<NSCALARS; n++) {
          for (int i=is; i<=ie+1; i++) {
            rl3d_(n,k,j,i) = rl_(n,i);
            rr3d_(n,k,j,i) = rr_(n,i);
          }
        }
      }
    }
  }

  if (order == 4) {
    // TODO(felker): assuming uniform mesh with dx1f=dx2f=dx3f, so this should factor out
    // TODO(felker): also, this may need to be dx1v, since Laplacian is cell-centered
    Real h = pmb->pcoord->dx1f(is);  // pco->dx1f(i); inside loop
    Real C = (h*h)/24.0;

    // construct Laplacian from x1flux
    pmb->pcoord->LaplacianX1All(x1flux, laplacian_all_fc, 0, NHYDRO-1,
                                kl, ku, jl, ju, is, ie+1);
    if (NSCALARS) {
      pmb->pcoord->LaplacianX1All(s_x1flux, s_laplacian_all_fc, 0, NSCALARS-1,
                                  kl, ku, jl, ju, is, ie+1);
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        // Compute Laplacian of primitive Riemann states on x1 faces
        for (int n=0; n<NWAVE; ++n) {
          pmb->pcoord->LaplacianX1(wl3d_, laplacian_l_fc_, n, k, j, is, ie+1);
          pmb->pcoord->LaplacianX1(wr3d_, laplacian_r_fc_, n, k, j, is, ie+1);
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            wl_(n,i) = wl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
            wr_(n,i) = wr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
          }
        }
        // repeat for scalars
        for (int n=0; n<NSCALARS; ++n) {
          pmb->pcoord->LaplacianX1(rl3d_, ps->laplacian_l_fc_, n, k, j, is, ie+1);
          pmb->pcoord->LaplacianX1(rr3d_, ps->laplacian_r_fc_, n, k, j, is, ie+1);
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            rl_(n,i) = rl3d_(n,k,j,i) - C*ps->laplacian_l_fc_(i);
            rr_(n,i) = rr3d_(n,k,j,i) - C*ps->laplacian_r_fc_(i);
          }
        }
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          pmb->peos->ApplyPrimitiveFloors(wl_, rl_, k, j, i);
          pmb->peos->ApplyPrimitiveFloors(wr_, rr_, k, j, i);
        }

        // Compute x1 interface fluxes from face-centered primitive variables
        // TODO(felker): check that e3x1,e2x1 arguments added in late 2017 work here
        pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
        RiemannSolver(k, j, is, ie+1, IVX, wl_, wr_, flux_fc, dxw_, rl_, rr_, s_flux_fc);
#else  // MHD:
        RiemannSolver(k, j, is, ie+1, IVX, b1, wl_, wr_, flux_fc, e3x1, e2x1,
                      w_x1f, dxw_, rl_, rr_, s_flux_fc);
#endif
        // Apply Laplacian of second-order accurate face-averaged flux on x1 faces
        for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
          for (int i=is; i<=ie+1; i++) {
            x1flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
          }
        }
        for (int n=0; n<NSCALARS; ++n) {
#pragma omp simd
          for (int i=is; i<=ie+1; i++)
            s_x1flux(n,k,j,i) = s_flux_fc(n,k,j,i) + C*s_laplacian_all_fc(n,k,j,i);
        }
      }
    }
  } // end if (order == 4)
  //------------------------------------------------------------------------------
  // end x1 fourth-order hydro

  //--------------------------------------------------------------------------------------
  // j-direction

  if (pmb->pmy_mesh->f2) {
    AthenaArray<Real> &x2flux = flux[X2DIR];
    AthenaArray<Real> &s_x2flux = ps->s_flux[X2DIR];
    // set the loop limits
    il = is-1, iu = ie+1, kl = ks, ku = ke;
    if (MAGNETIC_FIELDS_ENABLED || order == 4) {
      if (pmb->block_size.nx3 == 1) // 2D
        kl = ks, ku = ke;
      else // 3D
        kl = ks-1, ku = ke+1;
    }

    for (int k=kl; k<=ku; ++k) {
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js-1, il, iu, w, bcc, wl_, wr_);
        if (NSCALARS) pmb->precon->DonorCellX2(k, js-1, il, iu, r, rl_, rr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, w, bcc, wl_, wr_, r, rl_, rr_);
      } else {
        pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, w, bcc, wl_, wr_, r, rl_, rr_);
      }
      for (int j=js; j<=je+1; ++j) {
        // reconstruct L/R states at j
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, il, iu, w, bcc, wlb_, wr_);
          if (NSCALARS) pmb->precon->DonorCellX2(k, j, il, iu, r, rlb_, rr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, il, iu, w, bcc, wlb_, wr_, r, rlb_, rr_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, j, il, iu, w, bcc, wlb_, wr_, r, rlb_,rr_);
        }

        pmb->pcoord->CenterWidth2(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
        RiemannSolver(k, j, il, iu, IVY, wl_, wr_, x2flux, dxw_, rl_, rr_, s_x2flux);
#else  // MHD:
        // flx(IBY) = (v2*b3 - v3*b2) = -EMFX
        // flx(IBZ) = (v2*b1 - v1*b2) =  EMFZ
        RiemannSolver(k, j, il, iu, IVY, b2, wl_, wr_, x2flux, e1x2, e3x2, w_x2f, dxw_,
                      rl_, rr_, s_x2flux);
#endif

        if (order == 4) {
          for (int n=0; n<NWAVE; n++) {
            for (int i=il; i<=iu; i++) {
              wl3d_(n,k,j,i) = wl_(n,i);
              wr3d_(n,k,j,i) = wr_(n,i);
            }
          }
          for (int n=0; n<NSCALARS; n++) {
            for (int i=il; i<=iu; i++) {
              rl3d_(n,k,j,i) = rl_(n,i);
              rr3d_(n,k,j,i) = rr_(n,i);
            }
          }
        }

        // swap the arrays for the next step
        wl_.SwapAthenaArray(wlb_);
        if (NSCALARS) rl_.SwapAthenaArray(rlb_);
      }
    }
    if (order == 4) {
      // TODO(felker): assuming uniform mesh with dx1f=dx2f=dx3f, so factor this out
      // TODO(felker): also, this may need to be dx2v, since Laplacian is cell-centered
      Real h = pmb->pcoord->dx2f(js);  // pco->dx2f(j); inside loop
      Real C = (h*h)/24.0;

      // construct Laplacian from x2flux
      pmb->pcoord->LaplacianX2All(x2flux, laplacian_all_fc, 0, NHYDRO-1,
                                  kl, ku, js, je+1, il, iu);
      if (NSCALARS) {
        pmb->pcoord->LaplacianX2All(s_x2flux, s_laplacian_all_fc, 0, NSCALARS-1,
                                    kl, ku, js, je+1, il, iu);
      }

      // Approximate x2 face-centered primitive Riemann states
      for (int k=kl; k<=ku; ++k) {
        for (int j=js; j<=je+1; ++j) {
          // Compute Laplacian of primitive Riemann states on x2 faces
          for (int n=0; n<NWAVE; ++n) {
            pmb->pcoord->LaplacianX2(wl3d_, laplacian_l_fc_, n, k, j, il, iu);
            pmb->pcoord->LaplacianX2(wr3d_, laplacian_r_fc_, n, k, j, il, iu);
#pragma omp simd
            for (int i=il; i<=iu; ++i) {
              wl_(n,i) = wl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
              wr_(n,i) = wr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
            }
          }
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            pmb->peos->ApplyPrimitiveFloors(wl_, rl_, k, j, i);
            pmb->peos->ApplyPrimitiveFloors(wr_, rr_, k, j, i);
          }
          for (int n=0; n<NSCALARS; ++n) {
            pmb->pcoord->LaplacianX2(rl3d_, ps->laplacian_l_fc_, n, k, j, il, iu);
            pmb->pcoord->LaplacianX2(rr3d_, ps->laplacian_r_fc_, n, k, j, il, iu);
#pragma omp simd
            for (int i=il; i<=iu; ++i) {
              rl_(n,i) = rl3d_(n,k,j,i) - C*ps->laplacian_l_fc_(i);
              rr_(n,i) = rr3d_(n,k,j,i) - C*ps->laplacian_r_fc_(i);
              pmb->peos->ApplyPassiveScalarFloors(rl_, k, j, i);
              pmb->peos->ApplyPassiveScalarFloors(rr_, k, j, i);
            }
          }

          // Compute x2 interface fluxes from face-centered primitive variables
          // TODO(felker): check that e1x2,e3x2 arguments added in late 2017 work here
          pmb->pcoord->CenterWidth2(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
          RiemannSolver(k, j, il, iu, IVY, wl_, wr_, flux_fc, dxw_, rl_, rr_, s_flux_fc);
#else  // MHD:
          RiemannSolver(k, j, il, iu, IVY, b2, wl_, wr_, flux_fc, e1x2, e3x2,
                        w_x2f, dxw_, rl_, rr_, s_flux_fc);
#endif

          // Apply Laplacian of second-order accurate face-averaged flux on x1 faces
          for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
            for (int i=il; i<=iu; i++) {
              x2flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
            }
          }
          for (int n=0; n<NSCALARS; ++n) {
#pragma omp simd
            for (int i=il; i<=iu; i++)
              s_x2flux(n,k,j,i) = s_flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
            }
        }
      }
    } // end if (order == 4)
  }

  //--------------------------------------------------------------------------------------
  // k-direction

  if (pmb->pmy_mesh->f3) {
    AthenaArray<Real> &x3flux = flux[X3DIR];
    AthenaArray<Real> &s_x3flux = ps->s_flux[X3DIR];
    // set the loop limits
    il = is, iu = ie, jl = js, ju = je;
    if (MAGNETIC_FIELDS_ENABLED || order == 4) {
      il = is-1, iu = ie+1, jl = js-1, ju = je+1;
    }

    for (int j=jl; j<=ju; ++j) { // this loop ordering is intentional
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
        if (NSCALARS) pmb->precon->DonorCellX3(ks-1, j, il, iu, r, rl_, rr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, w, bcc, wl_, wr_, r, rl_, rr_);
      } else {
        pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, w, bcc, wl_, wr_, r, rl_, rr_);
      }
      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k
        if (order == 1) {
          pmb->precon->DonorCellX3(k, j, il, iu, w, bcc, wlb_, wr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(k, j, il, iu, w, bcc, wlb_, wr_, r, rlb_, rr_);
        } else {
          pmb->precon->PiecewiseParabolicX3(k, j, il, iu, w, bcc, wlb_, wr_, r, rlb_,rr_);
        }

        pmb->pcoord->CenterWidth3(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
        RiemannSolver(k, j, il, iu, IVZ, wl_, wr_, x3flux, dxw_, rl_, rr_, s_x3flux);
#else  // MHD:
        // flx(IBY) = (v3*b1 - v1*b3) = -EMFY
        // flx(IBZ) = (v3*b2 - v2*b3) =  EMFX
        RiemannSolver(k, j, il, iu, IVZ, b3, wl_, wr_, x3flux, e2x3, e1x3, w_x3f, dxw_,
                      rl_, rr_, s_x3flux);
#endif
        if (order == 4) {
          for (int n=0; n<NWAVE; n++) {
            for (int i=il; i<=iu; i++) {
              wl3d_(n,k,j,i) = wl_(n,i);
              wr3d_(n,k,j,i) = wr_(n,i);
            }
          }
          for (int n=0; n<NSCALARS; n++) {
            for (int i=il; i<=iu; i++) {
              rl3d_(n,k,j,i) = rl_(n,i);
              rr3d_(n,k,j,i) = rr_(n,i);
            }
          }
        }

        // swap the arrays for the next step
        wl_.SwapAthenaArray(wlb_);
        if (NSCALARS) rl_.SwapAthenaArray(rlb_);
      }
    }
    if (order == 4) {
      // TODO(felker): assuming uniform mesh with dx1f=dx2f=dx3f, so factor this out
      // TODO(felker): also, this may need to be dx3v, since Laplacian is cell-centered
      Real h = pmb->pcoord->dx3f(ks);  // pco->dx3f(j); inside loop
      Real C = (h*h)/24.0;

      // construct Laplacian from x3flux
      pmb->pcoord->LaplacianX3All(x3flux, laplacian_all_fc, 0, NHYDRO-1,
                                  ks, ke+1, jl, ju, il, iu);
      if (NSCALARS) {
        pmb->pcoord->LaplacianX3All(s_x3flux, s_laplacian_all_fc, 0, NSCALARS-1,
                                    ks, ke+1, jl, ju, il, iu);
      }

      // Approximate x3 face-centered primitive Riemann states
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=jl; j<=ju; ++j) {
          // Compute Laplacian of primitive Riemann states on x3 faces
          for (int n=0; n<NWAVE; ++n) {
            pmb->pcoord->LaplacianX3(wl3d_, laplacian_l_fc_, n, k, j, il, iu);
            pmb->pcoord->LaplacianX3(wr3d_, laplacian_r_fc_, n, k, j, il, iu);
#pragma omp simd
            for (int i=il; i<=iu; ++i) {
              wl_(n,i) = wl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
              wr_(n,i) = wr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
            }
          }
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            pmb->peos->ApplyPrimitiveFloors(wl_, rl_, k, j, i);
            pmb->peos->ApplyPrimitiveFloors(wr_, rr_, k, j, i);
          }
          for (int n=0; n<NSCALARS; ++n) {
            pmb->pcoord->LaplacianX3(rl3d_, laplacian_l_fc_, n, k, j, il, iu);
            pmb->pcoord->LaplacianX3(rr3d_, laplacian_r_fc_, n, k, j, il, iu);
#pragma omp simd
            for (int i=il; i<=iu; ++i) {
              rl_(n,i) = rl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
              rr_(n,i) = rr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
              pmb->peos->ApplyPassiveScalarFloors(rl_, k, j, i);
              pmb->peos->ApplyPassiveScalarFloors(rr_, k, j, i);
            }
          }

          // Compute x3 interface fluxes from face-centered primitive variables
          // TODO(felker): check that e2x3,e1x3 arguments added in late 2017 work here
          pmb->pcoord->CenterWidth3(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
          RiemannSolver(k, j, il, iu, IVZ, wl_, wr_, flux_fc, dxw_, rl_, rr_, s_flux_fc);
#else  // MHD:
          RiemannSolver(k, j, il, iu, IVZ, b3, wl_, wr_, flux_fc, e2x3, e1x3,
                        w_x3f, dxw_, rl_, rr_, s_flux_fc);
#endif
          // Apply Laplacian of second-order accurate face-averaged flux on x3 faces
          for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
            for (int i=il; i<=iu; i++) {
              x3flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
            }
          }
          for (int n=0; n<NSCALARS; ++n) {
#pragma omp simd
            for (int i=il; i<=iu; i++)
              s_x3flux(n,k,j,i) = s_flux_fc(n,k,j,i) + C*s_laplacian_all_fc(n,k,j,i);
          }
        }
      }
    } // end if (order == 4)
  }

  if (!STS_ENABLED) {
    AddDiffusionFluxes();
    if (NSCALARS) ps->AddDiffusionFluxes();
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes_STS
//! \brief Calculate Hydrodynamic Diffusion Fluxes for STS

void Hydro::CalculateFluxes_STS() {
  AddDiffusionFluxes();
}

void Hydro::AddDiffusionFluxes() {
  Field *pf = pmy_block->pfield;
  // add diffusion fluxes
  if (hdif.hydro_diffusion_defined) {
    if (hdif.nu_iso > 0.0 || hdif.nu_aniso > 0.0)
      hdif.AddDiffusionFlux(hdif.visflx,flux);
    if (NON_BAROTROPIC_EOS) {
      if (hdif.kappa_iso > 0.0 || hdif.kappa_aniso > 0.0)
        hdif.AddDiffusionEnergyFlux(hdif.cndflx,flux);
    }
  }
  if (MAGNETIC_FIELDS_ENABLED && NON_BAROTROPIC_EOS) {
    if (pf->fdif.field_diffusion_defined)
      pf->fdif.AddPoyntingFlux(pf->fdif.pflux);
  }
  return;
}
