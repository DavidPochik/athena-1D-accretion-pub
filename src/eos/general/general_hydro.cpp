
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file general_hydro.cpp
//! \brief implements most but not all of the functions in class EquationOfState
//!        for general EOS hydrodynamics`
//!
//! These functions MUST be implemented in an additional file.
//!
//! Real EquationOfState::PresFromRhoEg(Real rho, Real egas)
//! Real EquationOfState::EgasFromRhoP(Real rho, Real pres)
//! Real EquationOfState::AsqFromRhoP(Real rho, Real pres)
//! void EquationOfState::InitEosConstants(ParameterInput *pin) // can be empty


// C headers

// C++ headers
#include <cmath>   // sqrt()
#include <sstream>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../field/field.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../eos.hpp"

// EquationOfState constructor

EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin) :
  ptable{pmb->pmy_mesh->peos_table},
  pmy_block_{pmb},
  gamma_{pin->GetOrAddReal("hydro", "gamma", 2.)},
  density_floor_{pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024*float_min))}, 
  scalar_floor_{pin->GetOrAddReal("hydro", "sfloor", std::sqrt(1024*float_min))},
  rho_unit_{pin->GetOrAddReal("hydro", "eos_rho_unit", 1.0)},
  inv_rho_unit_{1.0/rho_unit_},
  egas_unit_{pin->GetOrAddReal("hydro", "eos_egas_unit", 1.0)},
  inv_egas_unit_{1.0/egas_unit_},
  vsqr_unit_{egas_unit_/rho_unit_},
  inv_vsqr_unit_{1.0/vsqr_unit_}
  { 
  if (pin->DoesParameterExist("hydro", "efloor")) {
    energy_floor_ = pin->GetReal("hydro", "efloor");
    pressure_floor_ = energy_floor_*(pin->GetOrAddReal("hydro", "gamma", 2.) - 1.);
    pressure_floor_ = pin->GetOrAddReal("hydro", "pfloor", pressure_floor_);
  } else {
    pressure_floor_ = pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024*float_min));
    energy_floor_ = pressure_floor_/(pin->GetOrAddReal("hydro", "gamma", 2.) - 1.);
    pin->SetReal("hydro", "efloor", energy_floor_);
  }
  if (EOS_TABLE_ENABLED) {
    if (!ptable) {
      std::stringstream msg;
      msg << "### FATAL ERROR in EquationOfState::EquationOfState" << std::endl
          << "EOS table data uninitialized. Should be initialized by Mesh." << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  InitEosConstants(pin);
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
//!           const AthenaArray<Real> &prim_old, const FaceField &b,
//!           AthenaArray<Real> &prim, AthenaArray<Real> &bcc, AthenaArray<Real> &s,
//!           const AthenaArray<Real> &r_old, AthenaArray<Real> &r,
//!           Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku)
//! \brief Converts conserved into primitive variables.

void EquationOfState::ConservedToPrimitive(
    AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old, const FaceField &b,
    AthenaArray<Real> &prim, AthenaArray<Real> &bcc, AthenaArray<Real> &s,
    const AthenaArray<Real> &r_old, AthenaArray<Real> &r,
    Coordinates *pco, int il,int iu, int jl,int ju, int kl,int ku) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real& u_d  = cons(IDN,k,j,i);
        Real& u_m1 = cons(IM1,k,j,i);
        Real& u_m2 = cons(IM2,k,j,i);
        Real& u_m3 = cons(IM3,k,j,i);
        Real& u_e  = cons(IEN,k,j,i);

        Real& w_d  = prim(IDN,k,j,i);
        Real& w_vx = prim(IVX,k,j,i);
        Real& w_vy = prim(IVY,k,j,i);
        Real& w_vz = prim(IVZ,k,j,i);
        Real& w_p  = prim(IPR,k,j,i);

        // apply density floor, without changing momentum or energy
        if (u_d < density_floor_) {
          for (int n=0; n<NSCALARS; ++n) {
            Real& s_n = s(n,k,j,i);
            s_n *= u_d / density_floor_;
            s_n = (s_n < scalar_floor_*density_floor_) ?
                  scalar_floor_*density_floor_ : s_n;
          }
          u_d = density_floor_;
        }
        w_d = u_d;

        Real di = 1.0/u_d;
        w_vx = u_m1*di;
        w_vy = u_m2*di;
        w_vz = u_m3*di;

        Real ke = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
        Real s_cell[NSCALARS];
        for (int n=0; n<NSCALARS; ++n) {
          s_cell[n] = s(n,k,j,i);
        }

        // apply pressure/energy floor, correct total energy
        Real efloor = std::max(energy_floor_, GetEgasFloor(u_d, s_cell));
        u_e = (u_e - ke > efloor) ?  u_e : efloor + ke;
        // MSBC: if ke >> energy_floor_ then u_e - ke may still be zero at this point due
        //       to floating point errors/catastrophic cancellation
        w_p = PresFromRhoEg(u_d, u_e - ke, s_cell);

        for (int n=0; n<NSCALARS; ++n) {
          Real& s_n = s(n,k,j,i);
          Real& r_n = r(n,k,j,i);
          // apply passive scalars floor to conserved variable first, then transform:
          // (multi-D fluxes may have caused it to drop below floor)
          s_n = (s_cell[n] < scalar_floor_ * u_d) ? scalar_floor_ * u_d : s_cell[n];
          r_n = s_n * di;
        }
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
//!           const AthenaArray<Real> &bc, AthenaArray<Real> &cons,
//!           const AthenaArray<Real> &r,AthenaArray<Real> &s, Coordinates *pco,
//!           int il, int iu, int jl, int ju, int kl, int ku);
//! \brief Converts primitive variables into conservative variables

void EquationOfState::PrimitiveToConserved(
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bc,
    AthenaArray<Real> &cons, const AthenaArray<Real> &r,
    AthenaArray<Real> &s, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {
  // Force outer-loop vectorization
#pragma omp simd
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      //#pragma omp simd
#pragma novector
      for (int i=il; i<=iu; ++i) {
        Real& u_d  = cons(IDN,k,j,i);
        Real& u_m1 = cons(IM1,k,j,i);
        Real& u_m2 = cons(IM2,k,j,i);
        Real& u_m3 = cons(IM3,k,j,i);
        Real& u_e  = cons(IEN,k,j,i);

        const Real& w_d  = prim(IDN,k,j,i);
        const Real& w_vx = prim(IVX,k,j,i);
        const Real& w_vy = prim(IVY,k,j,i);
        const Real& w_vz = prim(IVZ,k,j,i);
        const Real& w_p  = prim(IPR,k,j,i);

        u_d = w_d;
        u_m1 = w_vx*w_d;
        u_m2 = w_vy*w_d;
        u_m3 = w_vz*w_d;

        Real r_cell[NSCALARS];
        for (int n=0; n<NSCALARS; ++n) {
          r_cell[n] = r(n,k,j,i);
        }
        // cellwise conversion
        u_e = EgasFromRhoP(u_d, w_p, r_cell) + 0.5*w_d*(SQR(w_vx)+SQR(w_vy)+SQR(w_vz));
        for (int n=0; n<NSCALARS; ++n) {
          Real& s_n = s(n,k,j,i);
          r_cell[n] = (r_cell[n] < scalar_floor_) ? scalar_floor_ : r_cell[n];
          s_n = r_cell[n] * u_d;
        }
      }
    }
  }

  return;
}

// overload eos calls without tracers for backward compatibility with pgens
void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
    const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {
  if (NSCALARS > 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::PrimitiveToConserved" << std::endl
        << "When NSCALARS>0, scalars must be passed to this function." << std::endl;
    ATHENA_ERROR(msg);
    return;
  }
  AthenaArray<Real> empty;
  PrimitiveToConserved(prim, bc, cons, empty, empty, pco, il, iu, jl, ju, kl, ku);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::SoundSpeed(Real prim[NHYDRO])
//! \brief returns adiabatic sound speed given vector of primitive variables

Real EquationOfState::SoundSpeed(const Real prim[NHYDRO+NSCALARS]) {
  return std::sqrt(AsqFromRhoP(prim[IDN], prim[IPR], prim+NHYDRO));
}

//---------------------------------------------------------------------------------------
//! \fn void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, int k, int j,
//!                                                 int i)
//! \brief Apply density and pressure floors to reconstructed L/R cell interface states

void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, AthenaArray<Real> &r,
                                           int k, int j, int i) {
  Real& w_d  = prim(IDN,i);
  Real& w_p  = prim(IPR,i);

  // apply density floor
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;
  // apply pressure floor
  Real r_cell[NSCALARS];
  for (int n=0; n<NSCALARS; ++n) {
     r_cell[n] = r(n,i);

  }
  Real pfloor = std::max(pressure_floor_, GetPresFloor(w_d,r_cell));
  w_p = (w_p > pfloor) ?  w_p : pfloor;
  for (int n=0; n<NSCALARS; ++n) {
    Real& r_n  = r(n,i);
    r_n = (r_n > scalar_floor_) ?  r_n : scalar_floor_;
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::ApplyPrimitiveConservedFloors(AthenaArray<Real> &prim,
//!           AthenaArray<Real> &cons, FaceField &b, int k, int j, int i) {
//! \brief Apply pressure (prim) floor and correct energy (cons) (typically after W(U))
void EquationOfState::ApplyPrimitiveConservedFloors(
    AthenaArray<Real> &prim, AthenaArray<Real> &cons, AthenaArray<Real> &bcc,
    AthenaArray<Real> &r, AthenaArray<Real> &s, int k, int j, int i) {
  Real& w_d  = prim(IDN,k,j,i);
  Real& w_p  = prim(IPR,k,j,i);

  Real& u_d  = cons(IDN,k,j,i);
  Real& u_e  = cons(IEN,k,j,i);
  // apply (prim) density floor, without changing momentum or energy
  if (w_d < density_floor_) {
    w_d = density_floor_;
    for (int n=0; n<NSCALARS; ++n) {
      // should be u_d instead but since this is non-rel u_d=w_d
      s(n,k,j,i) = r(n,k,j,i) * w_d;
    }
  }
  // ensure cons density matches
  u_d = w_d;

  Real e_k = 0.5*w_d*(SQR(prim(IVX,k,j,i)) + SQR(prim(IVY,k,j,i)) + SQR(prim(IVZ,k,j,i)));
  // apply pressure floor, correct total energy
  Real s_cell[NSCALARS];
  for (int n=0; n<NSCALARS; ++n) {
     s_cell[n] = s(n,k,j,i);
  }
  Real r_cell[NSCALARS];
  for (int n=0; n<NSCALARS; ++n) {
     r_cell[n] = r(n,k,j,i);
  }
  Real efloor = std::max(energy_floor_, GetEgasFloor(u_d,s_cell));
  Real pfloor = std::max(pressure_floor_, GetPresFloor(w_d,r_cell));
  u_e = std::max(u_e, efloor + e_k); //(w_p > energy_floor_) ? u_e : energy_floor_ + e_k;
  w_p = (w_p > pfloor) ? w_p : pfloor;
  if (NSCALARS) {
    Real di = 1.0/w_d;
    for (int n=0; n<NSCALARS; ++n) {
      Real& s_n  = s(n,k,j,i);
      Real& r_n  = r(n,k,j,i);

      // should be u_d instead but since this is non-rel u_d=w_d
      s_n = (s_n < scalar_floor_*w_d) ?  scalar_floor_*w_d : s_n;

      // this next line, when applied indiscriminately, erases the accuracy gains
      // performed in the 4th order stencils, since <r> != <s>*<1/di>, in general
      r_n = s_n*di;
      // however, if r_n is riding the variable floor, it probably should be applied so
      // that s_n = rho*r_n is consistent (more concerned with conservation than order of
      // accuracy when quantities are floored)
    }
  }

  return;
}

// overload eos calls without tracers for backward compatibility

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::PresFromRhoEg(Real rho, Real egas)
//! \brief Return gas pressure
Real EquationOfState::PresFromRhoEg(Real rho, Real egas) {
  if (NSCALARS > 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::PresFromRhoEg" << std::endl
        << "When NSCALARS>0, scalars (s) must be passed as a 3rd argument" << std::endl;
    ATHENA_ERROR(msg);
    return -1.0;
  }
  return PresFromRhoEg(rho, egas, nullptr);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::EgasFromRhoP(Real rho, Real pres)
//! \brief Return internal energy density
Real EquationOfState::EgasFromRhoP(Real rho, Real pres) {
  if (NSCALARS > 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::EgasFromRhoP" << std::endl
        << "When NSCALARS>0, scalars (r) must be passed as a 3rd argument" << std::endl;
    ATHENA_ERROR(msg);
    return -1.0;
  }
  return EgasFromRhoP(rho, pres, nullptr);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::AsqFromRhoP(Real rho, Real pres)
//! \brief Return adiabatic sound speed squared
Real EquationOfState::AsqFromRhoP(Real rho, Real pres) {
  if (NSCALARS > 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::AsqFromRhoP" << std::endl
        << "When NSCALARS>0, scalars (r) must be passed as a 3rd argument" << std::endl;
    ATHENA_ERROR(msg);
    return -1.0;
  }
  return AsqFromRhoP(rho, pres, nullptr);
}

Real EquationOfState::GetGamma() {
  std::stringstream msg;
  msg << "GetGamma is not defined for general EOS." << std::endl;
  ATHENA_ERROR(msg);
}
