//==============================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//==============================================================================
//! \file accretion.cpp
//  \brief Problem generator for steady-state accretion, developed from a Parker wind model using the QW EoS.
//

#define COMP_DT 1    // Compute hydro time-steps and save to uov

// C/C++ headers
#include <algorithm>  // min, max
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>   // std::numeric_limits<float>::epsilon()
#include <sstream>
#include <string>
// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"
#include "../scalars/scalars.hpp"

// Static Variables
// r=r_PNS equilibrium quantities: inner pressure, inner temperature, inner electron fraction, inner sound speed squared, GZ pressure modifier, Mach number.
static Real T_eq, Ye_eq, dpdd_0;
// Outer boundary quantities
static Real v_outer, rho_outer, ye_outer, T_outer, abar_outer, P_outer;
// Problem quantities: gravitational mass, inner density, outer density, inner temperature, outer velocity, outer pressure, mass accretion rate, inner radius,
// inverse inner radius squared (qdotQW), Avogadro's number, cutoff temperature (qdotQW), Inner B-field strength, B-field parameter (angle?), rho0 integration start time.
static Real mu, rho_0, rho_f, v_f, p_f, r_0, inv_r2, Na, B_0, alpha;
// EoS quantities: initial Lnu*Enu^2 (qdotQW), initial Lnubar*Enubar^2 (qdotQW), initial perturbation time, final perturbation time,
// perturbation time coefficient, Coeff_nu perturbation parameter (qdotQW), Coeff_nubar perturbation parameter (qdotQW),
// electrion neutrino average energy, electron antineutrino average energy, electron neutrino luminosity, electron antineutrino luminosity,
// floor temperature.
static Real Coeff_nu_0, Coeff_nubar_0, t_L_0, t_L_1, t_coeff, dCoeff_nu, dCoeff_nubar, eps_nue, eps_nueb, L_nue, L_nueb, T_floor;
// Machine epsilon.
static const Real float_eps = std::numeric_limits<Real>::epsilon();
// Number of rows in the IC file
static int rows;
// User-ouput variable (UOV) indices
static int IDT1, IDT2, IDT3, IDT4, IDT5, IDT6, IDT7, IDT8, IDT9, IDT10, IDT11, IDT12, IDT13, IDT14, IDT15, IDT16;
// IC file logical parameter
static bool use_IC_file;
// Scalar Ye parameters: Ye index in r & s arrays, Ye index in prim & cons arrays, t index for scalar arrays, number of scalars
static int ye_index, t_index, abar_index, nscalar_size;
// Scalar Ye parameters: inner boundary ye, outer boundary ye, constant inner Ye.
// Scalar Ye pointers: electron number density, electron fraction (passed into EoS)
Real* edens;
Real* efrac;
// NR parameters: guess temperature, guess electron fraction, temperature increment, electron fraction increment,
// tolerance, modifier, maximum iterations, derivative tolerance.
static Real tgs, dtg, dyg, tolsNR, modsNR, maxcsNR;
// cross section, tolerance for finding tau, tolerance for adjusting rho_0
static Real g_a, delta_np, mcsq, sigma_0, tau_v, tau_epsilon, rho_epsilon;
// The number of cycles that must pass before base density is perturbed
static int nthcycle;
// Optical depth
Real global_tau;
// Final active zone velocity
Real Vr_FinalActiveZone;
// Final Ghozt zone pressure
Real P_GZ;

// Vector Potential
static Real A3(const Real x1, const Real x2, const Real x3) {
  Real a3 = 0.5 * B_0 * r_0 * std::pow(r_0/x1,2) *
    (std::sin(x2)*std::cos(alpha) - std::cos(x2)*std::cos(x3)*std::sin(alpha));
  return a3;
}

static Real A2(const Real x1, const Real x2, const Real x3) {
  Real a2 = -0.5 * B_0 * r_0 * std::pow(r_0/x1,2) * std::sin(x3) * std::sin(alpha);
  return a2;
}

static Real A1(const Real x1, const Real x2, const Real x3) {
  Real a1=0.0;
  return a1;
}

// Estimates Fermi integrals up to n=5 for a given eta parameter.
Real fermi(Real n, Real eta){
  if (n==0) {
    Real fermi = log(1.0 + exp(eta));
    return fermi;
  } else if (n==1) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,2) / 2.0 + 1.6449341;
    Real ff = a - std::pow(a,2) / 4.0 + std::pow(a,3) / 9.0 - std::pow(a,4) / 16.0 + std::pow(a,5) / 25.0 - std::pow(a,6) / 36.0 + std::pow(a,7) / 49.0;
    if (eta < 0) {
        Real fermi = ff;
        return fermi;
    } else if (eta == 0) {
        Real fermi = s - ff;
        return fermi;
    } else {
        Real fermi = s - ff;
        return fermi;
    }
  } else if (n==2) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,3) / 3.0 + 3.2898681 * eta;
    Real ff = 2.0 * (a - std::pow(a,2) / 8.0 + std::pow(a,3) / 27.0 - std::pow(a,4) / 64.0 + std::pow(a,5) / 125.0 - std::pow(a,6) / 216.0);
    if (eta<0) {
        Real fermi = ff;
        return fermi;
    } else if (eta==0) {
        Real fermi = s + ff;
        return fermi;
    } else {
        Real fermi = s + ff;
        return fermi;
    }
  } else if (n==3) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,4) / 4.0 + 4.9348022 * std::pow(eta,2) + 11.351273;
    Real ff = 6.0 * (a - std::pow(a,2) / 16.0 + std::pow(a,3) / 81.0 - std::pow(a,4) / 256.0);
    if (eta<0) {
        Real fermi = ff;
        return fermi;
    } else if (eta==0) {
        Real fermi = s - ff;
        return fermi;
    } else {
        Real fermi = s - ff;
        return fermi;
    }
  } else if (n==4) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,5) / 5.0 + 6.5797363 * std::pow(eta,3) + 45.457576 * eta;
    Real ff = 24.0 * (a - std::pow(a,2) / 32.0 + std::pow(a,3) / 243.0);
    if (eta<0) {
        Real fermi = ff;
        return fermi;
    } else if (eta==0) {
        Real fermi = s + ff;
        return fermi;
    } else {
        Real fermi = s + ff;
        return fermi;
    }
  } else if (n==5) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,6) / 6.0 + 8.2246703 * std::pow(eta,4) + 113.64394 * std::pow(eta,2) + 236.53226;
    Real ff = 120.0 * (a - std::pow(a,2) / 64.0 + std::pow(a,3) / 729.0);
    if(eta<0) {
        Real fermi = ff;
        return fermi;
    } else if(eta==0) {
        Real fermi = s - ff;
        return fermi;
    } else {
        Real fermi = s - ff;
        return fermi;
    }
  } else {
      std::cout << " (accretion.cpp) \n n for fermi_approx wasnt between 0 and 5. \n Something might not be correct";
      return 0.0;
  }
}

// Finds the QW EoS electron chemical potential.
Real QWEta(Real rho, Real T, Real Ye) {
  // Returns eta = mu_e / kbT
  // Physical constants
  Real third         = 1.0 / 3.0;
  Real c             = 2.99792458e10;              // Speed of light in cm/s
  Real k             = 1.380649e-16;               // Boltzmann constant in erg/K
  Real mn            = 1.6726e-24;                 // Baryon mass in g
  Real hbar          = 6.62607015e-27/(2.0*PI);    // Reduced Planck's constant in MeV s
  Real c3            = std::pow(k/(hbar*c),3);
  Real eta_den_const = std::pow(6.0,2.0*third);
  Real root3         = std::sqrt(3.0);
  Real eta_den_A     = std::pow(2.0,third) / eta_den_const;
  Real eta_den_B     = 2.0*std::pow(3.0,third)/eta_den_const;

  Real vol       = mn/rho;
  Real T3        = T*T*T;
  Real T4        = T*T3;
  Real a         = c3*std::pow(T,3)*vol*(PI/3.0);
  Real a2        = SQR(a);
  Real a4        = SQR(a2);
  Real a6        = a2 * a4;
  Real y2        = SQR(Ye);
  Real b         = std::sqrt(4.0*a6+27.0*a4*y2);
  Real term      = std::pow(9.0*a2*Ye+root3*b, third);
  Real eta_by_pi = (eta_den_A)*term/a - (eta_den_B)*a/term; // actually eta/pi
  Real eta       = eta_by_pi * PI;
  return eta;
}


Real YeTejas(Real temp, Real ye, Real x, Real Rho, Real etaele) {
  // Returns Sources for Yedot=vdYe/dr in s^-1 units
  // Temperature argument is in MeV
  // Define constants
  Real alpha    = 1.26;                                // Coupling Coefficient
  Real G_F      = 1.16637e-11;                         // Fermi Constant in MeV^-2
  Real Delta    = 1.2935;                              // neutron-proton mass difference in MeV
  Real hbar     = 6.582119569e-22;                     // Planck's constant in MeV s
  Real c        = 2.99792458e10;                       // Light speed in cm/s
  Real erg2MeV  = 6.24151e5;                           // convert erg to MeV
  Real kbol_MeV  = 8.61733326e-11;
  //Luminosity and epsilon
  Real lnu_y      = L_nue;
  Real epsnu_y    = eps_nue;
  Real lnubar_y   = L_nueb;
  Real epsnubar_y = eps_nueb;

  Real ferm5 = fermi(5,0); //118.266;
  Real ferm4 = fermi(4,0); //23.3309;
  Real ferm3 = fermi(3,0); //5.6822;
  Real ferm2 = fermi(2,0); //1.80309;

  Real ktemp_nubar = epsnubar_y  * ferm2 / ferm3;
  Real eavg1_nubar = ktemp_nubar * ferm3 / ferm2;
  Real eavg2_nubar = ktemp_nubar * ktemp_nubar * ferm4 / ferm2;

  Real ktemp_nu = epsnu_y  * ferm2 / ferm3;
  Real eavg1_nu = ktemp_nu * ferm3 / ferm2;
  Real eavg2_nu = ktemp_nu * ktemp_nu * ferm4 / ferm2;

  Real lambda_nue_n  = std::pow(hbar*c,2.0)*(1.0+3.0*alpha*alpha)/(2*PI*PI)*G_F*G_F*lnu_y*1e51*erg2MeV/(r_0*r_0)
                      *(eavg2_nu/eavg1_nu+2.0*Delta+Delta*Delta/eavg1_nu)*(1.0-x); // s^-1 units
  Real lambda_nuebar_p  = std::pow(hbar*c,2.0)*(1.0+3.0*alpha*alpha)/(2*PI*PI)*G_F*G_F*lnubar_y*1e51*erg2MeV/(r_0*r_0)
                      *(eavg2_nubar/eavg1_nubar-2.0*Delta+Delta*Delta/eavg1_nubar)*(1.0-x); // s^-1 units

   // Calculate Eta, send in T in K 
  Real Eta = etaele; // peos->EtaFromRhoT(Rho, temp / kbol_MeV, r_scalars);

  // Reverse reaction rates
  Real f4_0         = fermi(4,0.0);
  Real f4_eta       = fermi(4,Eta);
  Real f4_negeta    = fermi(4,-1.0*Eta);
  Real lambda_ele_p = 0.448*std::pow(temp,5.0)*f4_eta/f4_0;    // s^-1 units (e- p)
  Real lambda_pos_n = 0.448*std::pow(temp,5.0)*f4_negeta/f4_0; // s^-1 units (e+ n)

  // Source term
  return lambda_nue_n+lambda_pos_n-(lambda_nue_n+lambda_pos_n+lambda_nuebar_p+lambda_ele_p)*ye;
}

Real YeTejas_Source(Real temp, Real ye, Real x, Real Rho) {
  // Returns Sources for Yedot=vdYe/dr in s^-1 units
  // Temperature argument is in MeV
  // Define constants
  Real alpha    = 1.26;                                // Coupling Coefficient
  Real G_F      = 1.16637e-11;                         // Fermi Constant in MeV^-2
  Real Delta    = 1.2935;                              // neutron-proton mass difference in MeV
  Real hbar     = 6.582119569e-22;                     // Planck's constant in MeV s
  Real c        = 2.99792458e10;                       // Light speed in cm/s
  Real erg2MeV  = 6.24151e5;                           // convert erg to MeV
  Real kbol_MeV  = 8.61733326e-11;

  //Luminosity and epsilon
  Real lnu_y      = L_nue;
  Real epsnu_y    = eps_nue;
  Real lnubar_y   = L_nueb;
  Real epsnubar_y = eps_nueb;

  Real ferm5 = fermi(5,0); //118.266;
  Real ferm4 = fermi(4,0); //23.3309;
  Real ferm3 = fermi(3,0); //5.6822;
  Real ferm2 = fermi(2,0); //1.80309;

  Real ktemp_nubar = epsnubar_y  * ferm2 / ferm3;
  Real eavg1_nubar = ktemp_nubar * ferm3 / ferm2;
  Real eavg2_nubar = ktemp_nubar * ktemp_nubar * ferm4 / ferm2;

  Real ktemp_nu = epsnu_y  * ferm2 / ferm3;
  Real eavg1_nu = ktemp_nu * ferm3 / ferm2;
  Real eavg2_nu = ktemp_nu * ktemp_nu * ferm4 / ferm2;

  Real lambda_nue_n  = std::pow(hbar*c,2.0)*(1.0+3.0*alpha*alpha)/(2*PI*PI)*G_F*G_F*lnu_y*1e51*erg2MeV/(r_0*r_0)
                      *(eavg2_nu/eavg1_nu+2.0*Delta+Delta*Delta/eavg1_nu)*(1.0-x); // s^-1 units
  Real lambda_nuebar_p  = std::pow(hbar*c,2.0)*(1.0+3.0*alpha*alpha)/(2*PI*PI)*G_F*G_F*lnubar_y*1e51*erg2MeV/(r_0*r_0)
                      *(eavg2_nubar/eavg1_nubar-2.0*Delta+Delta*Delta/eavg1_nubar)*(1.0-x); // s^-1 units

  // Calculate Eta, send in T in K
  Real Eta = QWEta(Rho,temp/8.6173e-11,ye);

  // Reverse reaction rates
  Real f4_0         = fermi(4,0.0);
  Real f4_eta       = fermi(4,Eta);
  Real f4_negeta    = fermi(4,-1.0*Eta);
  Real lambda_ele_p = 0.448*std::pow(temp,5.0)*f4_eta/f4_0;    // s^-1 units (e- p)
  Real lambda_pos_n = 0.448*std::pow(temp,5.0)*f4_negeta/f4_0; // s^-1 units (e+ n)
 
  return lambda_nue_n+lambda_pos_n;

}

Real YeTejas_Sink(Real temp, Real ye, Real x, Real Rho) {
  // Returns Sources for Yedot=vdYe/dr in s^-1 units
  // Temperature argument is in MeV
  // Define constants
  Real alpha    = 1.26;                                // Coupling Coefficient
  Real G_F      = 1.16637e-11;                         // Fermi Constant in MeV^-2
  Real Delta    = 1.2935;                              // neutron-proton mass difference in MeV
  Real hbar     = 6.582119569e-22;                     // Planck's constant in MeV s
  Real c        = 2.99792458e10;                       // Light speed in cm/s
  Real erg2MeV  = 6.24151e5;                           // convert erg to MeV
  Real kbol_MeV  = 8.61733326e-11;

  //Luminosity and epsilon
  Real lnu_y      = L_nue;
  Real epsnu_y    = eps_nue;
  Real lnubar_y   = L_nueb;
  Real epsnubar_y = eps_nueb;

  Real ferm5 = fermi(5,0); //118.266;
  Real ferm4 = fermi(4,0); //23.3309;
  Real ferm3 = fermi(3,0); //5.6822;
  Real ferm2 = fermi(2,0); //1.80309;

  Real ktemp_nubar = epsnubar_y  * ferm2 / ferm3;
  Real eavg1_nubar = ktemp_nubar * ferm3 / ferm2;
  Real eavg2_nubar = ktemp_nubar * ktemp_nubar * ferm4 / ferm2;

  Real ktemp_nu = epsnu_y  * ferm2 / ferm3;
  Real eavg1_nu = ktemp_nu * ferm3 / ferm2;
  Real eavg2_nu = ktemp_nu * ktemp_nu * ferm4 / ferm2;

  Real lambda_nue_n  = std::pow(hbar*c,2.0)*(1.0+3.0*alpha*alpha)/(2*PI*PI)*G_F*G_F*lnu_y*1e51*erg2MeV/(r_0*r_0)
                      *(eavg2_nu/eavg1_nu+2.0*Delta+Delta*Delta/eavg1_nu)*(1.0-x); // s^-1 units
  Real lambda_nuebar_p  = std::pow(hbar*c,2.0)*(1.0+3.0*alpha*alpha)/(2*PI*PI)*G_F*G_F*lnubar_y*1e51*erg2MeV/(r_0*r_0)
                      *(eavg2_nubar/eavg1_nubar-2.0*Delta+Delta*Delta/eavg1_nubar)*(1.0-x); // s^-1 units

  // Calculate Eta, send in T in K
  Real Eta = QWEta(Rho,temp/8.6173e-11,ye);

  // Reverse reaction rates
  Real f4_0         = fermi(4,0.0);
  Real f4_eta       = fermi(4,Eta);
  Real f4_negeta    = fermi(4,-1.0*Eta);
  Real lambda_ele_p = 0.448*std::pow(temp,5.0)*f4_eta/f4_0;    // s^-1 units (e- p)
  Real lambda_pos_n = 0.448*std::pow(temp,5.0)*f4_negeta/f4_0; // s^-1 units (e+ n)
 
  return (lambda_nue_n+lambda_pos_n+lambda_nuebar_p+lambda_ele_p)*ye;

}

// Calculates the quasi equilibrium (nQSE) condition for Ye (equation 77 in QW1996)
// This function exists mostly for diagnostic purposes.
Real Ye_f(Real temp, Real ye, Real x){
  // Returns Ye_f (equation 77 in Qian & Woosley 1996)
  // Temperature argument is in MeV
  // Define constants
  Real alpha   = 1.26;                                // Coupling Coefficient
  Real G_F     = 1.1663787e-5 * std::pow(1.0e3,-2);   // Fermi Constant in MeV^-2
  Real e_nue   = 1.2 * eps_nue;                       // Electron neutrino energy defined in QW1996 appendix in MeV
  Real e_nueb  = 1.2 * eps_nueb;                      // Electron antineutrino energy defined in QW1996 appendix in MeV
  Real Delta   = 1.293;                               // neutron-proton mass difference in MeV
  Real hbar    = 6.582119569e-22;                     // Planck's constant in MeV s
  Real c       = 2.99792458e10;                       // Light speed in cm/s
  Real erg2MeV = 6.24151e5;                           // convert erg to MeV
  Real lum     = std::pow(10,51);                     // 10^51 erg/s unit

  // Forward reaction rates
  Real lambda_nue_n  = std::pow(hbar*c,2.0)*(1.0+3.0*alpha*alpha)/(2*PI*PI)*G_F*G_F*L_nue*erg2MeV*lum/(r_0*r_0)
                      *(e_nue+2.0*Delta+1.2*Delta*Delta/e_nue)*(1.0-x); // s^-1 units
  Real lambda_nueb_p = std::pow(hbar*c,2.0)*(1.0+3.0*alpha*alpha)/(2*PI*PI)*G_F*G_F*L_nueb*erg2MeV*lum/(r_0*r_0)
                      *(e_nueb-2.0*Delta+1.2*Delta*Delta/e_nueb)*(1.0-x); // s^-1 units

  // Ye_f term
  Real out = (lambda_nue_n) / (lambda_nue_n + lambda_nueb_p);
  return out;
}

Real Heating_QW_Modified(Real x, Real Ye) {
  // returns heating term from QW 1996
  Real sigma_0   = 1.76e-44;                                     // Weak interaction cross section in cm^2, taken from Scheck et al 2006
  Real alpha     = 1.26;                                         // Axial coupling coefficient used in QW 1996
  Real mn        = 1.0 / Na;                                     // Baryon mass in g
  Real E_e       = 0.511;                                        // Electrn energy in MeV
  Real F         = fermi(5,0.0) / fermi(3,0.0);                  // Fermi factors from neutrino energy moments
  Real kbT_nue   = eps_nue  * fermi(2,0.0) / fermi(3,0.0);       // Electron neutrino temperature in MeV
  Real kbT_nueb  = eps_nueb * fermi(2,0.0) / fermi(3,0.0);       // Electron antineutrino temperature in MeV
  Real erg2MeV   = 6.24151e5;                                    // convert erg to MeV
  Real MeV2erg   = 1.60218e-6;                                   // convert MeV to erg
  Real prefactor = 9.65 * Na * 1.0e12;                           // QW heating prefactor, Rnu is in units of 10^6 cm
  Real Xn        = 1.0 - Ye;                                     // Neutron fraction
  Real Xp        = Ye;                                           // Proton fraction

  Real H = MeV2erg * prefactor * F * (Xn * L_nue * std::pow(kbT_nue,2.0) + Xp * L_nueb * std::pow(kbT_nueb,2.0)) * (1.0 - x) / std::pow(r_0,2.0); // erg/s/g
  return H;
}

Real Cooling_QW_Modified(Real Rho, Real T, Real Ye, Real etaele) {
  // T is in MeV
  // Returns cooling term from QW 1996
  Real kbol_MeV  = 8.61733326e-11;                               // Boltzmann constant in MeV / K
  Real sigma_0   = 1.76e-44;                                     // Weak interaction cross section in cm^2, taken from Scheck et al 2006
  Real alpha     = 1.26;                                         // Axial coupling coefficient used in QW 1996
  Real mn        = 1.0 / Na;                                     // Baryon mass in g
  Real E_e       = 0.511;                                        // Electrn energy in MeV
  Real h         = 4.1356677e-21;                                // Planck's constant in MeV s
  Real c         = 2.99792458e10;                                // Speed of light in cm/s
  Real prefactor = sigma_0 * (1.0 + 3.0 * std::pow(alpha,2.0)) * PI * c /
                  (std::pow(E_e,2.0) * mn * pow(h * c,3.0));     // Cooling prefactor in MeV^-5 s^-1 g^-1
  Real MeV2erg   = 1.60218e-6;                                   // convert MeV to erg
  Real Xn        = 1.0 - Ye;                                     // Neutron fraction
  Real Xp        = Ye;                                           // Proton fraction

 
  // Calculate Eta, send in T in K 
  Real eta_e = etaele; //peos->EtaFromRhoT(Rho, T / kbol_MeV, r_scalars);
  Real C = MeV2erg * 2.27 * Na / fermi(5,0) * std::pow(T,6.0) * (Xp * fermi(5.0,eta_e) + Xn * fermi(5.0, -1.0 * eta_e)); // erg/s/g;
  return C;
}

Real X_nucleon(Real rho, Real T) {
  // Temperature argument is in MeV
  Real rho_8    = rho / 1.0e8;    // Mass density in 10^8 g/cm^3
  //Real xn       = 828.0 * std::pow(T,9.0/8.0) / std::pow(rho_8,3.0/4.0) * std::exp(-7.074 / T); // Free nucleon mass fraction
  Real xn       = 828.0 * std::pow(T,9.0/8.0) / std::pow(rho_8,3.0/4.0) * std::exp(-8.8 / T); // Free nucleon mass fraction
  return xn;
}

// exists_test1 from https://stackoverflow.com/a/12774387/2275975
inline bool exists (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

//Heating/cooling source term
void heat_cool(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
               const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
               AthenaArray<Real> &cons_scalar);

// Boundary Condition
void InflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh);
void InflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh);

//------------------------------------------------------------------------------
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.

void Mesh::InitUserMeshData(ParameterInput *pin) {
  EnrollUserExplicitSourceFunction(heat_cool);
  // set outer BC
  if (pin->GetString("mesh", "ox1_bc").compare("user") == 0) {
    if (Globals::my_rank == 0)
      printf("Using USER outer BC.\n");
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, InflowOuterX1);
  }
  // set inner BC
  int inner_BC_choice = pin->GetOrAddInteger("problem", "inner_BC_choice", 0);
  if (inner_BC_choice == 0) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, InflowInnerX1);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
        << "Invalid inner BC choice (" << inner_BC_choice << ")." << std::endl;
    ATHENA_ERROR(msg);
  }
  printf("Using USER inner BC %d.\n", inner_BC_choice);

  Real inf     = std::numeric_limits<Real>::infinity();

  // Problem quantities
  mu       = pin->GetReal("problem","GM");
  rho_0    = pin->GetReal("problem","rho_0"); 
  r_0      = pin->GetReal("mesh","x1min");
  inv_r2   = std::pow(r_0, -2);
  Na       = pin->GetReal("problem","Na"); 
  B_0      = pin->GetReal("problem","B_0");
  B_0      = B_0/(std::pow(4.0*PI,0.5)); //convert to Lorentz-Heaviside units
  alpha    = pin->GetReal("problem","alpha");  

  // Outer boundary quantities
  v_outer    = pin->GetReal("problem", "v_outer");
  rho_outer  = pin->GetReal("problem", "rho_outer");
  ye_outer   = pin->GetReal("problem", "ye_outer");
  T_outer    = pin->GetReal("problem", "T_outer");
  abar_outer = pin->GetReal("problem", "abar_outer");
  P_outer    = pin->GetReal("problem", "P_outer");

  // Passive scalar quantities
  ye_index     = pin->GetInteger("hydro","helm_ye_index");
  t_index      = pin->GetInteger("hydro","helm_temp_index");
  abar_index   = pin->GetInteger("hydro","helm_abar_index");
  nscalar_size = pin->GetInteger("hydro","nsSize"); 

  // Quantities used for calculating optical depth
  g_a         = pin->GetReal("problem","Ga");
  delta_np    = pin->GetReal("problem","Delta");
  mcsq        = pin->GetReal("problem","MeCsq");
  sigma_0     = pin->GetReal("problem","Sigma0");
  tau_v       = pin->GetReal("problem","Tau");
  tau_epsilon = pin->GetReal("problem","Tau_Eps");
  rho_epsilon = pin->GetReal("problem","Rho_Eps");
  nthcycle    = pin->GetInteger("problem","ModuloNumber");

  // Single NR parameters
  tgs     = pin->GetReal("problem","Tg_NR");
  dtg     = pin->GetReal("problem","DeltaTg_NR");
  dyg     = pin->GetReal("problem","DeltaYeg_NR");
  tolsNR  = pin->GetReal("problem","Tolerance_NR");
  maxcsNR = pin->GetReal("problem","maxC_NR");
  modsNR  = pin->GetReal("problem","Modifier_NR");

  // final lumonosity/energies
  Real L_nu      = pin->GetReal("problem","L_nu");
  Real L_nubar   = pin->GetReal("problem","L_nubar");
  L_nue          = pin->GetReal("problem","L_nu");
  L_nueb         = pin->GetReal("problem","L_nubar");
  Real eps_nu    = pin->GetReal("problem","eps_nu");
  Real eps_nubar = pin->GetReal("problem","eps_nubar");
  eps_nue        = pin->GetReal("problem","eps_nu");
  eps_nueb       = pin->GetReal("problem","eps_nubar");
  Coeff_nu_0     = L_nu * SQR(eps_nu);
  Coeff_nubar_0  = L_nubar * SQR(eps_nubar);

  // finial lumonosity/energies
  L_nu         = pin->GetOrAddReal("problem","L_nu_f",L_nu);
  L_nubar      = pin->GetOrAddReal("problem","L_nubar_f",L_nubar);
  eps_nu       = pin->GetOrAddReal("problem","eps_nu_f",eps_nu);
  eps_nubar    = pin->GetOrAddReal("problem","eps_nubar_f",eps_nubar);
  Real coeff   = L_nu * SQR(eps_nu);
  dCoeff_nu    = coeff - Coeff_nu_0;
  coeff        = L_nubar * SQR(eps_nubar);
  dCoeff_nubar = coeff - Coeff_nubar_0;
  t_L_0        = pin->GetOrAddReal("problem","l_transition_start",inf);
  t_L_1        = pin->GetOrAddReal("problem","l_transition_end",inf);
  if (t_L_1 < inf && t_L_1 <= t_L_0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
        << "l_transition_end <= l_transition_start" << std::endl;
    ATHENA_ERROR(msg);
  }
  t_coeff = 0.5 * PI / (t_L_1 - t_L_0);
  T_floor = pin->GetOrAddReal("hydro", "T_floor", float_eps);

  // Parse IC choice
  std::string file;
  bool has_file = pin->DoesParameterExist("problem", "file");
  if (has_file) {
    file = pin->GetString("problem", "file");
  }
  bool use_IC_specified = pin->DoesParameterExist("problem", "use_IC_file");
  use_IC_file = pin->GetOrAddBoolean("problem", "use_IC_file", has_file);

  if (use_IC_specified && use_IC_file) {
    if (!has_file) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
          << "No IC file specified in input file." << std::endl;
      ATHENA_ERROR(msg);
    }
    if (!exists(file)) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
          << "Specified IC file " << file << "does not exits." << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  if (has_file) {
    if (!exists(file)) {
      use_IC_file = false;
      if (Globals::my_rank == 0) {
        std::cout << "Unable to locate IC file " << file << ", reverting to code IC."
                  << std::endl;
      }
    }
  }

  // Read ICs from data file
  if (use_IC_file) {
    rows        = pin->GetInteger("problem", "rows");
    int cols    = pin->GetInteger("problem", "cols");
    int col_rho = pin->GetInteger("problem", "col_rho");
    int col_v   = pin->GetInteger("problem", "col_v");
    int col_T   = pin->GetInteger("problem", "col_T");
    int col_Ye  = pin->GetInteger("problem", "col_Ye");

    // Prepare arrays to hold profile
    AllocateRealUserMeshDataField(11);

    ruser_mesh_data[0].NewAthenaArray(rows);
    ruser_mesh_data[1].NewAthenaArray(rows);
    ruser_mesh_data[2].NewAthenaArray(rows);
    ruser_mesh_data[3].NewAthenaArray(rows);
    ruser_mesh_data[4].NewAthenaArray(rows);
    ruser_mesh_data[5].NewAthenaArray(rows);
    ruser_mesh_data[6].NewAthenaArray(rows);
    ruser_mesh_data[7].NewAthenaArray(rows);
    ruser_mesh_data[8].NewAthenaArray(rows);
    ruser_mesh_data[9].NewAthenaArray(rows);
    AthenaArray<Real>& r_in{ruser_mesh_data[0]};
    AthenaArray<Real>& rho_in{ruser_mesh_data[1]};
    AthenaArray<Real>& v_in{ruser_mesh_data[2]};
    AthenaArray<Real>& T_in{ruser_mesh_data[3]};
    AthenaArray<Real>& Ye_in{ruser_mesh_data[4]};

    // Data structure for calculating optical depth [IS THIS CORRECT? I FEEL LIKE SOMETHING IS OFF HERE]
    AthenaArray<Real>& lsum{ruser_mesh_data[5]};
    AthenaArray<Real>& oldrho{ruser_mesh_data[6]};
    AthenaArray<Real>& Teq{ruser_mesh_data[7]};
    AthenaArray<Real>& PGZ{ruser_mesh_data[8]};

    // Data stucture for tracking cycle number
    AllocateIntUserMeshDataField(1);
    iuser_mesh_data[0].NewAthenaArray(rows);
    AthenaArray<int>& counter{iuser_mesh_data[0]};

    if (Globals::my_rank == 0)
      std::cout<< "Using IC file: " << file << "\n";

    std::string line;
    std::ifstream stream;
    stream.open(file);
    Real s_vals[cols];

    for (int n = 0; n < rows; ++n) {
      std::getline(stream, line);
      std::string word;
      std::stringstream iss(line);
      int m=0;
      while (iss >> word) {
        s_vals[m]=std::stof(word);
        m=m+1;
      }
      r_in(n)   = s_vals[0];
      rho_in(n) = s_vals[col_rho+1];
      v_in(n)   = s_vals[col_v+1];
      T_in(n)   = s_vals[col_T+1];
      Ye_in(n)  = s_vals[col_Ye+1];
    }
  }
  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {

  if (COMP_DT) {
      int i = 0;

      IDT1 = i;
      IDT2 = i + 1;
      IDT3 = i + 2;
      IDT4 = i + 3;
      IDT5 = i + 4;
      IDT6 = i + 5;
      IDT7 = i + 6;
      IDT8 = i + 7;
      IDT9 = i + 8;
      IDT10 = i + 9;
      IDT11 = i + 10;
      IDT12 = i + 11;
      IDT13 = i + 12;
      IDT14 = i + 13;
      IDT15 = i + 14;
      IDT16 = i + 15;
      i += 16; //takes index of last one

      AllocateUserOutputVariables(i);

      SetUserOutputVariableName(IDT1, "dt1");
      SetUserOutputVariableName(IDT2, "dt2");
      SetUserOutputVariableName(IDT3, "dt3");
      SetUserOutputVariableName(IDT4, "x1flux");
      SetUserOutputVariableName(IDT5, "dflx_vol");
      SetUserOutputVariableName(IDT6, "coord_src1");
      SetUserOutputVariableName(IDT7, "extra1");
      SetUserOutputVariableName(IDT8, "extra2");
      SetUserOutputVariableName(IDT9, "extra3");
      SetUserOutputVariableName(IDT10, "extra4");
      SetUserOutputVariableName(IDT11, "extra5");
      SetUserOutputVariableName(IDT12, "extra6");
      SetUserOutputVariableName(IDT13, "extra7");
      SetUserOutputVariableName(IDT14, "extra8");
      SetUserOutputVariableName(IDT15, "extra9");
      SetUserOutputVariableName(IDT16, "extra10");

  }
}

//------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Parker wind

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  if (use_IC_file) {
    // define references for MeshBlock::ProblemGenerator
    AthenaArray<Real>& r_in{pmy_mesh->ruser_mesh_data[0]};
    AthenaArray<Real>& rho_in{pmy_mesh->ruser_mesh_data[1]};
    AthenaArray<Real>& v_in{pmy_mesh->ruser_mesh_data[2]};
    AthenaArray<Real>& T_in{pmy_mesh->ruser_mesh_data[3]};
    AthenaArray<Real>& Ye_in{pmy_mesh->ruser_mesh_data[4]};

    for (int k=ks; k<=ke; k++) {
      // Real phi = pcoord->x3v(k);
      for (int j=js; j<=je; j++) {
        Real theta = pcoord->x2v(j);
        for (int i=is; i<=ie; i++) {
          Real r  = pcoord->x1v(i);
          Real r0 = 5e6;
          Real mn = 1.6749286e-24;                // baryon mass in g
          Real rho, v, temp, ye;

          int index=0;
          Real min = 1e15;
          Real diff;

          for (int f=0; f<rows; f++) {
            diff = r-r_in(f);
            if(diff>=0.0) {
              if(diff<min) {
                min   = diff;
                index = f;
              }
            }
          }
          //linear interpolation when r values in ICs and simulation are different
          if(r<2.1e6 and rho_0>1.5e12) {
            Real qp = std::pow(r_0/r,20.0);
            rho     = rho_0*qp;
            // Real mdot= 4.0*3.14*r*r*rho_in(index)*v_in(index);
            v = (4.34e-4)*2e33/(4.0*3.14*r*r*rho);
          }
          else {
            rho = rho_in(index)+(r-r_in(index))*(rho_in(index+1)-rho_in(index))
                      /(r_in(index+1)-r_in(index));
            v   = v_in(index)+(r-r_in(index))*(v_in(index+1)-v_in(index))
                    /(r_in(index+1)-r_in(index));
          }

          temp = T_in(index)+(r-r_in(index))*(T_in(index+1)-T_in(index))
                    /(r_in(index+1)-r_in(index));
          ye   = Ye_in(index)+(r-r_in(index))*(Ye_in(index+1)-Ye_in(index))
                    /(r_in(index+1)-r_in(index));

          Real kbol_MeV = 8.61733326e-11;  // Boltzmann constant in MeV / K 
          phydro->u(IDN,k,j,i)        = rho;
          phydro->u(IM1,k,j,i)        = v * rho;
          phydro->u(IM2,k,j,i)        = 0.0;
          phydro->u(IM3,k,j,i)        = 0.0;
          pscalars->s(ye_index,k,j,i) = ye * rho;
          pscalars->s(t_index,k,j,i)  = temp * rho;
          // I put r > 5.0e6 here because X_nucleon acts weird at low radii, so I want to skip over that
          if(X_nucleon(rho, temp * kbol_MeV) < 1.0 && r > 5.0e6) {
            // Abar is 56 at r>r_shock, otherwise it is 1. X_alpha < 1 at r>r_shock usually.
            pscalars->s(abar_index,k,j,i) = 56. * rho;
            pscalars->r(abar_index,k,j,i) = 56.;
          } else {
            pscalars->s(abar_index,k,j,i) = 1. * rho;
            pscalars->r(abar_index,k,j,i) = 1.;
          }
          pscalars->r(ye_index,k,j,i) = ye;
          pscalars->r(t_index,k,j,i)  = temp;
          Real r_scalar[nscalar_size];
          Real s_scalar[nscalar_size];
          for (int ns=0; ns<nscalar_size; ns++) {
            r_scalar[ns] = pscalars->r(ns,k,j,i);
            s_scalar[ns] = pscalars->s(ns,k,j,i);
          }

          if(std::isnan(pscalars->r(ye_index,k,j,i)) || pscalars->r(ye_index,k,j,i) < 1.0e-4) {
            std::cout << "(ProblemGenerator) Ye is nan or super low in value (0):" << std::endl;
            std::cout << "Ye  = " << pscalars->r(ye_index,k,j,i) << std::endl;
            std::cout << "rho = " << phydro->u(IDN,k,j,i) << " g/cm^3" << std::endl;
            std::cout << "T   = " << temp << " K" << std::endl;
          }
          if(std::isnan(temp)) {
            std::cout << "(ProblemGenerator) temp is nan (1):" << std::endl;
            std::cout << "Ye  = " << pscalars->r(ye_index,k,j,i) << std::endl;
            std::cout << "rho = " << phydro->u(IDN,k,j,i) << " g/cm^3" << std::endl;
            std::cout << "T   = " << temp << " K" << std::endl;
          }
          if(std::isnan(phydro->u(IDN,k,j,i))) {
            std::cout << "(ProblemGenerator) rho is nan (2):" << std::endl;
            std::cout << "Ye  = " << pscalars->r(ye_index,k,j,i) << std::endl;
            std::cout << "rho = " << phydro->u(IDN,k,j,i) << " g/cm^3" << std::endl;
            std::cout << "T   = " << temp << " K" << std::endl;
          }
          if (NON_BAROTROPIC_EOS) {
            if (GENERAL_EOS) {
              Real pressure        = peos->PresFromRhoT(rho, temp, r_scalar);
              phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(rho, pressure, r_scalar); 
            }
            phydro->u(IEN,k,j,i) += 0.5 * (SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i))
                                        +  SQR(phydro->u(IM3,k,j,i))) / rho;

          }
        }
      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    // if root processor and zeroth block
    if ((Globals::my_rank == 0) && (lid == 0)){
      std::cout<<"YES ENTER B field\n";
    }
    AthenaArray<Real> a1,a2,a3;
    int nx1 = (ie-is)+1 + 2*(NGHOST);
    int nx2 = (je-js)+1 + 2*(NGHOST);
    int nx3 = (ke-ks)+1 + 2*(NGHOST);
    a1.NewAthenaArray(nx3,nx2,nx1);
    a2.NewAthenaArray(nx3,nx2,nx1);
    a3.NewAthenaArray(nx3,nx2,nx1);

    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          a1(k,j,i) = A1(pcoord->x1v(i), pcoord->x2f(j), pcoord->x3f(k));
          a2(k,j,i) = A2(pcoord->x1f(i), pcoord->x2v(j), pcoord->x3f(k));
          a3(k,j,i) = A3(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3v(k));
        }
      }
    }

    // Initialize interface fields
    AthenaArray<Real> area,len,len_p1;
    area.NewAthenaArray(nx1);
    len.NewAthenaArray(nx1);
    len_p1.NewAthenaArray(nx1);

    // for 1,2,3-D
    int jl=js; int ju=je+1;
    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary

      if ((pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) ||
        (pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge))
        jl=js+1;
      if ((pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) ||
        (pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge))
        ju=je;
      for (int j=jl; j<=ju; ++j) {
        pcoord->Face2Area(k,j,is,ie,area);
        pcoord->Edge3Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = -1.0*(len(i+1)*a3(k,j,i+1) - len(i)*a3(k,j,i))/area(i);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        pcoord->Face3Area(k,j,is,ie,area);
        pcoord->Edge2Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = (len(i+1)*a2(k,j,i+1) - len(i)*a2(k,j,i))/area(i);
        }
      }
    }

    // for 2D and 3D
    if (block_size.nx2 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face1Area(k,j,is,ie+1,area);
          pcoord->Edge3Length(k,j  ,is,ie+1,len);
          pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) = (len_p1(i)*a3(k,j+1,i) - len(i)*a3(k,j,i))/area(i);
          }
        }
      }
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face3Area(k,j,is,ie,area);
          pcoord->Edge1Length(k,j  ,is,ie,len);
          pcoord->Edge1Length(k,j+1,is,ie,len_p1);
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) -= (len_p1(i)*a1(k,j+1,i) - len(i)*a1(k,j,i))/area(i);
          }
        }
      }
    }
    // for 3D only
    if (block_size.nx3 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face1Area(k,j,is,ie+1,area);
          pcoord->Edge2Length(k  ,j,is,ie+1,len);
          pcoord->Edge2Length(k+1,j,is,ie+1,len_p1);
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) -= (len_p1(i)*a2(k+1,j,i) - len(i)*a2(k,j,i))/area(i);
          }
        }
      }
      for (int k=ks; k<=ke; ++k) {
        // reset loop limits for polar boundary
        int jl=js; int ju=je+1;
        if ((pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) ||
          (pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge))
          jl=js+1;
        if ((pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) ||
          (pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge))
          ju=je;
        for (int j=jl; j<=ju; ++j) {
          pcoord->Face2Area(k,j,is,ie,area);
          pcoord->Edge1Length(k  ,j,is,ie,len);
          pcoord->Edge1Length(k+1,j,is,ie,len_p1);
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) += (len_p1(i)*a1(k+1,j,i) - len(i)*a1(k,j,i))/area(i);
          }
        }
      }
    }
    // Calculate cell-centered magnetic field
    AthenaArray<Real> bb;
    bb.NewAthenaArray(3, ke+1, je+1, ie+NGHOST+1);
    pfield->CalculateCellCenteredField(pfield->b, bb, pcoord, is-NGHOST, ie+NGHOST, js,
                                       je, ks, ke);

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real& bcc1 = bb(IB1,k,j,i);
          Real& bcc2 = bb(IB2,k,j,i);
          Real& bcc3 = bb(IB3,k,j,i);

          phydro->u(IEN,k,j,i) += 0.5*(SQR(bcc1)+SQR(bcc2)+SQR(bcc3));
         }
      }
    }
    a1.DeleteAthenaArray();
    a2.DeleteAthenaArray();
    a3.DeleteAthenaArray();
    area.DeleteAthenaArray();
    len.DeleteAthenaArray();
    len_p1.DeleteAthenaArray();
    bb.DeleteAthenaArray();
  } // end if MAGNETIC_FIELDS_ENABLED
}

// Source Terms
void heat_cool(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
               const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
               AthenaArray<Real> &cons_scalar) {
  Real qdot, vdYedr;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real erg2MeV  = 6.24151e5;                    // Erg to MeV conversion factor
        Real r        = pmb->pcoord->x1v(i);          // Radial coordinate
        Real theta    = pmb->pcoord->x2v(i);          // angular coordinate
        Real p        = prim(IPR,k,j,i);              // Pressure
        Real rho      = prim(IDN,k,j,i);              // Density
        Real kbol_MeV = 8.61733326*std::pow(10,-11);  // Boltzmann constant in MeV / K
        Real Ye_i     = prim_scalar(ye_index,k,j,i);  // Electron fraction before being updated in this function

        Real r_scalar[nscalar_size];
        for (int ns=0; ns<nscalar_size; ns++) {
          r_scalar[ns] = prim_scalar(ns,k,j,i); 
        }
        Real cs = std::sqrt(pmb->peos->AsqFromRhoP(prim(IDN,k,j,i), prim(IPR,k,j,i), r_scalar));
        Real etae = pmb->peos->EtaFromRhoT(prim(IDN,k,j,i), prim_scalar(t_index,k,j,i), r_scalar);
        Real Mach = fabs(prim(IVX,k,j,i)) / cs;
      
        Real temp = prim_scalar(t_index,k,j,i) * kbol_MeV; // MeV
        if(std::isnan(Ye_i)) {
          std::cout << "(heat_cool) NANS (0) (Ye_i)" << std::endl;
          std::cout << "Ye_i  = " << Ye_i << std::endl;
          std::cout << "p     = " << p    << "erg/cm^3" <<std::endl;
          std::cout << "temp  = " << temp << "MeV" <<std::endl;
          std::cout << "rho   = " << rho  << "g/cm^3" <<std::endl;
          std::cout << " " << std::endl;
        }
        if(std::isnan(temp)) {
          std::cout << "(heat_cool) NANS (1) (temp)" << std::endl;
          std::cout << "Ye_i  = " << Ye_i << std::endl;
          std::cout << "p     = " << p    << "erg/cm^3" <<std::endl;
          std::cout << "temp  = " << temp << "MeV" <<std::endl;
          std::cout << "rho   = " << rho  << "g/cm^3" <<std::endl;
          std::cout << " " << std::endl;
        }
        if(rho<0.0) {
          std::cout << "density is negative in heat_cool:" << std::endl;
          std::cout << "r     = " << r/1.0e5 << " km" << std::endl;
          std::cout << "theta = " << theta * 180/PI << " degrees" << std::endl;
        }
        Real X_n     = X_nucleon(rho, temp);
        Real x       = std::sqrt(1.0-(r_0*r_0)/(r*r));
        Real ExpSupp = std::exp(-1.0 * rho / rho_0);
        if(Mach > 1.0) { // Using Mach instead of Xn because Mach clearly greatly exceeds 1 over the shock, whereas Xn is close to 1 across the shock.
          cons_scalar(abar_index,k,j,i) = 56. * prim(IDN,k,j,i); 
          vdYedr = 0.0; //ExpSupp * X_n * YeTejas(temp, Ye_i, x, rho, etae);
          qdot   = 0.0; //ExpSupp * X_n * (Heating_QW_Modified(x, Ye_i) - Cooling_QW_Modified(rho, temp, Ye_i, etae)); // erg/s/g  
        } else {
          vdYedr = ExpSupp * YeTejas(temp, Ye_i, x, rho, etae);
          qdot   = ExpSupp * (Heating_QW_Modified(x, Ye_i) - Cooling_QW_Modified(rho, temp, Ye_i, etae)); // erg/s/g 
          cons_scalar(abar_index,k,j,i) = 1. * prim(IDN,k,j,i);
        } 

        // Find electron fraction and energy increments
        Real dYe = dt*vdYedr;
        Real de  = dt*prim(IDN,k,j,i) * qdot; // erg/cm^-3

        // Update conserved composition and energy fields
        cons_scalar(ye_index,k,j,i) += dYe*prim(IDN,k,j,i); // g cm^-3
        if (cons(IEN,k,j,i)>1.0e50) {
          std::cout << "(helm_eos) before e update:" << std::endl;
          std::cout << "e = " << cons(IEN,k,j,i) << std::endl;
          std::cout << "r = " << r << " cm" << std::endl;
        }
        cons(IEN,k,j,i) += de;                              // erg/cm^3
        if(cons(IEN,k,j,i)>1.0e50) {
          std::cout << "(heat_cool) dt   = " << dt << " s" << std::endl;
          std::cout << "            r    = " << r    << " cm" << std::endl;
          std::cout << "            E    = " << cons(IEN,k,j,i) << std::endl;
          std::cout << "            rho  = " << rho  << " g/cm^3" << std::endl;
          std::cout << "            P    = " << p    << " erg/cm^3" << std::endl;
          std::cout << "         T (MeV) = " << temp    << " MeV" << std::endl;
          std::cout << "           T (K) = " << temp / kbol_MeV << " K" << std::endl;
          std::cout << "            Ye   = " << Ye_i << std::endl;
          std::cout << "           qdot  = " << qdot << " erg/s/g" << std::endl;

        }
        if(std::isnan(cons_scalar(ye_index,k,j,i))) {
          std::cout << "(heat_cool) NANS (2) (cons_scalar)" << std::endl;
          std::cout << "Ye_i  = " << Ye_i << std::endl;
          std::cout << "p     = " << p    << "erg/cm^3" <<std::endl;
          std::cout << "temp  = " << temp << "MeV" <<std::endl;
          std::cout << "rho   = " << rho  << "g/cm^3" <<std::endl;
          std::cout << " " << std::endl;
        }
      }
    }
  }
  return;
}

Real PerturbRho(int Roundncycle, Real ncycle, Real time, Real Vr_FZ, Real rho_active, Real rho_old, Real rho_0) {
  Real rho_new = rho_old;
  // Global tau is calcated in UserWorkInLoop
  Real tau = global_tau;

  // perturb base density every n'th cycle
  if (Roundncycle == ncycle) {
    if (Roundncycle % nthcycle == 0) {
      if (time != 0.0) {
        if ((tau < tau_v * (1.0 + tau_epsilon)) && (tau > tau_v * (1.0 - tau_epsilon))){
          rho_new = rho_old;
          std::cout << "converged" << std::endl;
          std::cout << "ncycle  = " << ncycle << std::endl;
          std::cout << "rho_new = " << rho_new << " g/cm^3" << std::endl;
          std::cout << "active  = " << rho_active << " g/cm^3" << std::endl;
          std::cout << "tau     = " << tau << std::endl;
          std::cout << "time    = " << time << " s" << std::endl;
        } else if (tau < tau_v * (1.0 - tau_epsilon)) {
          // if in an explosion regime, dont update rho
          if(Vr_FZ <= 0.0) {
            rho_new = rho_old * (1.0 + rho_epsilon);
          } else {
            std::cout << "Vr(Last Active Zone) = " << Vr_FZ << " cm/s" << std::endl;
            rho_new = rho_old;
          }
          std::cout << "too small" << std::endl;
          std::cout << "ncycle  = " << ncycle << std::endl;
          std::cout << "rho_new = " << rho_new << " g/cm^3" << std::endl;
          std::cout << "active  = " << rho_active << " g/cm^3" << std::endl;
          std::cout << "tau     = " << tau << std::endl;
          std::cout << "time    = " << time << " s" << std::endl;
        } else if (tau > tau_v * (1.0 + tau_epsilon)) {
          // if in an explosion regime, dont update rho
          if(Vr_FZ <= 0.0) {
            rho_new = rho_old * (1.0 - rho_epsilon);
          } else {
            std::cout << "Vr(Last Active Zone) = " << Vr_FZ << " cm/s" << std::endl;
            rho_new = rho_old;
          }
          std::cout << "too big" << std::endl;
          std::cout << "ncycle  = " << ncycle << std::endl;
          std::cout << "rho_new = " << rho_new << " g/cm^3" << std::endl;
          std::cout << "active  = " << rho_active << " g/cm^3" << std::endl;
          std::cout << "tau     = " << tau << std::endl;
          std::cout << "time    = " << time << " s" << std::endl;
        } else {
          std::cout << "(accretion.cpp) Somehow I've avoided the tau conditional statements. Something isn't right." << std::endl;
        }
      } else {
        std::cout << "t = 0, using default rho0" << std::endl;
        rho_new = rho_0;
        rho_old = rho_new;
      }
    } else {
      if (time == 0.0 || ncycle < 1.0) {
        // Use user-specified base density at t=0
        std::cout << "(1) time = 0s or ncycle < 1, using default rho0" << std::endl;
        std::cout << "time = " << time << " s" << std::endl;
        std::cout << "ncycle = " << ncycle << std::endl;
        rho_new = rho_0;
        rho_old = rho_new;
      } else {
        // Use ghost zone value from previous cycle (this assumes that an HSE profile isn't being set in the ghost zones).
        rho_new = rho_old;
      }
    }
  } else {
    if (time == 0.0 || ncycle < 1.0) {
      // Use user-specified base density at t=0
      std::cout << "(2) time = 0s or ncycle < 1, using default rho0" << std::endl;
      std::cout << "time = " << time << " s" << std::endl;
      std::cout << "ncycle = " << ncycle << std::endl;
 
      rho_new = rho_0;
      rho_old = rho_new;
    } else {
      rho_new = rho_old;
    }
  }
  return rho_new;
}

// Inflow Boundary Condition
void InflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh) {

  AthenaArray<int>& my_data1 = pmb->pmy_mesh->iuser_mesh_data[0];
  AthenaArray<Real>& my_data2 = pmb->pmy_mesh->ruser_mesh_data[6];
  AthenaArray<Real>& my_data3 = pmb->pmy_mesh->ruser_mesh_data[7];

  Real &rho_old{my_data2(0)};
  Real &T_equilibrium{my_data3(0)};
  // Ncycle is updated twice per time cycle because there are two steps in 1 time cycle
  int &ncycle_x2{my_data1(0)};
  if (time == 0.0) {
    ncycle_x2 = 0;
  } else {
    ncycle_x2 += 1;
  }

  // Divide ncycle_x2 by 2 to get the actual cycle number
  Real ncycle = ncycle_x2 / 2.0;
  int Roundncycle = round(ncycle);
 
  // Final zone velocity, used to determine if an explosion has occured
  Real Vr_FZ = Vr_FinalActiveZone;

  // Compute perturbed density based on the value of tau
  Real rho_new = PerturbRho(Roundncycle, ncycle, time, Vr_FZ, prim(IDN,ks,js,is), rho_old, rho_0);
  //std::cout << "(InflowInnerX1) rho_new = " << rho_new << " g/cm^3" << std::endl;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      Real theta = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
        prim(IDN,k,j,is-i) = rho_new;
        rho_old            = prim(IDN,k,j,is-i);
        prim(IVX,k,j,is-i) = prim(IVX,k,j,is);
        prim(IVY,k,j,is-i) = 0.0;
        prim(IVZ,k,j,is-i) = 0.0;
          
        pmb->pscalars->r(t_index,k,j,is-i)    = tgs; // Just using input T at inner boundary
        pmb->pscalars->r(ye_index,k,j,is-i)   = pmb->pscalars->r(ye_index,k,j,is);
        pmb->pscalars->r(abar_index,k,j,is-i) = 1.;

        Real r_scalar[nscalar_size];
        Real s_scalar[nscalar_size];
        for (int ns=0; ns<nscalar_size; ns++) {
          r_scalar[ns] = pmb->pscalars->r(ns,k,j,is-i);
          s_scalar[ns] = pmb->pscalars->s(ns,k,j,is-i);
        }

        if (NON_BAROTROPIC_EOS) {
          prim(IPR,k,j,is-i) = pmb->peos->PresFromRhoT(prim(IDN,k,j,is-i), pmb->pscalars->r(t_index,k,j,is-i), r_scalar); //This one works for tau
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
        }
      }
    }

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(is-i)) = b.x2f(k,j,is);
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
        }
      }
    }
  }
}

// Inflow Boundary Condition at outer boundary
void InflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh) {
  AthenaArray<Real>& DATA = pmb->pmy_mesh->ruser_mesh_data[8]; 
  Real &P_Ghost{DATA(0)};
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      Real theta = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) { 
        Real r = pco->x1v(ie+i);

        // Velocity
        prim(IVX,k,j,ie+i) = v_outer;
        prim(IVY,k,j,ie+i) = 0.0;
        prim(IVZ,k,j,ie+i) = 0.0;

        // Density
        prim(IDN,k,j,ie+i) = rho_outer;

        // Passive Scalar Data
        pmb->pscalars->r(ye_index,k,j,ie+i)   = ye_outer;
        pmb->pscalars->r(t_index,k,j,ie+i)    = T_outer;
        pmb->pscalars->r(abar_index,k,j,ie+i) = abar_outer;

        Real r_scalar[nscalar_size]; 
        for (int ns=0; ns<nscalar_size; ns++) {
          r_scalar[ns] = pmb->pscalars->r(ns,k,j,is-i);
        } 

        // Pressure 
        if (NON_BAROTROPIC_EOS) {
          prim(IPR,k,j,ie+i) = P_outer;
          P_Ghost = P_outer;
        } 
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
        }
      }
    }
  }
}

void MeshBlock::UserWorkInLoop() {
  AthenaArray<Real>& my_data = pmy_mesh->ruser_mesh_data[5];
  Real &local_sum{my_data(0)};
  Real t = pmy_mesh->time;
  Real fermi2 = fermi(2,0);
  Real fermi3 = fermi(3,0);
  Real fermi4 = fermi(4,0);
  // Squared energy term from Scheck 2006 eqn. D39
  Real energy_term = (std::pow(eps_nue * fermi2 / fermi3, 2.0) * fermi4 / fermi2 + 2.0 * delta_np * (eps_nue * fermi2 / fermi3) * fermi3 / fermi2 + std::pow(delta_np,2.0)) / std::pow(mcsq, 2.0);
  Real sigma_nue_n = sigma_0 * ((1.0 + 3.0 * std::pow(g_a,2)) / (4.0)) * energy_term;
  if (lid == 0) {
    local_sum = 0.0;
  }
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) { 
        Real rho     = phydro->w(IDN,k,j,i);
        Real r_ip1   = pcoord->x1f(i+1);
        Real r_ip    = pcoord->x1f(i); 
        Real ye      = pscalars->r(ye_index,k,j,i); 
        Real deltaR  = r_ip1 - r_ip; 
        local_sum   += deltaR * sigma_nue_n * Na * (1.0 - ye) * rho;
      }
    }
  } 
  if (lid == pmy_mesh->nblocal-1) {
    Real global_integral;
    MPI_Allreduce(&local_sum, &global_integral, 1,
                  MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    global_tau = global_integral;
    Vr_FinalActiveZone = phydro->w(IVX,ke,je,ie);
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        efrac         = &pscalars->r(ye_index,k,j,i);
        Real temp     = pscalars->r(t_index,k,j,i);
        Real r        = pcoord->x1v(i);
        Real x        = std::pow((1.0-(r_0*r_0)/(r*r)),0.5);
        Real kbol_MeV = 8.61733326*std::pow(10,-11);                                      // Boltzmann constant in MeV / K
        Real ye       = pscalars->r(ye_index,k,j,i);
        Real rho      = phydro->w(IDN,k,j,i);
        Real erg2MeV  = 6.24151e5; 

        Real temp_MeV = temp * kbol_MeV;
        Real t        = pmy_mesh->time;
        Real X_n      = X_nucleon(rho, temp_MeV);
        Real ExpSupp  = std::exp(-1.0 * rho / rho_0);

        // Read in P_Ghost data
        AthenaArray<Real>& DATA    = pmy_mesh->ruser_mesh_data[8];
        Real &P_Ghost{DATA(0)};

        Real r_scalar[nscalar_size];
        for (int ns=0; ns<nscalar_size; ns++) {
          r_scalar[ns] = pscalars->r(ns,k,j,i); 
        }	
        Real cs   = std::sqrt(peos->AsqFromRhoP(phydro->w(IDN,k,j,i), phydro->w(IPR,k,j,i), r_scalar));
        Real Mach = fabs(phydro->w(IVX,k,j,i)) / cs;
        Real eta_ele  = peos->EtaFromRhoT(rho, temp, r_scalar);

        // Define user output variables:
        // Temperature
        user_out_var(0,k,j,i) = temp;
        // Qdot 
        if(Mach > 1.0) { 
          user_out_var(1,k,j,i) = 0.0;
        } else {
          user_out_var(1,k,j,i) = ExpSupp * (Heating_QW_Modified(x, ye) - Cooling_QW_Modified(rho, temp_MeV, ye, eta_ele));
        } 

        // Sound Speed
        user_out_var(2,k,j,i) = cs;

        // Tau
        user_out_var(6,k,j,i) = global_tau;

        // Abar
        user_out_var(7,k,j,i)  = pscalars->r(abar_index,k,j,i);
        user_out_var(8,k,j,i)  = P_Ghost;
        user_out_var(9,k,j,i)  = eta_ele;
      }
    }
  }
}
