/**
 *  Constants used for SMC grid propagation program.  JGLi16Oct2026
 */

#ifndef CPP_SMCONSTS_H
#define CPP_SMCONSTS_H

namespace SMConsts {

// Mathematical constants
constexpr double Pie = 3.14159265358979323846;
constexpr double TPie = 2.0 * Pie;
constexpr double HPie = 0.5 * Pie;
constexpr double Rad2D = 180.0 / Pie;
constexpr double D2Rad = Pie / 180.0;

// Physical constants
constexpr double Gravty = 9.806;
constexpr double Temp0C = 273.16;
constexpr double RDryAr = 287.05;
constexpr double REarth = 6.371E6; 
constexpr double AngulS = 7.2921E-5;
constexpr double Ang36H = 4.8481E-5;
constexpr double CLight = 2.99792458E8;

/** @brief Density of water (kg/m^3).  */
constexpr double DWAT = 1000.0;
/** @brief Density of air (kg/m^3).  */
constexpr double DAIR = 1.225;
/** @brief Kinematic viscosity of air (m^2/s). */
constexpr double NU_AIR = 1.4e-5;
/** @brief Kinematic viscosity of water (m^2/s).  */
constexpr double NU_WATER = 1.31e-6;
/** @brief Specific gravity of sediments (N.D.). */
constexpr double SED_SG = 2.65;
/** @brief von Karman's constant (N.D.). */
constexpr double KAPPA = 0.40;

// Derived physical constants
/** @brief Inverse of gravity * 2 * Pi. */
constexpr double G1PI1I = 1.0 / (Gravty * TPie);

// Model constants
/** @brief Undefined value. */
constexpr double UNDEF = -999.9;
/** @brief Minimum value for spectral density (log10). */
constexpr double ABMIN = -1.0;
/** @brief Maximum value for spectral density (log10). */
constexpr double ABMAX = 8.0;
/** @brief Maximum value for k*d in dispersion calculations. */
constexpr double KDMAX = 20.0;
/** @brief Physics factor for JONSWAP spectrum (g^2 / (2*PI)^4). */
constexpr double JONSWAP_FACTOR = 0.06175;
/** @brief Minimum water depth for refraction and group speed. */
constexpr double WtrDepthMn = 10.0;
/** @brief Maximum advection CFL number. */
constexpr double CFL_Advctn = 0.625;
/** @brief Maximum diffusion Fourier number. */
constexpr double Fourier_Df = 0.25;
/** @brief Maximum spectral rotation angle in deg. */
constexpr double Rotatn_DMx = 36.0;

} // namespace SMConsts 

#endif // CPP_SMCONSTS_H

