/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka
 *
 * RAPTOR uses cgs units for light transport calculations.
 * Entries marked [BL] are only applicable to Boyer-Lindquist coordinates.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

typedef double real;


#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <png.h>
#include <omp.h>
#include <time.h>
#include "constants.h"

// OpenACC or OMP
////////////////////

#define ACC (0)
#define OMP (1)
#define COMP OMP


#include <math.h>

//#include <accel.h>
//#include <accelmath.h>
// OUTPUT SWITCHES
////////////////////

#define VTKFILE (0)
#define IMGFILE (1)
#define SPECFILE (1)

// Radiative Transfer Switch
#define RAD_TRANS (1)

//GEOD FILE, print geodesics
#define GEOD (0)

// CAMERA SWITCHES
////////////////////
// Maximum number of pixels to compute
#define maxsize 1000000

#define LOG_IMPACT_CAM (0)
#define LINEAR_IMPACT_CAM (1)

// VR parameter ??
real Ucam[4];
real Xcam[4];
real tcam;

#define sign(x) (((x) < 0) ? -1 : ((x) > 0))

// GLOBAL VARIABLES
////////////////////
real L_unit;
real T_unit;
real RHO_unit;
real U_unit;
real B_unit;
real Ne_unit;
real Thetae_unit;

// Number of dimensions, 1 time, 3 space
#define NDIM  (4)
// Defines number of primitive variables: KHRO, UU, U1, U2, U3, B1, B2, B3 [8 in total]
#define NPRIM    (8)
// What is N1?? Resolution in r-direction
#define N1 (256)
// What is N2?? Resolution in theta-direction
#define N2 (256)
// What is N3?? Resolution in phi-direction
#define N3 (1)

// What are these variables? What are they used for?
real R0, Rin, Rh, Rout, Rms;
// What is a? Read from HARM data!
real a;
// What is hslope? Read from HARM data!
real hslope;
// What is startx, stopx and dx? Also Read from HARM data!
real startx[NDIM], stopx[NDIM], dx[NDIM];
// These variables are not used anywhere? What are they for?
real dlE, lE0;
// Read from HARM data! Used to calculate Theate Unit!
real gam;
// These variables are not used anywhere? What are they for?
real dMsim;
// Captures plasma properties from HARM data,
// KHRO(Density), UU(Internal Energy),
// U1(Radial Velocity), U2(Theta Velocity), U3(Phi Velocity),
// B1(Radial Magnetic Field), B2(Theta Magnetic Field), B3(Phi Magnetic Field)
real ****p;


// METRIC PARAMETERS
////////////////////

// These are used for geodesic integration; BH mass normalized

//coordinate and metric choices
#define CAR      (0)        //  Minkowski
#define BL       (1)        // Boyer-Lindquist,               x1=r, x2=th, and x3=phi
#define MBL      (2)        // modified Boyer-Lindquist, x1=log(r), x2=th, and x3=phi
#define KS       (3)        // Kerr-Schild,                   x1=r, x2=th, and x3=phi
#define MKS      (4)        // Proper MKS coords

// Not used anywhere
#define TT      (0)
#define RR      (1)
#define TH      (2)
#define PH      (3)

extern real a;
extern real R0;

// Metric
#define metric   (MKS)
#if(metric == CAR || metric == BL || metric == KS)
#define logscale (0)    // Standard BL/KS coordinates; no logarithmic radius
#elif(metric == MBL || metric == MKS)
#define logscale (1)    // Modified BL/KS coordinates; logarithmic radius
#endif

// MODEL PARAMETERS
///////////////////

// These are used for light transport computations; BH now has a specific mass

// GRMHD data file

// Path to the GRMHD data file (dump file from HARM)
char GRMHD_FILE[256];
// This is not used anywhere? What is it for?
char OUTPUT_FILE[256];
// This is used in a pragma call in the end of this file?
// What is it for? Not used anywhere!!
int SPHERICAL_ACC;
// This is not used anywhere? What is it for?
char TEMP_MODEL[100];
// Read from model.in! What is it for? Not used
int ABSORPTION;

// Toggle light transport calculation on/off for integration debugging
// Not used anywhere
#define LIGHT_TRANSPORT     (1)
//	#define RT_OUTER_CUTOFF     (2000.)//1.01*rcam) // Outer boundary of radiative transfer computation
// Outer boundary of radiative transfer computation,
// rcam is the camera distance from the singularity,
// define in parameters.h, LINE 194
#define RT_OUTER_CUTOFF     (1.01*rcam)
// Black hole mass
double MBH;

//sets all units, calculated in set_units() function in raptor_harm_model.c
real M_UNIT;

real R_HIGH;
real R_LOW;
#define source_dist    (2.6228263e22) // Distance to Sgr A* (cm)
//#define source_dist    (5.061e25) // Distance to M87 (cm)

// These are for simple analytical model - move to different .c file!! WHAT??
// These are NEVER used in the current version of the code!
#define n_e0           (4.5e6)   // Electron density normalization
#define B0             (100.)    // Magnetic field scalar
#define THETA_e_0      (80.)     // Dimensionless temperature scalar
#define nblobs         (1)       // number of blobs?

// Dark matter spike parameters, used in Modified Schwarzschild metric with Dark Matter present
#define qfactor      (0.1) // M_DM/M_BH, see Lacroix & Silk 2013

// OBSERVER PARAMETERS
//////////////////////

// This is not used anywhere? What is it for?
real CAM_FREQ;
// Read from Command Line Input! What is it for?
real TIME_INIT;
// Read from Command Line Input! What is it for?
real INCLINATION;

// SED parameters, Read from model.in, need to fill in the meanings here.
int  FREQS_PER_DEC;
real FREQ_MIN;
real FREQ_MAX;


//	#define rcam         (15.)//(500.)    // Camera distance from the sing.(units of Rg)
#define rcam         (1.e4)
// Image Height and Width, taken from model.in
int IMG_WIDTH;
int IMG_HEIGHT;
// Camera size in X and Y, taken from model.in
real CAM_SIZE_X;
real CAM_SIZE_Y;
#define max_order    (100)       // Maximimum order of lensed image computed (0 = direct only)

// INTEGRATOR PARAMETERS
////////////////////////
// Read from model.in, used as a param for variable step size integrator
real STEPSIZE;

#define delta_num    (1.e-7)    // Used for numerical derivatives
#define max_steps    (1e5)   // Maximum number of integration steps

#define cutoff_outer (rcam*1.01)    // Outer cutoff, near flat spacetime, in M
#define horizon_marg (1e-1)     // Stop tracing at this distance from E.H. [BL]
#define VER       (1)        //
#define RK4       (2)        //
#define RK2      (3)
#define int_method   (RK2)     // method of integration 2=Verlet, 4-RK4

// MACROS
/////////

#define DIM 4
#define LOOP_i    for(i = 0; i < DIM; i++)
#define LOOP_ij   for(i = 0; i < DIM; i++) \
for(j = 0; j < DIM; j++)

#define LOOP_kl    for(k = 0; k < DIM; k++) \
for(l = 0; l < DIM; l++)


#define LOOP_ijk2   for( i=NDIM; i--; )\
for( j=NDIM; j--; )\
for( k=NDIM; k--; )
#define LOOP_ijk    for(i = 0; i < DIM; i++) \
for(j = 0; j < DIM; j++) \
for(k = 0; k < DIM; k++)
#define LOOP_ijkl for(i = 0; i < DIM; i++) \
for(j = 0; j < DIM; j++) \
for(k = 0; k < DIM; k++) \
for(l = 0; l < DIM; l++)

// Number of frequencies for raytracing
// #define num_indices 1
// int num_indices = (int) (log(FREQ_MAX/FREQ_MIN) * FREQS_PER_DEC + 1)
extern int num_indices;

// These variables are not used anywhere
#define N_kappa 11
#define N_theta 11
#define N_theta_e 21
#define N_nuratio 56

#pragma acc copyin(IMG_WIDTH,IMG_HEIGHT,CAM_SIZE_X,CAM_SIZE_Y,cutoff_inner,ABSORPTION,MBH,INCLINATION,STEPSIZE,SPHERICAL_ACC,a,R_LOW,R_HIGH,hslope)

#endif // PARAMETERS_H
