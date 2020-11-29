#include<stdio.h>
#include<stddef.h>
#include<math.h>
#include<cuComplex.h>
#include <cub/cub.cuh>
#include </usr/local/cuda-9.2/include/curand.h>
#include </usr/local/cuda-9.2/include/curand_kernel.h>
#include </usr/local/cuda-9.2/include/cuda.h>
#include </usr/local/cuda-9.2/include/cufft.h>
#include "/usr/local/cuda-9.2/samples/common/inc/helper_cuda.h"

#define PI acos(-1.0) 
#define Tolerance 1.0e-08
#define COMPERR 1.0e-08
#define Re 0
#define Im 1 



// Variables : composition field, derivative of bulk free energy w.r.t. composition


//NEW VARIABLES 
cufftDoubleComplex *comp, *dfdc, *dfdcx, *dfdcy, *dfdcz;
cufftDoubleComplex *comp_d, *dfdc_d, *dfdcx_d, *dfdcy_d, *dfdcz_d;
//Plan required in main.cu and evolve so global declaration
cufftHandle plan;
int blocks;



//p_up for forward transform and p_dn for backward transform

//total number of simulation steps
int num_steps;

//Configuration to be initialized or to be read
int initcount, initflag;

// Alloy composition, amplitude of white noise to be added to the system
double alloycomp, noise_level; 

//Step size along x and y, timestep
double dx, dy, dz, dt;

//Total simulation time (nondimensional)
double sim_time, total_time;

//System dimensions along x and y
int nx, ny,nz, nx_half, ny_half, nz_half;

//Bulk free energy coefficients
double A;

//Gradient energy coefficients associated with composition and structural
//order parameter fields
double kappa_c;

//Mobility of solute required for CH equation (mobility)
//Relaxation coefficient for CA equation (relax_coeff)
double P;
 
//Required for scaling
double one_by_nxnynz;



FILE *fpout;


