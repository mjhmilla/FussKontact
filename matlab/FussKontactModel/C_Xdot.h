#include <math.h>
#include <stdlib.h>
 
#ifndef DFPEXTERNALS
//EXTERNAL FUNCTION DEFINITIONS
//Required External Function

    double GrabnerDiskContact (   double *r20,
                                  double *r10,
                                  double *v20,
                                  double *v10,
                                  double *R20,
                                  double *R10,
                                  double *w20,
                                  double *w10,
                                  double *leParams,
                                  double *aOut);

#define DFPEXTERNALS 
#endif
//MAIN FUNCTION DEFINITION
#ifndef DFPPROC_Xdot

    void C_Xdot(double t, 
                double* vX_in, 
                double* vParam, 
                double* vInput, 
                double* XDOT);

#define DFPPROC_Xdot
#endif
 
