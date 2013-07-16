#include "C_Xdot.h"

void xDotMex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    
extern double GrabnerDiskContact (   double *r20,
                              double *r10,
                              double *v20,
                              double *v10,
                              double *R20,
                              double *R10,
                              double *w20,
                              double *w10,
                              double *leParams,
                              double *aOut);



extern void C_Xdot(double t, 
            double* vX_in, 
            double* vParam, 
            double* vInput, 
            double* XDOT);