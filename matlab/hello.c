#include "mex.h" /* Always include this */
void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{
    mexPrintf("Hello, world!\n\n"); /* Do something interesting */
    return;
}
