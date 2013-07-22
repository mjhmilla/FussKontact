/*
================================================================================
 xDotMex : A wrapper for exported C code from DynaflexPro, a Maple Package
 
 @author M.Millard
 @date 2013/07/16
 
================================================================================
 */
#include <math.h>
#include "mex.h"
//#include "xdotmex.h"

#define IS_REAL_2D_DOUBLE_VECTOR(P)( !mxIsComplex(P) && \
                                     mxGetNumberOfDimensions(P) == 2 && \
                                     mxIsDouble(P) && \
                                     mxGetN(P) == 1 )
#define IS_REAL_SCALAR(P) (IS_REAL_2D_DOUBLE_VECTOR(P) && \
                           mxGetNumberOfElements(P) == 1)

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

void C_Xdot(double t, 
            double* vX_in, 
            double* vParam, 
            double* vInput, 
            double* XDOT);                              
                           
                     
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Macros for the ouput and input arguments */
    #define xDot_OUT plhs[0]        //Output argument 1    
    #define t_IN prhs[0]            //Input argument 1 (time)
    #define vXin_IN prhs[1]         //               2 (state vector)
    #define vParam_IN prhs[2]       //               3 (model parameters)
    #define vInput_IN prhs[3]       //               4 (model inputs)
    
    double *xDot, *t, *vXin, *vParam, *vInput;
    int m;
    int mContact = 12;

    //============================================================================
    //Check number of input arguments
    if(nrhs != 4 )      
        mexErrMsgTxt("Wrong number of input arguments; "
                     "4 required: t,vX_in, vParam, vInput");
    else if(nlhs > 1)   //Check number of output arguments
        mexErrMsgTxt("Too many output arguments; one allowed: XDOT");

        
   
   //============================================================================    
   //Check the type and dimensions each input argument   
    if(!IS_REAL_SCALAR(t_IN))
        mexErrMsgTxt("time must be a real double scalar.");            
    if(!IS_REAL_2D_DOUBLE_VECTOR(vXin_IN) || mxGetNumberOfElements(vXin_IN) != 14 )
        mexErrMsgTxt("vXin must be a real vector with 14 elements.");
    if(!IS_REAL_2D_DOUBLE_VECTOR(vParam_IN) || mxGetNumberOfElements(vParam_IN) != 50) 
        mexErrMsgTxt("vParam must be a real vector with 50 elements.");
    if(!IS_REAL_2D_DOUBLE_VECTOR(vInput_IN) || mxGetNumberOfElements(vInput_IN) != 7) 
        mexErrMsgTxt("vInput must be a real vector with 7 elements.");        

   //============================================================================    
   //Get the pointers to each input data field
    
    t       = mxGetPr(t_IN);
    vXin    = mxGetPr(vXin_IN);
    m       = mxGetM(vXin_IN);
    vParam  = mxGetPr(vParam_IN);
    vInput  = mxGetPr(vInput_IN);

   //============================================================================    
   //Create the output data array
    xDot_OUT    = mxCreateDoubleMatrix(m+mContact,1,mxREAL);
    xDot        = mxGetPr(xDot_OUT);
            
   //============================================================================    
   //Call xdot
    C_Xdot((double)t[0],vXin,vParam,vInput,xDot);        
            
    
    
    return;
}    


 

double GrabnerDiskContact (
  double *r20,
  double *r10,
  double *v20,
  double *v10,
  double *R20,
  double *R10,
  double *w20,
  double *w10,
  double *leParams,
  double *aOut)
{
  double sx;
  double sy;
  double sz;
  double nx;
  double ny;
  double nz;
  double t1x;
  double t1y;
  double t1z;
  double t2x;
  double t2y;
  double t2z;
  double r;
  double a;
  double k;
  double p;
  double c;
  double mus;
  double mud;
  double stVel;
  double dyVel;
  double alpha;
  double radius;
  double tmpf1;
  double tmpf2;
  double epsRoot;
  double gN;
  int i;
  //double x;
  double velNormal;
  double mu;
  double forceNormal;
  double velTangentMag;
  double velEps;
  double forceTangentMag;
  //double daFactor;
  double delta;
  double rO[3];
  double rM[3];
  double rMO[3];
  double rPM[3];
  double rPMuv[3];
  double rPMO[3];
  double rQO[3];
  double rQM[3];
  double nC[3];
  double nP[3];
  double tP1[3];
  double tP2[3];
  double vQ10[3];
  //double vQ20[3];
  //double vP10[3];
  double vP20[3];
  double vPQ[3];
  double velTangent[3];
  double forceTangentUV[3];
  double tmp1[3];
  double tmp2[3];
  for (i = 1; i <= 15; i++) //edited
    aOut[i - 1] = 0;
  epsRoot = 0.1490116119384766e-7;
  sx = leParams[0];
  sy = leParams[1];
  sz = leParams[2];
  nx = leParams[3];
  ny = leParams[4];
  nz = leParams[5];
  t1x = leParams[6];
  t1y = leParams[7];
  t1z = leParams[8];
  t2x = leParams[9];
  t2y = leParams[10];
  t2z = leParams[11];
  r = leParams[12];
  a = leParams[13];
  k = leParams[14];
  p = leParams[15];
  c = leParams[16];
  mus = leParams[17];
  mud = leParams[18];
  stVel = leParams[19];
  dyVel = leParams[20];
  nP[0] = R10[0] * nx + R10[1] * ny + R10[2] * nz;
  nP[1] = R10[3] * nx + R10[4] * ny + R10[5] * nz;
  nP[2] = R10[6] * nx + R10[7] * ny + R10[8] * nz;
  tP1[0] = R10[0] * t1x + R10[1] * t1y + R10[2] * t1z;
  tP1[1] = R10[3] * t1x + R10[4] * t1y + R10[5] * t1z;
  tP1[2] = R10[6] * t1x + R10[7] * t1y + R10[8] * t1z;
  tP2[0] = R10[0] * t2x + R10[1] * t2y + R10[2] * t2z;
  tP2[1] = R10[3] * t2x + R10[4] * t2y + R10[5] * t2z;
  tP2[2] = R10[6] * t2x + R10[7] * t2y + R10[8] * t2z;
  rO[0] = R10[0] * sx + R10[1] * sy + R10[2] * sz + r10[0];
  rO[1] = R10[3] * sx + R10[4] * sy + R10[5] * sz + r10[1];
  rO[2] = R10[6] * sx + R10[7] * sy + R10[8] * sz + r10[2];
  rM[0] = r20[0];
  rM[1] = r20[1];
  rM[2] = r20[2];
  rMO[0] = rM[0] - rO[0];
  rMO[1] = rM[1] - rO[1];
  rMO[2] = rM[2] - rO[2];
  nC[0] = R20[6];
  nC[1] = R20[7];
  nC[2] = R20[8];
  tmp1[0] = nC[1] * nP[2] - nC[2] * nP[1];
  tmp1[1] = -nC[0] * nP[2] + nC[2] * nP[0];
  tmp1[2] = nC[0] * nP[1] - nC[1] * nP[0];
  tmp2[0] = nC[1] * tmp1[2] - nC[2] * tmp1[1];
  tmp2[1] = -nC[0] * tmp1[2] + nC[2] * tmp1[0];
  tmp2[2] = nC[0] * tmp1[1] - nC[1] * tmp1[0];
  tmpf1 = nP[0] * tmp2[0] + nP[1] * tmp2[1] + nP[2] * tmp2[2];
  if (0.0e0 < tmpf1)
  {
    tmp2[0] = -0.10e1 * tmp2[0];
    tmp2[1] = -0.10e1 * tmp2[1];
    tmp2[2] = -0.10e1 * tmp2[2];
  }
  tmpf1 = pow(pow(tmp1[0], 0.2e1) + pow(tmp1[1], 0.2e1) + pow(tmp1[2], 0.2e1), 0.5e0);
  if (tmpf1 < epsRoot){
    tmpf1 = epsRoot;
  }
  tmpf2 = fabs(nC[0] * nP[0] + nC[1] * nP[1] + nC[2] * nP[2]);
  if(tmpf2 > 1.0)
    tmpf2 = 1.0;
  
  alpha = acos(tmpf2);
  radius = r * (0.1e1 - exp(- (a * sin(alpha))) );
  rPMuv[0] = tmp2[0] / tmpf1;
  rPMuv[1] = tmp2[1] / tmpf1;
  rPMuv[2] = tmp2[2] / tmpf1;    
  rPM[0] = radius * rPMuv[0];
  rPM[1] = radius * rPMuv[1];
  rPM[2] = radius * rPMuv[2];
  rPMO[0] = rMO[0] + rPM[0];
  rPMO[1] = rMO[1] + rPM[1];
  rPMO[2] = rMO[2] + rPM[2];
  gN = rPMO[0] * nP[0] + rPMO[1] * nP[1] + rPMO[2] * nP[2];
  if (gN < 0.0e0)
  {
    vP20[0] = -rPM[1] * w20[2] + rPM[2] * w20[1] + v20[0];
    vP20[1] = rPM[0] * w20[2] - rPM[2] * w20[0] + v20[1];
    vP20[2] = -rPM[0] * w20[1] + rPM[1] * w20[0] + v20[2];
    tmpf1 = rPMO[0] * tP1[0] + rPMO[1] * tP1[1] + rPMO[2] * tP1[2];
    tmpf2 = rPMO[0] * tP2[0] + rPMO[1] * tP2[1] + rPMO[2] * tP2[2];
    rQO[0] = tmpf1 * tP1[0] + tmpf2 * tP2[0];
    rQO[1] = tmpf1 * tP1[1] + tmpf2 * tP2[1];
    rQO[2] = tmpf1 * tP1[2] + tmpf2 * tP2[2];
    vQ10[0] = -rQO[1] * w10[2] + rQO[2] * w10[1] + v10[0];
    vQ10[1] = rQO[0] * w10[2] - rQO[2] * w10[0] + v10[1];
    vQ10[2] = -rQO[0] * w10[1] + rQO[1] * w10[0] + v10[2];
    vPQ[0] = vP20[0] - vQ10[0];
    vPQ[1] = vP20[1] - vQ10[1];
    vPQ[2] = vP20[2] - vQ10[2];
    velNormal = nP[0] * vPQ[0] + nP[1] * vPQ[1] + nP[2] * vPQ[2];
    velTangent[0] = vPQ[0] * (tP1[0] + tP2[0]);
    velTangent[1] = vPQ[1] * (tP1[1] + tP2[1]);
    velTangent[2] = vPQ[2] * (tP1[2] + tP2[2]);
    velTangentMag = pow(pow(velTangent[0], 0.2e1) + pow(velTangent[1], 0.2e1) + pow(velTangent[2], 0.2e1), 0.5e0);
    mu = 0.0e0;
    delta = 0.0e0;
    if (0.0e0 < velTangentMag)
      if (dyVel < velTangentMag)
        mu = mud;
      else if (stVel < velTangentMag && velTangentMag < dyVel)
      {
        delta = (velTangentMag - stVel) / (dyVel - stVel);
        mu = mus + (mud - mus) * delta * delta * (0.30e1 - 0.20e1 * delta);
      }
      else
      {
        delta = (velTangentMag + stVel) * 0.5000000000e0 / stVel;
        mu = -mus + 0.2e1 * mus * delta * delta * (0.30e1 - 0.20e1 * delta);
      }
    forceNormal = k * pow(fabs(gN), p) * (0.1e1 - c * velNormal);
    if (forceNormal < 0.0e0)
      forceNormal = 0.0e0;
    forceTangentMag = -mu * forceNormal;
    velEps = 0.10000e4 * epsRoot;
    if (velEps < velTangentMag)
    {
      forceTangentUV[0] = velTangent[0] / velTangentMag;
      forceTangentUV[1] = velTangent[1] / velTangentMag;
      forceTangentUV[2] = velTangent[2] / velTangentMag;
    }
    else
    {
      tmpf1 = velTangent[0] / velEps;
      forceTangentUV[0] = tmpf1 * (0.1500000000e1 * fabs(tmpf1) - 0.5000000000e0 * fabs(pow(tmpf1, 0.3e1)));
      tmpf1 = velTangent[1] / velEps;
      forceTangentUV[1] = tmpf1 * (0.1500000000e1 * fabs(tmpf1) - 0.5000000000e0 * fabs(pow(tmpf1, 0.3e1)));
      tmpf1 = velTangent[2] / velEps;
      forceTangentUV[2] = tmpf1 * (0.1500000000e1 * fabs(tmpf1) - 0.5000000000e0 * fabs(pow(tmpf1, 0.3e1)));
    }
    rQM[0] = rQO[0] - rMO[0];
    rQM[1] = rQO[1] - rMO[1];
    rQM[2] = rQO[2] - rMO[2];
    aOut[0] = (forceNormal * nP[0] + forceTangentMag * forceTangentUV[0]);
    aOut[1] = (forceNormal * nP[1] + forceTangentMag * forceTangentUV[1]);
    aOut[2] = (forceNormal * nP[2] + forceTangentMag * forceTangentUV[2]);
    aOut[3] = (rPM[1] * (double) aOut[2] - rPM[2] * (double) aOut[1]);
    aOut[4] = (-rPM[0] * (double) aOut[2] + rPM[2] * (double) aOut[0]);
    aOut[5] = (rPM[0] * (double) aOut[1] - rPM[1] * (double) aOut[0]);
    aOut[6] = -aOut[0];
    aOut[7] = -aOut[1];
    aOut[8] = -aOut[2];
    aOut[9] =  (rPMO[1] * (double) aOut[8] - rPMO[2] * (double) aOut[7]);
    aOut[10] = (-rPMO[0] * (double) aOut[8] + rPMO[2] * (double) aOut[6]);
    aOut[11] = (rPMO[0] * (double) aOut[7] - rPMO[1] * (double) aOut[6]);
    

    
  }
  
      //edited - only available in this copy
    aOut[12] = rPMO[0] + rO[0];
    aOut[13] = rPMO[1] + rO[1];
    aOut[14] = rPMO[2] + rO[2];
  
  return(0.0e0);
}

void C_Xdot(double t, double* vX_in, double* vParam, double* vInput, double* XDOT)
{
  double zz[227];
double CFHeel[15]; //edited
double R1Heel[9];
double R2Heel[9];
double r1Heel[3];
double r2Heel[3];
double v1Heel[3];
double v2Heel[3];
double w1Heel[3];
double w2Heel[3];
double CFForefoot[15]; //edited
double CParamsForefoot[21];
double CParamsHeel[21];
double R1Forefoot[9];
double R2Forefoot[9];
double r1Forefoot[3];
double r2Forefoot[3];
double v1Forefoot[3];
double v2Forefoot[3];
double w1Forefoot[3];
double w2Forefoot[3];
 
CParamsForefoot[0] = vParam[41];
CParamsForefoot[1] = vParam[42];
CParamsForefoot[2] = vParam[43];
CParamsForefoot[3] = vParam[33];
CParamsForefoot[4] = vParam[34];
CParamsForefoot[5] = vParam[35];
CParamsForefoot[6] = vParam[44];
CParamsForefoot[7] = vParam[45];
CParamsForefoot[8] = vParam[46];
CParamsForefoot[9] = vParam[47];
CParamsForefoot[10] = vParam[48];
CParamsForefoot[11] = vParam[49];
CParamsForefoot[12] = vParam[38];
CParamsForefoot[13] = vParam[22];
CParamsForefoot[14] = vParam[27];
CParamsForefoot[15] = vParam[36];
CParamsForefoot[16] = vParam[24];
CParamsForefoot[17] = vParam[32];
CParamsForefoot[18] = vParam[31];
CParamsForefoot[19] = vParam[40];
CParamsForefoot[20] = vParam[26];
zz[0] = cos(vX_in[13]);
zz[24] = cos(vX_in[12]);
zz[35] = cos(vParam[15]);
zz[46] = sin(vX_in[13]);
zz[57] = sin(vParam[15]);
zz[68] = zz[24] * (zz[0] * zz[35] - 0.1e1 * zz[46] * zz[57]);
R1Forefoot[0] = zz[68];
zz[79] = sin(vX_in[12]);
zz[90] = sin(vX_in[11]);
zz[101] = cos(vX_in[11]);
zz[1] = zz[0] * zz[90];
zz[12] = zz[46] * zz[101];
zz[16] = zz[1] * zz[79] + zz[12];
zz[17] = zz[46] * zz[90];
zz[18] = zz[0] * zz[101];
zz[19] = -0.1e1 * zz[17] * zz[79] + zz[18];
zz[20] = zz[16] * zz[35] + zz[19] * zz[57];
R1Forefoot[1] = zz[20];
zz[17] = -0.1e1 * zz[18] * zz[79] + zz[17];
zz[1] = zz[12] * zz[79] + zz[1];
zz[12] = zz[1] * zz[57] + zz[35] * zz[17];
R1Forefoot[2] = zz[12];
zz[18] = cos(vX_in[10]);
zz[21] = sin(vX_in[10]);
zz[22] = zz[0] * zz[57] + zz[35] * zz[46];
zz[23] = zz[24] * zz[18];
zz[25] = zz[21] * zz[79] - 0.1e1 * zz[23] * zz[22];
R1Forefoot[3] = zz[25];
zz[26] = -0.1e1 * zz[16] * zz[57] + zz[19] * zz[35];
zz[27] = zz[90] * zz[24];
zz[28] = zz[18] * zz[26] - 0.1e1 * zz[27] * zz[21];
R1Forefoot[4] = zz[28];
zz[29] = zz[1] * zz[35] - 0.1e1 * zz[17] * zz[57];
zz[30] = zz[101] * zz[24];
zz[31] = zz[18] * zz[29] + zz[30] * zz[21];
R1Forefoot[5] = zz[31];
zz[32] = zz[21] * zz[22] * zz[24] + zz[18] * zz[79];
R1Forefoot[6] = zz[32];
zz[33] = -0.1e1 * zz[21] * zz[26] - 0.1e1 * zz[23] * zz[90];
R1Forefoot[7] = zz[33];
zz[23] = -0.1e1 * zz[21] * zz[29] + zz[23] * zz[101];
R1Forefoot[8] = zz[23];
zz[34] = vX_in[3] * zz[35] + vX_in[4];
zz[36] = vX_in[3] * zz[57] + vX_in[5];
zz[37] = zz[79] * vX_in[6];
w1Forefoot[0] = zz[24] * (zz[0] * zz[34] - 0.1e1 * zz[36] * zz[46]) + zz[37];
zz[38] = zz[27] * vX_in[6];
w1Forefoot[1] = zz[16] * zz[34] + zz[19] * zz[36] - 0.1e1 * zz[38];
zz[39] = zz[30] * vX_in[6];
w1Forefoot[2] = zz[1] * zz[36] + zz[17] * zz[34] + zz[39];
R2Forefoot[0] = 0.1e1;
R2Forefoot[1] = 0.0e0;
R2Forefoot[2] = 0.0e0;
R2Forefoot[3] = 0.0e0;
R2Forefoot[4] = 0.1e1;
R2Forefoot[5] = 0.0e0;
R2Forefoot[6] = 0.0e0;
R2Forefoot[7] = 0.0e0;
R2Forefoot[8] = 0.1e1;
w2Forefoot[0] = 0.0e0;
w2Forefoot[1] = 0.0e0;
w2Forefoot[2] = 0.0e0;
zz[40] = -0.1e1 * vParam[9] + vParam[18];
zz[41] = -0.1e1 * vParam[7] + vParam[16];
zz[42] = vParam[8] - 0.1e1 * vParam[17];
zz[43] = zz[68] * vParam[19];
zz[44] = zz[25] * vParam[20];
zz[45] = zz[32] * vParam[21];
r1Forefoot[0] = zz[24] * (zz[0] * zz[41] + zz[42] * zz[46]) + zz[40] * zz[79] - 0.1e1 * zz[43] - 0.1e1 * zz[44] - 0.1e1 * zz[45] + vX_in[7];
zz[47] = zz[20] * vParam[19];
zz[48] = zz[28] * vParam[20];
zz[49] = zz[33] * vParam[21];
r1Forefoot[1] = zz[16] * zz[41] - 0.1e1 * zz[19] * zz[42] - 0.1e1 * zz[27] * zz[40] - 0.1e1 * zz[47] - 0.1e1 * zz[48] - 0.1e1 * zz[49] + vX_in[8];
zz[50] = zz[12] * vParam[19];
zz[51] = zz[31] * vParam[20];
zz[52] = zz[23] * vParam[21];
r1Forefoot[2] = -0.1e1 * zz[1] * zz[42] + zz[17] * zz[41] + zz[30] * zz[40] - 0.1e1 * zz[50] - 0.1e1 * zz[51] - 0.1e1 * zz[52] + vX_in[9];
zz[40] = vParam[8] * vX_in[6] - 0.1e1 * vParam[9] * vX_in[5];
zz[41] = -0.1e1 * vParam[7] * vX_in[6] + vParam[9] * vX_in[4];
zz[42] = vParam[7] * vX_in[5] - 0.1e1 * vParam[8] * vX_in[4];
zz[53] = -0.1e1 * vParam[17] * vX_in[6] + vParam[18] * vX_in[5];
zz[54] = vParam[16] * vX_in[6] - 0.1e1 * vParam[18] * vX_in[4];
zz[55] = -0.1e1 * vParam[16] * vX_in[5] + vParam[17] * vX_in[4];
zz[56] = vParam[20] * zz[21] + vParam[21] * zz[18];
zz[18] = zz[18] * vParam[20] - 0.1e1 * zz[21] * vParam[21];
zz[21] = vParam[19] * zz[57] + zz[18] * zz[35];
zz[58] = zz[21] * vX_in[6] - 0.1e1 * zz[36] * zz[56];
zz[18] = -0.1e1 * zz[18] * zz[57] + zz[35] * vParam[19];
zz[56] = -0.1e1 * zz[18] * vX_in[6] + zz[34] * zz[56];
zz[18] = zz[36] * zz[18] - 0.1e1 * zz[21] * zz[34];
zz[21] = zz[42] + zz[55] + zz[18];
zz[59] = zz[40] + zz[53] + zz[58];
zz[60] = -0.1e1 * zz[41] - 0.1e1 * zz[54] - 0.1e1 * zz[56];
v1Forefoot[0] = zz[21] * zz[79] + zz[24] * (zz[0] * zz[59] + zz[46] * zz[60]) + vX_in[0];
v1Forefoot[1] = zz[16] * zz[59] - 0.1e1 * zz[19] * zz[60] - 0.1e1 * zz[27] * zz[21] + vX_in[1];
v1Forefoot[2] = -0.1e1 * zz[1] * zz[60] + zz[17] * zz[59] + zz[30] * zz[21] + vX_in[2];
r2Forefoot[0] = 0.0e0;
r2Forefoot[1] = 0.0e0;
r2Forefoot[2] = 0.0e0;
v2Forefoot[0] = 0.0e0;
v2Forefoot[1] = 0.0e0;
v2Forefoot[2] = 0.0e0;
zz[112] = GrabnerDiskContact(r1Forefoot, r2Forefoot, v1Forefoot, v2Forefoot, R1Forefoot, R2Forefoot, w1Forefoot, w2Forefoot, CParamsForefoot, CFForefoot);
zz[114] = CFForefoot[0];
zz[115] = CFForefoot[1];
zz[116] = CFForefoot[2];
zz[117] = CFForefoot[3];
zz[118] = CFForefoot[4];
zz[119] = CFForefoot[5];
CParamsHeel[0] = vParam[41];
CParamsHeel[1] = vParam[42];
CParamsHeel[2] = vParam[43];
CParamsHeel[3] = vParam[33];
CParamsHeel[4] = vParam[34];
CParamsHeel[5] = vParam[35];
CParamsHeel[6] = vParam[44];
CParamsHeel[7] = vParam[45];
CParamsHeel[8] = vParam[46];
CParamsHeel[9] = vParam[47];
CParamsHeel[10] = vParam[48];
CParamsHeel[11] = vParam[49];
CParamsHeel[12] = vParam[39];
CParamsHeel[13] = vParam[23];
CParamsHeel[14] = vParam[28];
CParamsHeel[15] = vParam[37];
CParamsHeel[16] = vParam[25];
CParamsHeel[17] = vParam[32];
CParamsHeel[18] = vParam[31];
CParamsHeel[19] = vParam[40];
CParamsHeel[20] = vParam[26];
zz[21] = cos(vParam[11]);
zz[59] = sin(vParam[11]);
zz[60] = sin(vParam[10]);
zz[61] = cos(vParam[10]);
zz[62] = zz[46] * zz[24];
zz[63] = -0.1e1 * zz[62] * zz[60] - 0.1e1 * zz[61] * zz[79];
zz[64] = zz[24] * zz[0];
R1Heel[0] = zz[21] * zz[64] + zz[59] * zz[63];
zz[65] = zz[19] * zz[60] + zz[27] * zz[61];
R1Heel[1] = zz[16] * zz[21] + zz[59] * zz[65];
zz[66] = zz[1] * zz[60] - 0.1e1 * zz[30] * zz[61];
R1Heel[2] = zz[17] * zz[21] + zz[59] * zz[66];
R1Heel[3] = zz[60] * zz[79] - 0.1e1 * zz[62] * zz[61];
R1Heel[4] = zz[19] * zz[61] - 0.1e1 * zz[60] * zz[27];
R1Heel[5] = zz[1] * zz[61] + zz[30] * zz[60];
R1Heel[6] = -0.1e1 * zz[21] * zz[63] + zz[64] * zz[59];
R1Heel[7] = zz[16] * zz[59] - 0.1e1 * zz[21] * zz[65];
R1Heel[8] = zz[17] * zz[59] - 0.1e1 * zz[21] * zz[66];
zz[21] = zz[0] * vX_in[4] - 0.1e1 * zz[46] * vX_in[5];
zz[37] = zz[21] * zz[24] + zz[37];
w1Heel[0] = zz[37];
zz[38] = zz[16] * vX_in[4] + zz[19] * vX_in[5] - 0.1e1 * zz[38];
w1Heel[1] = zz[38];
zz[39] = vX_in[4] * zz[17] + vX_in[5] * zz[1] + zz[39];
w1Heel[2] = zz[39];
R2Heel[0] = 0.1e1;
R2Heel[1] = 0.0e0;
R2Heel[2] = 0.0e0;
R2Heel[3] = 0.0e0;
R2Heel[4] = 0.1e1;
R2Heel[5] = 0.0e0;
R2Heel[6] = 0.0e0;
R2Heel[7] = 0.0e0;
R2Heel[8] = 0.1e1;
w2Heel[0] = w2Forefoot[0];
w2Heel[1] = w2Forefoot[1];
w2Heel[2] = w2Forefoot[2];
zz[59] = -0.1e1 * vParam[9] + vParam[14];
zz[60] = -0.1e1 * vParam[7] + vParam[12];
zz[61] = vParam[8] - 0.1e1 * vParam[13];
r1Heel[0] = zz[24] * (zz[0] * zz[60] + zz[46] * zz[61]) + zz[59] * zz[79] + vX_in[7];
r1Heel[1] = zz[16] * zz[60] - 0.1e1 * zz[19] * zz[61] - 0.1e1 * zz[27] * zz[59] + vX_in[8];
r1Heel[2] = -0.1e1 * zz[1] * zz[61] + zz[17] * zz[60] + zz[30] * zz[59] + vX_in[9];
zz[59] = -0.1e1 * vParam[13] * vX_in[6] + vParam[14] * vX_in[5] + zz[40];
zz[60] = -0.1e1 * vParam[12] * vX_in[5] + vParam[13] * vX_in[4] + zz[42];
zz[61] = -0.1e1 * vParam[12] * vX_in[6] + vParam[14] * vX_in[4] - 0.1e1 * zz[41];
v1Heel[0] = zz[24] * (zz[0] * zz[59] + zz[46] * zz[61]) + zz[60] * zz[79] + vX_in[0];
v1Heel[1] = zz[16] * zz[59] - 0.1e1 * zz[19] * zz[61] - 0.1e1 * zz[60] * zz[27] + vX_in[1];
v1Heel[2] = -0.1e1 * zz[1] * zz[61] + zz[17] * zz[59] + zz[30] * zz[60] + vX_in[2];
r2Heel[0] = r2Forefoot[0];
r2Heel[1] = r2Forefoot[1];
r2Heel[2] = r2Forefoot[2];
v2Heel[0] = v2Forefoot[0];
v2Heel[1] = v2Forefoot[1];
v2Heel[2] = v2Forefoot[2];
zz[113] = GrabnerDiskContact(r1Heel, r2Heel, v1Heel, v2Heel, R1Heel, R2Heel, w1Heel, w2Heel, CParamsHeel, CFHeel);
zz[120] = CFHeel[0];
zz[121] = CFHeel[1];
zz[122] = CFHeel[2];
zz[123] = CFHeel[3];
zz[124] = CFHeel[4];
zz[125] = CFHeel[5];
zz[59] = vParam[29] + vParam[30];
zz[170] = zz[59];
zz[50] = zz[50] + zz[51] + zz[52];
zz[47] = zz[47] + zz[48] + zz[49];
zz[48] = zz[12] * zz[47] - 0.1e1 * zz[20] * zz[50];
zz[49] = vParam[30] * zz[48];
zz[171] = zz[49];
zz[51] = vParam[7] * zz[17] + vParam[8] * zz[1] + vParam[9] * zz[30];
zz[52] = zz[16] * vParam[7] + zz[19] * vParam[8] - 0.1e1 * zz[27] * vParam[9];
zz[60] = -0.1e1 * zz[16] * zz[51] + zz[17] * zz[52];
zz[61] = vParam[16] * zz[17] + vParam[17] * zz[1] + vParam[18] * zz[30];
zz[63] = zz[16] * vParam[16] + zz[19] * vParam[17] - 0.1e1 * zz[27] * vParam[18];
zz[65] = zz[47] + zz[52] - 0.1e1 * zz[63];
zz[66] = -0.1e1 * zz[50] - 0.1e1 * zz[51] + zz[61];
zz[67] = zz[16] * zz[66] + zz[17] * zz[65];
zz[69] = vParam[29] * zz[60] + vParam[30] * zz[67];
zz[172] = zz[69];
zz[70] = zz[1] * zz[52] - 0.1e1 * zz[19] * zz[51];
zz[71] = zz[1] * zz[65] + zz[19] * zz[66];
zz[72] = vParam[29] * zz[70] + vParam[30] * zz[71];
zz[173] = zz[72];
zz[73] = zz[90] * zz[51];
zz[74] = zz[101] * zz[52];
zz[75] = zz[24] * (zz[90] * (zz[50] + zz[51] - 0.1e1 * zz[61]) + zz[101] * (zz[47] + zz[52] - 0.1e1 * zz[63]));
zz[76] = vParam[29] * zz[24] * (zz[73] + zz[74]) + vParam[30] * zz[75];
zz[174] = zz[76];
zz[175] = zz[59];
zz[43] = zz[43] + zz[44] + zz[45];
zz[44] = -0.1e1 * zz[12] * zz[43] + zz[50] * zz[68];
zz[45] = vParam[30] * zz[44];
zz[176] = zz[45];
zz[77] = zz[24] * (zz[0] * vParam[7] - 0.1e1 * zz[46] * vParam[8]) + zz[79] * vParam[9];
zz[78] = zz[24] * (zz[0] * vParam[16] - 0.1e1 * zz[46] * vParam[17]) + zz[79] * vParam[18];
zz[80] = -0.1e1 * zz[43] - 0.1e1 * zz[77] + zz[78];
zz[81] = zz[17] * zz[80] - 0.1e1 * zz[64] * zz[66];
zz[82] = vParam[29] * (-0.1e1 * zz[17] * zz[77] + zz[51] * zz[64]) + vParam[30] * zz[81];
zz[177] = zz[82];
zz[83] = zz[1] * zz[80] + zz[62] * zz[66];
zz[84] = vParam[29] * (-0.1e1 * zz[1] * zz[77] - 0.1e1 * zz[62] * zz[51]) + vParam[30] * zz[83];
zz[178] = zz[84];
zz[66] = zz[30] * zz[80] - 0.1e1 * zz[66] * zz[79];
zz[85] = vParam[29] * (-0.1e1 * zz[30] * zz[77] + zz[51] * zz[79]) + vParam[30] * zz[66];
zz[179] = zz[85];
zz[180] = zz[59];
zz[86] = zz[20] * zz[43] - 0.1e1 * zz[47] * zz[68];
zz[87] = vParam[30] * zz[86];
zz[181] = zz[87];
zz[88] = -0.1e1 * zz[16] * zz[80] - 0.1e1 * zz[64] * zz[65];
zz[89] = vParam[29] * (zz[16] * zz[77] - 0.1e1 * zz[64] * zz[52]) + vParam[30] * zz[88];
zz[182] = zz[89];
zz[91] = -0.1e1 * zz[19] * zz[80] + zz[62] * zz[65];
zz[92] = vParam[29] * (zz[19] * zz[77] + zz[52] * zz[62]) + vParam[30] * zz[91];
zz[183] = zz[92];
zz[65] = zz[27] * zz[80] - 0.1e1 * zz[65] * zz[79];
zz[80] = vParam[29] * (-0.1e1 * zz[27] * zz[77] - 0.1e1 * zz[52] * zz[79]) + vParam[30] * zz[65];
zz[184] = zz[80];
zz[185] = zz[49];
zz[186] = zz[45];
zz[187] = zz[87];
zz[45] = vParam[2] * pow(zz[68], 2) + vParam[4] * pow(zz[25], 2) + vParam[6] * pow(zz[32], 2);
zz[32] = zz[32] * vParam[6];
zz[25] = zz[25] * vParam[4];
zz[49] = zz[68] * vParam[2];
zz[87] = zz[20] * zz[49] + zz[25] * zz[28] + zz[32] * zz[33];
zz[25] = zz[12] * zz[49] + zz[23] * zz[32] + zz[25] * zz[31];
zz[32] = zz[12] * zz[25] + zz[20] * zz[87] + zz[45] * zz[68];
zz[49] = vParam[2] * pow(zz[20], 2) + vParam[4] * pow(zz[28], 2) + vParam[6] * pow(zz[33], 2);
zz[28] = vParam[2] * zz[12] * zz[20] + vParam[4] * zz[28] * zz[31] + vParam[6] * zz[23] * zz[33];
zz[33] = zz[12] * zz[28] + zz[20] * zz[49] + zz[68] * zz[87];
zz[23] = vParam[2] * pow(zz[12], 2) + vParam[4] * pow(zz[31], 2) + vParam[6] * pow(zz[23], 2);
zz[31] = zz[12] * zz[23] + zz[20] * zz[28] + zz[25] * zz[68];
zz[93] = vParam[30] * (zz[44] * zz[50] - 0.1e1 * zz[47] * zz[86]);
zz[94] = vParam[30] * (zz[43] * zz[86] - 0.1e1 * zz[48] * zz[50]);
zz[95] = vParam[30] * (-0.1e1 * zz[43] * zz[44] + zz[47] * zz[48]);
zz[188] = zz[12] * (zz[31] + zz[95]) + zz[20] * (zz[33] + zz[94]) + zz[68] * (zz[32] + zz[93]);
zz[96] = zz[16] * zz[87] + zz[17] * zz[25] + zz[45] * zz[64];
zz[97] = zz[16] * zz[49] + zz[17] * zz[28] + zz[64] * zz[87];
zz[98] = zz[16] * zz[28] + zz[17] * zz[23] + zz[25] * zz[64];
zz[99] = vParam[30] * (-0.1e1 * zz[47] * zz[88] + zz[50] * zz[81]);
zz[100] = vParam[30] * (zz[43] * zz[88] - 0.1e1 * zz[50] * zz[67]);
zz[102] = vParam[30] * (-0.1e1 * zz[43] * zz[81] + zz[47] * zz[67]);
zz[189] = zz[12] * (zz[98] + zz[102]) + zz[20] * (zz[97] + zz[100]) + zz[68] * (zz[96] + zz[99]);
zz[103] = zz[1] * zz[25] + zz[19] * zz[87] - 0.1e1 * zz[62] * zz[45];
zz[104] = zz[1] * zz[28] + zz[19] * zz[49] - 0.1e1 * zz[62] * zz[87];
zz[105] = zz[1] * zz[23] + zz[19] * zz[28] - 0.1e1 * zz[62] * zz[25];
zz[106] = vParam[30] * (-0.1e1 * zz[47] * zz[91] + zz[50] * zz[83]);
zz[107] = vParam[30] * (zz[43] * zz[91] - 0.1e1 * zz[50] * zz[71]);
zz[108] = vParam[30] * (-0.1e1 * zz[43] * zz[83] + zz[47] * zz[71]);
zz[190] = zz[12] * (zz[105] + zz[108]) + zz[20] * (zz[104] + zz[107]) + zz[68] * (zz[103] + zz[106]);
zz[109] = zz[24] * (zz[25] * zz[101] - 0.1e1 * zz[87] * zz[90]) + zz[45] * zz[79];
zz[110] = zz[24] * (zz[28] * zz[101] - 0.1e1 * zz[49] * zz[90]) + zz[79] * zz[87];
zz[111] = zz[24] * (zz[23] * zz[101] - 0.1e1 * zz[28] * zz[90]) + zz[25] * zz[79];
zz[2] = vParam[30] * (-0.1e1 * zz[47] * zz[65] + zz[50] * zz[66]);
zz[3] = vParam[30] * (zz[43] * zz[65] - 0.1e1 * zz[50] * zz[75]);
zz[4] = vParam[30] * (-0.1e1 * zz[43] * zz[66] + zz[47] * zz[75]);
zz[191] = zz[12] * (zz[111] + zz[4]) + zz[20] * (zz[110] + zz[3]) + zz[68] * (zz[109] + zz[2]);
zz[59] = -0.1e1 * zz[59];
zz[5] = -0.1e1 * zz[50] + zz[61];
zz[6] = zz[47] - 0.1e1 * zz[63];
zz[192] = -0.1e1 * zz[60] * zz[59] + (zz[5] * zz[16] + zz[6] * zz[17]) * vParam[30];
zz[60] = (-0.1e1 * zz[43] + zz[78]) * vParam[30];
zz[7] = zz[77] * zz[59];
zz[8] = zz[60] + zz[7];
zz[9] = zz[5] * vParam[30];
zz[10] = zz[51] * zz[59];
zz[11] = -0.1e1 * zz[9] - 0.1e1 * zz[10];
zz[193] = zz[8] * zz[17] + zz[11] * zz[64];
zz[13] = zz[6] * vParam[30];
zz[14] = zz[52] * zz[59];
zz[15] = -0.1e1 * zz[13] + zz[14];
zz[60] = -0.1e1 * zz[60] - 0.1e1 * zz[7];
zz[194] = zz[15] * zz[64] + zz[16] * zz[60];
zz[31] = zz[31] + zz[95] + vParam[30] * (-0.1e1 * zz[44] * zz[77] + zz[48] * zz[52]) + vParam[30] * (zz[44] * zz[78] - 0.1e1 * zz[48] * zz[63]);
zz[32] = zz[32] + zz[93] + vParam[30] * (zz[44] * zz[51] - 0.1e1 * zz[52] * zz[86]) + vParam[30] * (-0.1e1 * zz[61] * zz[44] + zz[63] * zz[86]);
zz[33] = zz[33] + zz[94] + vParam[30] * (-0.1e1 * zz[48] * zz[51] + zz[77] * zz[86]) + vParam[30] * (zz[48] * zz[61] - 0.1e1 * zz[78] * zz[86]);
zz[195] = zz[16] * zz[33] + zz[17] * zz[31] + zz[32] * zz[64];
zz[44] = pow(zz[24], 2);
zz[48] = pow(zz[79], 2) * vParam[5] + zz[44] * (vParam[1] * pow(zz[0], 2) + vParam[3] * pow(zz[46], 2));
zz[86] = (zz[0] * zz[16] * vParam[1] - 0.1e1 * zz[46] * zz[19] * vParam[3] - 0.1e1 * zz[79] * zz[90] * vParam[5]) * zz[24];
zz[93] = (zz[0] * zz[17] * vParam[1] - 0.1e1 * zz[46] * zz[1] * vParam[3] + zz[79] * zz[101] * vParam[5]) * zz[24];
zz[94] = vParam[5] * zz[44] * pow(zz[90], 2) + vParam[1] * pow(zz[16], 2) + vParam[3] * pow(zz[19], 2);
zz[95] = -0.1e1 * zz[44] * zz[90] * vParam[5] * zz[101] + zz[19] * vParam[3] * zz[1] + zz[16] * vParam[1] * zz[17];
zz[44] = vParam[5] * zz[44] * pow(zz[101], 2) + vParam[1] * pow(zz[17], 2) + vParam[3] * pow(zz[1], 2);
zz[97] = zz[16] * zz[94] + zz[17] * zz[95] - 0.1e1 * zz[51] * zz[69] + zz[64] * zz[86] + zz[77] * zz[89] + zz[97] + zz[100] + vParam[30] * (zz[61] * zz[67] - 0.1e1 * zz[78] * zz[88]);
zz[67] = zz[16] * zz[95] + zz[17] * zz[44] + zz[52] * zz[69] + zz[64] * zz[93] - 0.1e1 * zz[77] * zz[82] + zz[98] + zz[102] + vParam[30] * (-0.1e1 * zz[63] * zz[67] + zz[78] * zz[81]);
zz[69] = zz[16] * zz[86] + zz[17] * zz[93] + zz[64] * zz[48] + zz[51] * zz[82] - 0.1e1 * zz[52] * zz[89] + zz[96] + zz[99] + vParam[30] * (-0.1e1 * zz[61] * zz[81] + zz[63] * zz[88]);
zz[196] = zz[16] * zz[97] + zz[17] * zz[67] + zz[64] * zz[69];
zz[81] = zz[1] * zz[93] + zz[19] * zz[86] - 0.1e1 * zz[62] * zz[48] + zz[51] * zz[84] - 0.1e1 * zz[52] * zz[92] + zz[103] + zz[106] + vParam[30] * (-0.1e1 * zz[61] * zz[83] + zz[63] * zz[91]);
zz[82] = zz[1] * zz[95] + zz[19] * zz[94] - 0.1e1 * zz[51] * zz[72] - 0.1e1 * zz[62] * zz[86] + zz[77] * zz[92] + zz[104] + zz[107] + vParam[30] * (zz[61] * zz[71] - 0.1e1 * zz[78] * zz[91]);
zz[71] = zz[1] * zz[44] + zz[19] * zz[95] + zz[52] * zz[72] - 0.1e1 * zz[62] * zz[93] - 0.1e1 * zz[77] * zz[84] + zz[105] + zz[108] + vParam[30] * (-0.1e1 * zz[63] * zz[71] + zz[78] * zz[83]);
zz[197] = zz[16] * zz[82] + zz[17] * zz[71] + zz[64] * zz[81];
zz[72] = zz[110] + zz[3] + zz[24] * (-0.1e1 * zz[90] * zz[94] + zz[95] * zz[101]) + zz[79] * zz[86] - 0.1e1 * zz[51] * zz[76] + zz[77] * zz[80] + vParam[30] * (zz[61] * zz[75] - 0.1e1 * zz[65] * zz[78]);
zz[75] = zz[111] + zz[4] + zz[24] * (zz[44] * zz[101] - 0.1e1 * zz[90] * zz[95]) + zz[79] * zz[93] + zz[52] * zz[76] - 0.1e1 * zz[77] * zz[85] + vParam[30] * (-0.1e1 * zz[63] * zz[75] + zz[66] * zz[78]);
zz[65] = zz[109] + zz[2] + zz[24] * (-0.1e1 * zz[86] * zz[90] + zz[93] * zz[101]) + zz[48] * zz[79] + zz[51] * zz[85] - 0.1e1 * zz[52] * zz[80] + vParam[30] * (-0.1e1 * zz[61] * zz[66] + zz[63] * zz[65]);
zz[198] = zz[16] * zz[72] + zz[17] * zz[75] + zz[64] * zz[65];
zz[199] = -0.1e1 * zz[70] * zz[59] + (zz[1] * zz[6] + zz[5] * zz[19]) * vParam[30];
zz[200] = zz[1] * zz[8] + zz[62] * (zz[9] + zz[10]);
zz[201] = zz[19] * zz[60] + zz[62] * (zz[13] - 0.1e1 * zz[14]);
zz[202] = zz[1] * zz[31] + zz[19] * zz[33] - 0.1e1 * zz[62] * zz[32];
zz[203] = zz[1] * zz[67] + zz[19] * zz[97] - 0.1e1 * zz[62] * zz[69];
zz[204] = zz[1] * zz[71] + zz[19] * zz[82] - 0.1e1 * zz[62] * zz[81];
zz[205] = zz[1] * zz[75] + zz[19] * zz[72] - 0.1e1 * zz[62] * zz[65];
zz[206] = zz[24] * ((-0.1e1 * zz[74] - 0.1e1 * zz[73]) * zz[59] + (zz[90] * (zz[50] - 0.1e1 * zz[61]) + zz[101] * (zz[47] - 0.1e1 * zz[63])) * vParam[30]);
zz[207] = zz[8] * zz[30] + zz[11] * zz[79];
zz[208] = zz[8] * zz[27] + zz[15] * zz[79];
zz[209] = zz[24] * (zz[31] * zz[101] - 0.1e1 * zz[33] * zz[90]) + zz[32] * zz[79];
zz[210] = zz[24] * (zz[67] * zz[101] - 0.1e1 * zz[90] * zz[97]) + zz[69] * zz[79];
zz[211] = zz[24] * (zz[71] * zz[101] - 0.1e1 * zz[82] * zz[90]) + zz[79] * zz[81];
zz[212] = zz[24] * (-0.1e1 * zz[72] * zz[90] + zz[75] * zz[101]) + zz[65] * zz[79];
zz[31] = w2Heel[0] * zz[39];
zz[32] = w2Heel[2] * zz[37];
zz[33] = -0.1e1 * zz[31] + zz[32];
zz[59] = w2Heel[0] * zz[38];
zz[60] = w2Heel[1] * zz[37];
zz[65] = zz[59] - 0.1e1 * zz[60];
zz[66] = vX_in[6] * zz[41] - 0.1e1 * vX_in[5] * zz[42];
zz[42] = -0.1e1 * vX_in[6] * zz[40] + vX_in[4] * zz[42];
zz[40] = vX_in[5] * zz[40] - 0.1e1 * vX_in[4] * zz[41];
zz[35] = vX_in[3] * (-0.1e1 * zz[35] * vX_in[5] + zz[57] * vX_in[4]);
zz[57] = vX_in[6] * vX_in[3];
zz[26] = zz[57] * zz[26] - 0.1e1 * zz[27] * zz[35] - 0.1e1 * zz[31] + zz[32];
zz[29] = zz[57] * zz[29] + zz[30] * zz[35] + zz[59] - 0.1e1 * zz[60];
zz[31] = zz[36] * zz[18] - 0.1e1 * zz[54] * vX_in[6] + zz[55] * vX_in[5] - 0.1e1 * vX_in[6] * zz[56] - 0.1e1 * zz[66];
zz[32] = zz[34] * zz[56] - 0.1e1 * zz[36] * zz[58] - 0.1e1 * zz[53] * vX_in[5] + zz[54] * vX_in[4] - 0.1e1 * zz[40];
zz[18] = zz[18] * zz[34] - 0.1e1 * zz[53] * vX_in[6] + zz[55] * vX_in[4] - 0.1e1 * vX_in[6] * zz[58] + zz[42];
zz[34] = zz[52] - 0.1e1 * zz[63];
zz[36] = -0.1e1 * zz[51] + zz[61];
zz[41] = vParam[29] * (zz[24] * (-0.1e1 * zz[0] * zz[66] + zz[42] * zz[46]) - 0.1e1 * zz[33] * zz[51] - 0.1e1 * zz[40] * zz[79] + zz[52] * zz[65]);
zz[53] = vParam[30] * (zz[24] * (zz[0] * zz[31] + zz[18] * zz[46]) - 0.1e1 * zz[26] * zz[50] + zz[29] * zz[47] + zz[32] * zz[79] + zz[33] * zz[36] + zz[34] * zz[65]);
zz[218] = zz[114] + zz[120] + vInput[0] - 0.1e1 * zz[41] - 0.1e1 * zz[53];
zz[54] = w2Heel[1] * zz[39];
zz[55] = w2Heel[2] * zz[38];
zz[56] = zz[54] - 0.1e1 * zz[55];
zz[35] = -0.1e1 * zz[57] * zz[24] * zz[22] + zz[35] * zz[79] + zz[54] - 0.1e1 * zz[55];
zz[57] = -0.1e1 * zz[77] + zz[78];
zz[22] = vParam[29] * (-0.1e1 * zz[16] * zz[66] - 0.1e1 * zz[19] * zz[42] + zz[27] * zz[40] + zz[51] * zz[56] - 0.1e1 * zz[65] * zz[77]);
zz[36] = vParam[30] * (zz[16] * zz[31] - 0.1e1 * zz[18] * zz[19] - 0.1e1 * zz[27] * zz[32] - 0.1e1 * zz[29] * zz[43] + zz[35] * zz[50] - 0.1e1 * zz[36] * zz[56] + zz[57] * zz[65]);
zz[219] = zz[115] + zz[121] + vInput[1] - 0.1e1 * zz[22] - 0.1e1 * zz[36];
zz[40] = (zz[1] * zz[42] + zz[17] * zz[66] + zz[30] * zz[40] - 0.1e1 * zz[33] * zz[77] + zz[52] * zz[56] - 0.1e1 * vParam[0]) * vParam[29];
zz[57] = (zz[1] * zz[18] - 0.1e1 * zz[17] * zz[31] - 0.1e1 * zz[26] * zz[43] - 0.1e1 * zz[30] * zz[32] + zz[33] * zz[57] + zz[34] * zz[56] + zz[35] * zz[47] - 0.1e1 * vParam[0]) * vParam[30];
zz[220] = zz[116] + zz[122] + vInput[2] + zz[40] + zz[57];
zz[18] = zz[23] * w1Forefoot[2] + zz[25] * w1Forefoot[0] + zz[28] * w1Forefoot[1];
zz[31] = zz[28] * w1Forefoot[2] + zz[49] * w1Forefoot[1] + zz[87] * w1Forefoot[0];
zz[32] = zz[18] * w1Forefoot[1] + zz[25] * zz[29] + zz[26] * zz[87] - 0.1e1 * zz[31] * w1Forefoot[2] + zz[35] * zz[45];
zz[34] = zz[25] * w1Forefoot[2] + zz[45] * w1Forefoot[0] + zz[87] * w1Forefoot[1];
zz[18] = -0.1e1 * zz[18] * w1Forefoot[0] + zz[26] * zz[49] + zz[28] * zz[29] + zz[34] * w1Forefoot[2] + zz[35] * zz[87];
zz[35] = zz[23] * zz[29] + zz[25] * zz[35] + zz[26] * zz[28] + zz[31] * w1Forefoot[0] - 0.1e1 * zz[34] * w1Forefoot[1];
zz[23] = zz[116] + zz[57];
zz[25] = zz[115] - 0.1e1 * zz[36];
zz[26] = zz[23] * zz[47] - 0.1e1 * zz[25] * zz[50];
zz[28] = zz[114] - 0.1e1 * zz[53];
zz[29] = -0.1e1 * zz[23] * zz[43] + zz[28] * zz[50];
zz[31] = zz[25] * zz[43] - 0.1e1 * zz[28] * zz[47];
zz[221] = zz[12] * (zz[119] - 0.1e1 * zz[35] - 0.1e1 * zz[31]) + zz[20] * (zz[118] - 0.1e1 * zz[18] - 0.1e1 * zz[29]) + zz[68] * (zz[117] - 0.1e1 * zz[32] - 0.1e1 * zz[26]) + vInput[6];
zz[68] = zz[115] + zz[121] - 0.1e1 * zz[22] - 0.1e1 * zz[36];
zz[12] = zz[114] + zz[120] - 0.1e1 * zz[41] - 0.1e1 * zz[53];
zz[20] = -0.1e1 * zz[12] * zz[52] + zz[68] * zz[77];
zz[57] = zz[116] + zz[122] + zz[40] + zz[57];
zz[12] = zz[12] * zz[51] - 0.1e1 * zz[57] * zz[77];
zz[22] = -0.1e1 * zz[78] * zz[25] + zz[63] * zz[28];
zz[28] = zz[78] * zz[23] - 0.1e1 * zz[61] * zz[28];
zz[34] = zz[37] * zz[93] + zz[38] * zz[95] + zz[39] * zz[44];
zz[36] = zz[37] * zz[48] + zz[38] * zz[86] + zz[39] * zz[93];
zz[40] = zz[33] * zz[94] - 0.1e1 * zz[34] * zz[37] + zz[36] * zz[39] + zz[56] * zz[86] + zz[65] * zz[95];
zz[41] = zz[37] * zz[86] + zz[38] * zz[94] + zz[39] * zz[95];
zz[27] = zz[16] * vParam[12] + zz[19] * vParam[13] - 0.1e1 * zz[27] * vParam[14];
zz[30] = vParam[12] * zz[17] + vParam[13] * zz[1] + vParam[14] * zz[30];
zz[42] = zz[24] * (zz[0] * vParam[12] - 0.1e1 * zz[46] * vParam[13]) + zz[79] * vParam[14];
zz[43] = zz[27] * zz[120] - 0.1e1 * zz[42] * zz[121];
zz[42] = -0.1e1 * zz[30] * zz[120] + zz[42] * zz[122];
zz[36] = zz[33] * zz[95] - 0.1e1 * zz[36] * zz[38] + zz[37] * zz[41] + zz[44] * zz[65] + zz[56] * zz[93];
zz[57] = zz[63] * zz[23] - 0.1e1 * zz[61] * zz[25] + zz[27] * zz[122] - 0.1e1 * zz[30] * zz[121] - 0.1e1 * zz[33] * zz[86] - 0.1e1 * zz[34] * zz[38] + zz[39] * zz[41] - 0.1e1 * zz[48] * zz[56] + zz[51] * zz[68] - 0.1e1 * zz[52] * zz[57] - 0.1e1 * zz[65] * zz[93] - 0.1e1 * zz[26] - 0.1e1 * zz[32] + zz[117] + zz[123] + vInput[3];
zz[222] = zz[16] * (zz[118] - 0.1e1 * zz[18] - 0.1e1 * zz[12]) + zz[17] * (zz[119] + zz[125] - 0.1e1 * zz[35] - 0.1e1 * zz[31] - 0.1e1 * zz[20] - 0.1e1 * zz[22]) + zz[16] * (zz[124] - 0.1e1 * zz[29] - 0.1e1 * zz[28] - 0.1e1 * zz[40] - 0.1e1 * zz[42] + vInput[4]) + zz[17] * (-0.1e1 * zz[43] - 0.1e1 * zz[36] + vInput[5]) + zz[64] * zz[57];
zz[223] = zz[1] * (-0.1e1 * zz[35] - 0.1e1 * zz[31] - 0.1e1 * zz[20] - 0.1e1 * zz[22] - 0.1e1 * zz[36]) + zz[19] * (-0.1e1 * zz[18] - 0.1e1 * zz[29] - 0.1e1 * zz[12] - 0.1e1 * zz[28]) + zz[1] * (zz[119] + zz[125] - 0.1e1 * zz[43] + vInput[5]) + zz[19] * (zz[118] + zz[124] - 0.1e1 * zz[40] - 0.1e1 * zz[42] + vInput[4]) - 0.1e1 * zz[62] * zz[57];
zz[224] = zz[57] * zz[79] + zz[24] * (zz[90] * (-0.1e1 * zz[118] - 0.1e1 * zz[124] + zz[18] + zz[29] + zz[12] + zz[28] + zz[40] + zz[42] - 0.1e1 * vInput[4]) + zz[101] * (zz[119] + zz[125] - 0.1e1 * zz[35] - 0.1e1 * zz[31] - 0.1e1 * zz[20] - 0.1e1 * zz[22] - 0.1e1 * zz[43] - 0.1e1 * zz[36] + vInput[5]));
zz[225] = vX_in[0];
zz[226] = vX_in[1];
zz[213] = vX_in[2];
zz[214] = vX_in[3];
zz[35] = 0.1e1 / zz[24];
zz[215] = zz[21] * zz[35];
zz[216] = vX_in[4] * zz[46] + vX_in[5] * zz[0];
zz[217] = -0.1e1 * (zz[21] * zz[79] - 0.1e1 * zz[24] * vX_in[6]) * zz[35];
zz[0] = 0.1e1 / zz[170];
zz[137] = zz[171] * zz[0];
zz[138] = zz[172] * zz[0];
zz[139] = zz[173] * zz[0];
zz[140] = zz[174] * zz[0];
zz[149] = -0.1e1 * zz[185] * zz[137] + zz[188];
zz[150] = -0.1e1 * zz[185] * zz[138] + zz[189];
zz[151] = -0.1e1 * zz[185] * zz[139] + zz[190];
zz[152] = -0.1e1 * zz[185] * zz[140] + zz[191];
zz[153] = -0.1e1 * zz[192] * zz[137] + zz[195];
zz[154] = -0.1e1 * zz[192] * zz[138] + zz[196];
zz[155] = -0.1e1 * zz[192] * zz[139] + zz[197];
zz[156] = -0.1e1 * zz[192] * zz[140] + zz[198];
zz[157] = -0.1e1 * zz[199] * zz[137] + zz[202];
zz[158] = -0.1e1 * zz[199] * zz[138] + zz[203];
zz[159] = -0.1e1 * zz[199] * zz[139] + zz[204];
zz[160] = -0.1e1 * zz[199] * zz[140] + zz[205];
zz[161] = -0.1e1 * zz[206] * zz[137] + zz[209];
zz[162] = -0.1e1 * zz[206] * zz[138] + zz[210];
zz[163] = -0.1e1 * zz[206] * zz[139] + zz[211];
zz[164] = -0.1e1 * zz[206] * zz[140] + zz[212];
zz[24] = 0.1e1 / zz[175];
zz[141] = zz[176] * zz[24];
zz[142] = zz[177] * zz[24];
zz[143] = zz[178] * zz[24];
zz[144] = zz[179] * zz[24];
zz[149] = -0.1e1 * zz[186] * zz[141] + zz[149];
zz[150] = -0.1e1 * zz[186] * zz[142] + zz[150];
zz[151] = -0.1e1 * zz[186] * zz[143] + zz[151];
zz[152] = -0.1e1 * zz[186] * zz[144] + zz[152];
zz[153] = -0.1e1 * zz[193] * zz[141] + zz[153];
zz[154] = -0.1e1 * zz[193] * zz[142] + zz[154];
zz[155] = -0.1e1 * zz[193] * zz[143] + zz[155];
zz[156] = -0.1e1 * zz[193] * zz[144] + zz[156];
zz[157] = -0.1e1 * zz[200] * zz[141] + zz[157];
zz[158] = -0.1e1 * zz[200] * zz[142] + zz[158];
zz[159] = -0.1e1 * zz[200] * zz[143] + zz[159];
zz[160] = -0.1e1 * zz[200] * zz[144] + zz[160];
zz[161] = -0.1e1 * zz[207] * zz[141] + zz[161];
zz[162] = -0.1e1 * zz[207] * zz[142] + zz[162];
zz[163] = -0.1e1 * zz[207] * zz[143] + zz[163];
zz[164] = -0.1e1 * zz[207] * zz[144] + zz[164];
zz[35] = 0.1e1 / zz[180];
zz[145] = zz[181] * zz[35];
zz[146] = zz[182] * zz[35];
zz[147] = zz[183] * zz[35];
zz[148] = zz[184] * zz[35];
zz[149] = -0.1e1 * zz[187] * zz[145] + zz[149];
zz[150] = -0.1e1 * zz[187] * zz[146] + zz[150];
zz[151] = -0.1e1 * zz[187] * zz[147] + zz[151];
zz[152] = -0.1e1 * zz[187] * zz[148] + zz[152];
zz[153] = -0.1e1 * zz[194] * zz[145] + zz[153];
zz[154] = -0.1e1 * zz[194] * zz[146] + zz[154];
zz[155] = -0.1e1 * zz[194] * zz[147] + zz[155];
zz[156] = -0.1e1 * zz[194] * zz[148] + zz[156];
zz[157] = -0.1e1 * zz[201] * zz[145] + zz[157];
zz[158] = -0.1e1 * zz[201] * zz[146] + zz[158];
zz[159] = -0.1e1 * zz[201] * zz[147] + zz[159];
zz[160] = -0.1e1 * zz[201] * zz[148] + zz[160];
zz[161] = -0.1e1 * zz[208] * zz[145] + zz[161];
zz[162] = -0.1e1 * zz[208] * zz[146] + zz[162];
zz[163] = -0.1e1 * zz[208] * zz[147] + zz[163];
zz[164] = -0.1e1 * zz[208] * zz[148] + zz[164];
zz[46] = 0.1e1 / zz[149];
zz[150] = zz[150] * zz[46];
zz[151] = zz[151] * zz[46];
zz[152] = zz[152] * zz[46];
zz[154] = -0.1e1 * zz[153] * zz[150] + zz[154];
zz[155] = -0.1e1 * zz[153] * zz[151] + zz[155];
zz[156] = -0.1e1 * zz[153] * zz[152] + zz[156];
zz[158] = -0.1e1 * zz[157] * zz[150] + zz[158];
zz[159] = -0.1e1 * zz[157] * zz[151] + zz[159];
zz[160] = -0.1e1 * zz[157] * zz[152] + zz[160];
zz[162] = -0.1e1 * zz[161] * zz[150] + zz[162];
zz[163] = -0.1e1 * zz[161] * zz[151] + zz[163];
zz[164] = -0.1e1 * zz[161] * zz[152] + zz[164];
zz[57] = 0.1e1 / zz[154];
zz[155] = zz[155] * zz[57];
zz[156] = zz[156] * zz[57];
zz[159] = -0.1e1 * zz[158] * zz[155] + zz[159];
zz[160] = -0.1e1 * zz[158] * zz[156] + zz[160];
zz[163] = -0.1e1 * zz[162] * zz[155] + zz[163];
zz[164] = -0.1e1 * zz[162] * zz[156] + zz[164];
zz[68] = 0.1e1 / zz[159];
zz[160] = zz[160] * zz[68];
zz[164] = -0.1e1 * zz[163] * zz[160] + zz[164];
zz[130] = zz[218] * zz[0];
zz[126] = -0.1e1 * zz[130] * zz[185] + zz[221];
zz[127] = -0.1e1 * zz[130] * zz[192] + zz[222];
zz[128] = -0.1e1 * zz[130] * zz[199] + zz[223];
zz[129] = -0.1e1 * zz[130] * zz[206] + zz[224];
zz[131] = zz[219] * zz[24];
zz[126] = -0.1e1 * zz[131] * zz[186] + zz[126];
zz[127] = -0.1e1 * zz[131] * zz[193] + zz[127];
zz[128] = -0.1e1 * zz[131] * zz[200] + zz[128];
zz[129] = -0.1e1 * zz[131] * zz[207] + zz[129];
zz[132] = zz[220] * zz[35];
zz[126] = -0.1e1 * zz[132] * zz[187] + zz[126];
zz[127] = -0.1e1 * zz[132] * zz[194] + zz[127];
zz[128] = -0.1e1 * zz[132] * zz[201] + zz[128];
zz[129] = -0.1e1 * zz[132] * zz[208] + zz[129];
zz[133] = zz[126] * zz[46];
zz[127] = -0.1e1 * zz[133] * zz[153] + zz[127];
zz[128] = -0.1e1 * zz[133] * zz[157] + zz[128];
zz[129] = -0.1e1 * zz[133] * zz[161] + zz[129];
zz[134] = zz[127] * zz[57];
zz[128] = -0.1e1 * zz[134] * zz[158] + zz[128];
zz[129] = -0.1e1 * zz[134] * zz[162] + zz[129];
zz[135] = zz[128] * zz[68];
zz[129] = -0.1e1 * zz[135] * zz[163] + zz[129];
zz[136] = zz[129] / zz[164];
zz[0] = -0.1e1 * zz[136] * zz[160] + zz[135];
zz[169] = zz[0];
zz[24] = -0.1e1 * zz[136] * zz[156] - 0.1e1 * zz[169] * zz[155] + zz[134];
zz[168] = zz[24];
zz[35] = -0.1e1 * zz[136] * zz[152] - 0.1e1 * zz[168] * zz[150] - 0.1e1 * zz[169] * zz[151] + zz[133];
zz[167] = zz[35];
zz[46] = -0.1e1 * zz[136] * zz[148] - 0.1e1 * zz[167] * zz[145] - 0.1e1 * zz[168] * zz[146] - 0.1e1 * zz[169] * zz[147] + zz[132];
zz[166] = zz[46];
zz[57] = -0.1e1 * zz[136] * zz[144] - 0.1e1 * zz[167] * zz[141] - 0.1e1 * zz[168] * zz[142] - 0.1e1 * zz[169] * zz[143] + zz[131];
zz[165] = zz[57];
XDOT[0] = -0.1e1 * zz[136] * zz[140] - 0.1e1 * zz[167] * zz[137] - 0.1e1 * zz[168] * zz[138] - 0.1e1 * zz[169] * zz[139] + zz[130];
XDOT[1] = zz[57];
XDOT[2] = zz[46];
XDOT[3] = zz[35];
XDOT[4] = zz[24];
XDOT[5] = zz[0];
XDOT[6] = zz[136];
XDOT[7] = zz[225];
XDOT[8] = zz[226];
XDOT[9] = zz[213];
XDOT[10] = zz[214];
XDOT[11] = zz[215];
XDOT[12] = zz[216];
XDOT[13] = zz[217];

//edited
XDOT[14] = CFHeel[0]; //Heel Fx
XDOT[15] = CFHeel[1]; //Heel Fy
XDOT[16] = CFHeel[2]; //Heel Fz

XDOT[17] = CFHeel[12]; //Heel COP x
XDOT[18] = CFHeel[13]; //Heel COP y 
XDOT[19] = CFHeel[14]; //Heel COP z

XDOT[20] = CFForefoot[0]; //Forefoot Fx
XDOT[21] = CFForefoot[1]; //Forefoot Fy
XDOT[22] = CFForefoot[2]; //Forefoot Fz

XDOT[23] = CFForefoot[12]; //Forefoot COP x
XDOT[24] = CFForefoot[13]; //Forefoot COP x
XDOT[25] = CFForefoot[14]; //Forefoot COP x

}


