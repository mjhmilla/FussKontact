#include <math.h>
#include <stdlib.h>
 
#ifndef DFPEXTERNALS
//EXTERNAL FUNCTION DEFINITIONS
//Required External Function

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
  int *aOut)
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
  double x;
  double velNormal;
  double mu;
  double forceNormal;
  double velTangentMag;
  double velEps;
  double forceTangentMag;
  double daFactor;
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
  double vQ20[3];
  double vP10[3];
  double vP20[3];
  double vPQ[3];
  double velTangent[3];
  double forceTangentUV[3];
  double tmp1[3];
  double tmp2[3];
  for (i = 1; i <= 12; i++)
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
  if (tmpf1 < epsRoot)
    tmpf1 = epsRoot;
  alpha = acos(fabs(nC[0] * nP[0] + nC[1] * nP[1] + nC[2] * nP[2]));
  radius = r * (0.1e1 - exp(-a * sin(alpha)));
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
    aOut[0] = (int) (forceNormal * nP[0] + forceTangentMag * forceTangentUV[0]);
    aOut[1] = (int) (forceNormal * nP[1] + forceTangentMag * forceTangentUV[1]);
    aOut[2] = (int) (forceNormal * nP[2] + forceTangentMag * forceTangentUV[2]);
    aOut[3] = (int) (rPM[1] * (double) aOut[2] - rPM[2] * (double) aOut[1]);
    aOut[4] = (int) (-rPM[0] * (double) aOut[2] + rPM[2] * (double) aOut[0]);
    aOut[5] = (int) (rPM[0] * (double) aOut[1] - rPM[1] * (double) aOut[0]);
    aOut[6] = -aOut[0];
    aOut[7] = -aOut[1];
    aOut[8] = -aOut[2];
    aOut[9] = (int) (rPMO[1] * (double) aOut[8] - rPMO[2] * (double) aOut[7]);
    aOut[10] = (int) (-rPMO[0] * (double) aOut[8] + rPMO[2] * (double) aOut[6]);
    aOut[11] = (int) (rPMO[0] * (double) aOut[7] - rPMO[1] * (double) aOut[6]);
  }
  return(0.0e0);
}

#define DFPEXTERNALS 
#endif
//MAIN FUNCTION DEFINITION
#ifndef DFPPROC_Xdot
void C_Xdot(const double t, const double* vX_in, const double* vParam, const double* vInput, double* XDOT)
{
  double zz[78];
double CFContact11[12];
double R1Contact11[9];
double R2Contact11[9];
double r1Contact11[3];
double r2Contact11[3];
double v1Contact11[3];
double v2Contact11[3];
double w1Contact11[3];
double w2Contact11[3];
double CParamsContact11[21];
 
CParamsContact11[0] = vParam[17];
CParamsContact11[1] = vParam[18];
CParamsContact11[2] = vParam[19];
CParamsContact11[3] = vParam[11];
CParamsContact11[4] = vParam[12];
CParamsContact11[5] = vParam[13];
CParamsContact11[6] = vParam[20];
CParamsContact11[7] = vParam[21];
CParamsContact11[8] = vParam[22];
CParamsContact11[9] = vParam[23];
CParamsContact11[10] = vParam[24];
CParamsContact11[11] = vParam[25];
CParamsContact11[12] = vParam[15];
CParamsContact11[13] = vParam[4];
CParamsContact11[14] = vParam[7];
CParamsContact11[15] = vParam[14];
CParamsContact11[16] = vParam[5];
CParamsContact11[17] = vParam[10];
CParamsContact11[18] = vParam[9];
CParamsContact11[19] = vParam[16];
CParamsContact11[20] = vParam[6];
zz[0] = cos(vX_in[11]);
zz[11] = cos(vX_in[10]);
zz[22] = zz[0] * zz[11];
R1Contact11[0] = zz[22];
zz[25] = sin(vX_in[10]);
zz[26] = sin(vX_in[9]);
zz[27] = sin(vX_in[11]);
zz[28] = cos(vX_in[9]);
zz[29] = zz[0] * zz[26];
zz[30] = zz[27] * zz[28];
zz[1] = zz[25] * zz[29] + zz[30];
R1Contact11[1] = zz[1];
zz[2] = zz[0] * zz[28];
zz[3] = zz[27] * zz[26];
zz[4] = -0.1e1 * zz[2] * zz[25] + zz[3];
R1Contact11[2] = zz[4];
zz[5] = zz[27] * zz[11];
R1Contact11[3] = -0.1e1 * zz[5];
zz[2] = -0.1e1 * zz[3] * zz[25] + zz[2];
R1Contact11[4] = zz[2];
zz[29] = zz[25] * zz[30] + zz[29];
R1Contact11[5] = zz[29];
R1Contact11[6] = zz[25];
zz[30] = zz[11] * zz[26];
R1Contact11[7] = -0.1e1 * zz[30];
zz[3] = zz[11] * zz[28];
R1Contact11[8] = zz[3];
zz[6] = zz[0] * vX_in[3] - 0.1e1 * zz[27] * vX_in[4];
zz[7] = vX_in[5] * zz[25] + zz[6] * zz[11];
w1Contact11[0] = zz[7];
zz[30] = zz[1] * vX_in[3] + zz[2] * vX_in[4] - 0.1e1 * zz[30] * vX_in[5];
w1Contact11[1] = zz[30];
zz[3] = vX_in[3] * zz[4] + vX_in[4] * zz[29] + vX_in[5] * zz[3];
w1Contact11[2] = zz[3];
R2Contact11[0] = 0.1e1;
R2Contact11[1] = 0.0e0;
R2Contact11[2] = 0.0e0;
R2Contact11[3] = 0.0e0;
R2Contact11[4] = 0.1e1;
R2Contact11[5] = 0.0e0;
R2Contact11[6] = 0.0e0;
R2Contact11[7] = 0.0e0;
R2Contact11[8] = 0.1e1;
w2Contact11[0] = 0.0e0;
w2Contact11[1] = 0.0e0;
w2Contact11[2] = 0.0e0;
r1Contact11[0] = vX_in[6];
r1Contact11[1] = vX_in[7];
r1Contact11[2] = vX_in[8];
v1Contact11[0] = vX_in[0];
v1Contact11[1] = vX_in[1];
v1Contact11[2] = vX_in[2];
r2Contact11[0] = 0.0e0;
r2Contact11[1] = 0.0e0;
r2Contact11[2] = 0.0e0;
v2Contact11[0] = 0.0e0;
v2Contact11[1] = 0.0e0;
v2Contact11[2] = 0.0e0;
zz[31] = GrabnerDiskContact(r1Contact11, r2Contact11, v1Contact11, v2Contact11, R1Contact11, R2Contact11, w1Contact11, w2Contact11, CParamsContact11, CFContact11);
zz[32] = CFContact11[0];
zz[33] = CFContact11[1];
zz[34] = CFContact11[2];
zz[35] = CFContact11[3];
zz[36] = CFContact11[4];
zz[37] = CFContact11[5];
zz[54] = vParam[8];
zz[55] = vParam[8];
zz[56] = vParam[8];
zz[8] = pow(zz[11], 2);
zz[9] = pow(zz[25], 2) * vParam[3] + zz[8] * (vParam[1] * pow(zz[0], 2) + vParam[2] * pow(zz[27], 2));
zz[10] = (zz[0] * zz[1] * vParam[1] - 0.1e1 * zz[27] * zz[2] * vParam[2] - 0.1e1 * zz[25] * zz[26] * vParam[3]) * zz[11];
zz[12] = (zz[0] * zz[4] * vParam[1] + zz[25] * zz[28] * vParam[3] - 0.1e1 * zz[27] * zz[29] * vParam[2]) * zz[11];
zz[13] = zz[1] * zz[10] + zz[4] * zz[12] + zz[9] * zz[22];
zz[14] = vParam[3] * zz[8] * pow(zz[26], 2) + vParam[1] * pow(zz[1], 2) + vParam[2] * pow(zz[2], 2);
zz[15] = -0.1e1 * zz[8] * zz[26] * vParam[3] * zz[28] + zz[1] * vParam[1] * zz[4] + zz[2] * vParam[2] * zz[29];
zz[16] = zz[1] * zz[14] + zz[4] * zz[15] + zz[10] * zz[22];
zz[8] = vParam[3] * zz[8] * pow(zz[28], 2) + vParam[1] * pow(zz[4], 2) + vParam[2] * pow(zz[29], 2);
zz[17] = zz[1] * zz[15] + zz[4] * zz[8] + zz[12] * zz[22];
zz[57] = zz[1] * zz[16] + zz[4] * zz[17] + zz[13] * zz[22];
zz[18] = zz[2] * zz[10] - 0.1e1 * zz[5] * zz[9] + zz[12] * zz[29];
zz[19] = zz[2] * zz[14] - 0.1e1 * zz[5] * zz[10] + zz[15] * zz[29];
zz[20] = zz[2] * zz[15] - 0.1e1 * zz[5] * zz[12] + zz[8] * zz[29];
zz[58] = zz[1] * zz[19] + zz[4] * zz[20] + zz[18] * zz[22];
zz[21] = zz[9] * zz[25] + zz[11] * (-0.1e1 * zz[10] * zz[26] + zz[12] * zz[28]);
zz[23] = zz[10] * zz[25] + zz[11] * (-0.1e1 * zz[14] * zz[26] + zz[15] * zz[28]);
zz[24] = zz[11] * (zz[8] * zz[28] - 0.1e1 * zz[15] * zz[26]) + zz[12] * zz[25];
zz[59] = zz[1] * zz[23] + zz[4] * zz[24] + zz[21] * zz[22];
zz[60] = zz[2] * zz[16] - 0.1e1 * zz[5] * zz[13] + zz[17] * zz[29];
zz[61] = zz[2] * zz[19] - 0.1e1 * zz[5] * zz[18] + zz[20] * zz[29];
zz[62] = zz[2] * zz[23] - 0.1e1 * zz[5] * zz[21] + zz[24] * zz[29];
zz[63] = zz[11] * (-0.1e1 * zz[16] * zz[26] + zz[17] * zz[28]) + zz[13] * zz[25];
zz[64] = zz[11] * (-0.1e1 * zz[19] * zz[26] + zz[20] * zz[28]) + zz[18] * zz[25];
zz[65] = zz[11] * (-0.1e1 * zz[23] * zz[26] + zz[24] * zz[28]) + zz[21] * zz[25];
zz[69] = zz[32];
zz[70] = zz[33];
zz[71] = -0.1e1 * vParam[0] * vParam[8] + zz[34];
zz[13] = zz[3] * w2Contact11[1] - 0.1e1 * zz[30] * w2Contact11[2];
zz[3] = -0.1e1 * zz[3] * w2Contact11[0] + zz[7] * w2Contact11[2];
zz[30] = -0.1e1 * zz[7] * w2Contact11[1] + zz[30] * w2Contact11[0];
zz[7] = zz[8] * w1Contact11[2] + zz[12] * w1Contact11[0] + zz[15] * w1Contact11[1];
zz[16] = zz[10] * w1Contact11[0] + zz[14] * w1Contact11[1] + zz[15] * w1Contact11[2];
zz[17] = zz[9] * w1Contact11[0] + zz[10] * w1Contact11[1] + zz[12] * w1Contact11[2];
zz[8] = -0.1e1 * zz[3] * zz[15] - 0.1e1 * zz[8] * zz[30] - 0.1e1 * zz[12] * zz[13] - 0.1e1 * zz[16] * w1Contact11[0] + zz[17] * w1Contact11[1] + zz[37];
zz[9] = -0.1e1 * zz[3] * zz[10] - 0.1e1 * zz[7] * w1Contact11[1] - 0.1e1 * zz[9] * zz[13] - 0.1e1 * zz[12] * zz[30] + zz[16] * w1Contact11[2] + zz[35];
zz[30] = -0.1e1 * zz[3] * zz[14] + zz[7] * w1Contact11[0] - 0.1e1 * zz[10] * zz[13] - 0.1e1 * zz[15] * zz[30] - 0.1e1 * zz[17] * w1Contact11[2] + zz[36];
zz[72] = zz[1] * zz[30] + zz[4] * zz[8] + zz[9] * zz[22];
zz[73] = zz[2] * zz[30] - 0.1e1 * zz[5] * zz[9] + zz[8] * zz[29];
zz[74] = zz[9] * zz[25] + zz[11] * (zz[8] * zz[28] - 0.1e1 * zz[26] * zz[30]);
zz[75] = vX_in[0];
zz[76] = vX_in[1];
zz[77] = vX_in[2];
zz[22] = 0.1e1 / zz[11];
zz[66] = zz[6] * zz[22];
zz[67] = vX_in[3] * zz[27] + vX_in[4] * zz[0];
zz[68] = -0.1e1 * (zz[6] * zz[25] - 0.1e1 * zz[11] * vX_in[5]) * zz[22];
zz[0] = 0.1e1 / zz[57];
zz[46] = zz[58] * zz[0];
zz[47] = zz[59] * zz[0];
zz[48] = -0.1e1 * zz[60] * zz[46] + zz[61];
zz[49] = -0.1e1 * zz[60] * zz[47] + zz[62];
zz[50] = -0.1e1 * zz[63] * zz[46] + zz[64];
zz[51] = -0.1e1 * zz[63] * zz[47] + zz[65];
zz[11] = 0.1e1 / zz[48];
zz[49] = zz[49] * zz[11];
zz[51] = -0.1e1 * zz[50] * zz[49] + zz[51];
zz[40] = zz[69] / zz[54];
zz[41] = zz[70] / zz[55];
zz[42] = zz[71] / zz[56];
zz[43] = zz[72] * zz[0];
zz[38] = -0.1e1 * zz[43] * zz[60] + zz[73];
zz[39] = -0.1e1 * zz[43] * zz[63] + zz[74];
zz[44] = zz[38] * zz[11];
zz[39] = -0.1e1 * zz[44] * zz[50] + zz[39];
zz[45] = zz[39] / zz[51];
zz[0] = -0.1e1 * zz[45] * zz[49] + zz[44];
zz[53] = zz[0];
zz[11] = -0.1e1 * zz[45] * zz[47] - 0.1e1 * zz[53] * zz[46] + zz[43];
zz[52] = zz[11];
XDOT[0] = zz[40];
XDOT[1] = zz[41];
XDOT[2] = zz[42];
XDOT[3] = zz[11];
XDOT[4] = zz[0];
XDOT[5] = zz[45];
XDOT[6] = zz[75];
XDOT[7] = zz[76];
XDOT[8] = zz[77];
XDOT[9] = zz[66];
XDOT[10] = zz[67];
XDOT[11] = zz[68];
}
#define DFPPROC_Xdot
#endif
 
