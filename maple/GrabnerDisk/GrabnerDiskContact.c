#include <math.h>
#include <stdio.h>

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
  printf("sx, sy, sz  : %f,%f,%f \n", sx, sy, sz);
  printf("nx, ny, nz  : %f,%f,%f \n", nx, ny, nz);
  printf("t1x,t1y,t1z : %f,%f,%f \n", t1x, t1y, t1z);
  printf("t2x,t2y,t2z : %f,%f,%f \n", t2x, t2y, t2z);
  printf("  r,  a,    : %f,%f \n", r, a);
  printf("  k,  p,  c : %f,%f,%f \n", k, p, c);
  printf("mus,mud     : %f,%f \n", mus, mud);
  printf("stVel,dyVel : %f,%f \n\n", stVel, dyVel);
  nP[0] = R10[0] * nx + R10[1] * ny + R10[2] * nz;
  nP[1] = R10[3] * nx + R10[4] * ny + R10[5] * nz;
  nP[2] = R10[6] * nx + R10[7] * ny + R10[8] * nz;
  printf("nP: %f,%f,%f\n", nP[0], nP[1], nP[2]);
  tP1[0] = R10[0] * t1x + R10[1] * t1y + R10[2] * t1z;
  tP1[1] = R10[3] * t1x + R10[4] * t1y + R10[5] * t1z;
  tP1[2] = R10[6] * t1x + R10[7] * t1y + R10[8] * t1z;
  printf("tP1: %f,%f,%f\n", tP1[0], tP1[1], tP1[2]);
  tP2[0] = R10[0] * t2x + R10[1] * t2y + R10[2] * t2z;
  tP2[1] = R10[3] * t2x + R10[4] * t2y + R10[5] * t2z;
  tP2[2] = R10[6] * t2x + R10[7] * t2y + R10[8] * t2z;
  printf("tP2: %f,%f,%f\n", tP2[0], tP2[1], tP2[2]);
  rO[0] = R10[0] * sx + R10[1] * sy + R10[2] * sz + r10[0];
  rO[1] = R10[3] * sx + R10[4] * sy + R10[5] * sz + r10[1];
  rO[2] = R10[6] * sx + R10[7] * sy + R10[8] * sz + r10[2];
  printf("r0: %f,%f,%f\n\n", rO[0], rO[1], rO[2]);
  rM[0] = r20[0];
  rM[1] = r20[1];
  rM[2] = r20[2];
  printf("rM: %f,%f,%f\n", rM[0], rM[1], rM[2]);
  rMO[0] = rM[0] - rO[0];
  rMO[1] = rM[1] - rO[1];
  rMO[2] = rM[2] - rO[2];
  printf("rMO: %f,%f,%f\n\n", rMO[0], rMO[1], rMO[2]);
  nC[0] = R20[6];
  nC[1] = R20[7];
  nC[2] = R20[8];
  printf("nC: %f,%f,%f\n", nC[0], nC[1], nC[2]);
  tmp1[0] = nC[1] * nP[2] - nC[2] * nP[1];
  tmp1[1] = -nC[0] * nP[2] + nC[2] * nP[0];
  tmp1[2] = nC[0] * nP[1] - nC[1] * nP[0];
  printf("nC x nP: %f,%f,%f\n", tmp1[0], tmp1[1], tmp1[2]);
  tmp2[0] = nC[1] * tmp1[2] - nC[2] * tmp1[1];
  tmp2[1] = -nC[0] * tmp1[2] + nC[2] * tmp1[0];
  tmp2[2] = nC[0] * tmp1[1] - nC[1] * tmp1[0];
  printf("nC x nC x nP: %f,%f,%f\n\n", tmp2[0], tmp2[1], tmp2[2]);
  tmpf1 = tmp2[0] * nP[0] + tmp2[1] * nP[1] + tmp2[2] * nP[2];
  if (0.0e0 < tmpf1)
  {
    tmp2[0] = -0.10e1 * tmp2[0];
    tmp2[1] = -0.10e1 * tmp2[1];
    tmp2[2] = -0.10e1 * tmp2[2];
  }
  tmpf1 = pow(tmp1[0] * tmp1[0] + tmp1[1] * tmp1[1] + tmp1[2] * tmp1[2], 0.5e0);
  if (tmpf1 < epsRoot)
    tmpf1 = epsRoot;
  printf("tmpf1: %f\n", tmpf1);
  alpha = acos(fabs(nP[0] * nC[0] + nP[1] * nC[1] + nP[2] * nC[2]));
  printf("alpha : %f\n", alpha);
  radius = r * (0.1e1 - exp(-a * sin(alpha)));
  printf("radius: %f\n", radius);
  rPMuv[0] = tmp2[0] / tmpf1;
  rPMuv[1] = tmp2[1] / tmpf1;
  rPMuv[2] = tmp2[2] / tmpf1;
  printf("rPMuv: %f,%f,%f\n", rPMuv[0], rPMuv[1], rPMuv[2]);
  rPM[0] = radius * rPMuv[0];
  rPM[1] = radius * rPMuv[1];
  rPM[2] = radius * rPMuv[2];
  printf("rPM: %f,%f,%f\n", rPM[0], rPM[1], rPM[2]);
  rPMO[0] = rMO[0] + rPM[0];
  rPMO[1] = rMO[1] + rPM[1];
  rPMO[2] = rMO[2] + rPM[2];
  printf("rPMO: %f,%f,%f\n\n", rPMO[0], rPMO[1], rPMO[2]);
  gN = rPMO[0] * nP[0] + rPMO[1] * nP[1] + rPMO[2] * nP[2];
  printf("gN: %f\n\n", gN);
  if (gN < 0.0e0)
  {
    vP20[0] = -w20[2] * rPM[1] + w20[1] * rPM[2] + v20[0];
    vP20[1] = w20[2] * rPM[0] - w20[0] * rPM[2] + v20[1];
    vP20[2] = -w20[1] * rPM[0] + w20[0] * rPM[1] + v20[2];
    printf("vP20: %f,%f,%f\n", vP20[0], vP20[1], vP20[2]);
    tmpf1 = rPMO[0] * tP1[0] + rPMO[1] * tP1[1] + rPMO[2] * tP1[2];
    tmpf2 = rPMO[0] * tP2[0] + rPMO[1] * tP2[1] + rPMO[2] * tP2[2];
    rQO[0] = tmpf1 * tP1[0] + tmpf2 * tP2[0];
    rQO[1] = tmpf1 * tP1[1] + tmpf2 * tP2[1];
    rQO[2] = tmpf1 * tP1[2] + tmpf2 * tP2[2];
    printf("rQO: %f,%f,%f\n", rQO[0], rQO[1], rQO[2]);
    vQ10[0] = -w10[2] * rQO[1] + w10[1] * rQO[2] + v10[0];
    vQ10[1] = w10[2] * rQO[0] - w10[0] * rQO[2] + v10[1];
    vQ10[2] = -w10[1] * rQO[0] + w10[0] * rQO[1] + v10[2];
    printf("vQ10: %f,%f,%f\n", vQ10[0], vQ10[1], vQ10[2]);
    vPQ[0] = vP20[0] - vQ10[0];
    vPQ[1] = vP20[1] - vQ10[1];
    vPQ[2] = vP20[2] - vQ10[2];
    printf("vPQ: %f,%f,%f\n", vPQ[0], vPQ[1], vPQ[2]);
    velNormal = vPQ[0] * nP[0] + vPQ[1] * nP[1] + vPQ[2] * nP[2];
    printf("velNormal: %f\n", velNormal);
    velTangent[0] = vPQ[0] * (tP1[0] + tP2[0]);
    velTangent[1] = vPQ[1] * (tP1[1] + tP2[1]);
    velTangent[2] = vPQ[2] * (tP1[2] + tP2[2]);
    printf("velTangent: %f,%f,%f \n", velTangent[0], velTangent[1], velTangent[2]);
    velTangentMag = pow(velTangent[0] * velTangent[0] + velTangent[1] * velTangent[1] + velTangent[2] * velTangent[2], 0.5e0);
    printf("velTangentMag: %f\n\n", velTangentMag);
    mu = 0.0e0;
    delta = 0.0e0;
    if (0.0e0 < velTangentMag)
      if (dyVel < velTangentMag)
      {
        mu = mud;
        printf("Dynamic Region\n");
      }
      else if (stVel < velTangentMag && velTangentMag < dyVel)
      {
        delta = (velTangentMag - stVel) / (dyVel - stVel);
        mu = mus + (mud - mus) * delta * delta * (0.30e1 - 0.20e1 * delta);
        printf("Static to Dynamic Region\n");
      }
      else
      {
        delta = (velTangentMag + stVel) / 0.20e1 / stVel;
        mu = -mus + 0.2e1 * mus * delta * delta * (0.30e1 - 0.20e1 * delta);
        printf("Static Region\n");
      }
    printf("mu(%f), mus(%f), mud(%f), vels(%f), veld(%f)\n\n", mu, mus, mud, stVel, dyVel);
    forceNormal = k * pow(fabs(gN), p) * (0.1e1 - c * velNormal);
    if (forceNormal < 0.0e0)
      forceNormal = 0.0e0;
    printf("Normal Force: %f\n", forceNormal);
    forceTangentMag = -mu * forceNormal;
    printf("Tangential Force: %f\n", forceTangentMag);
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
      forceTangentUV[0] = tmpf1 * (0.30e1 / 0.20e1 * fabs(tmpf1) - 0.10e1 / 0.20e1 * fabs(tmpf1 * tmpf1 * tmpf1));
      tmpf1 = velTangent[1] / velEps;
      forceTangentUV[1] = tmpf1 * (0.30e1 / 0.20e1 * fabs(tmpf1) - 0.10e1 / 0.20e1 * fabs(tmpf1 * tmpf1 * tmpf1));
      tmpf1 = velTangent[2] / velEps;
      forceTangentUV[2] = tmpf1 * (0.30e1 / 0.20e1 * fabs(tmpf1) - 0.10e1 / 0.20e1 * fabs(tmpf1 * tmpf1 * tmpf1));
    }
    rQM[0] = rQO[0] - rMO[0];
    rQM[1] = rQO[1] - rMO[1];
    rQM[2] = rQO[2] - rMO[2];
    aOut[0] = (int) (forceNormal * nP[0] + forceTangentMag * forceTangentUV[0]);
    aOut[1] = (int) (forceNormal * nP[1] + forceTangentMag * forceTangentUV[1]);
    aOut[2] = (int) (forceNormal * nP[2] + forceTangentMag * forceTangentUV[2]);
    printf("Forces on disk: %f, %f, %f\n", aOut[0], aOut[1], aOut[2]);
    aOut[3] = (int) (rPM[1] * (double) aOut[2] - rPM[2] * (double) aOut[1]);
    aOut[4] = (int) (-rPM[0] * (double) aOut[2] + rPM[2] * (double) aOut[0]);
    aOut[5] = (int) (rPM[0] * (double) aOut[1] - rPM[1] * (double) aOut[0]);
    printf("Moments on disk: %f, %f, %f\n", aOut[3], aOut[4], aOut[5]);
    aOut[6] = -aOut[0];
    aOut[7] = -aOut[1];
    aOut[8] = -aOut[2];
    printf("Forces on Ground: %f, %f, %f\n", aOut[6], aOut[7], aOut[8]);
    aOut[9] = (int) (rPMO[1] * (double) aOut[8] - rPMO[2] * (double) aOut[7]);
    aOut[10] = (int) (-rPMO[0] * (double) aOut[8] + rPMO[2] * (double) aOut[6]);
    aOut[11] = (int) (rPMO[0] * (double) aOut[7] - rPMO[1] * (double) aOut[6]);
    printf("Moments on Ground: %f, %f, %f\n", aOut[9], aOut[10], aOut[11]);
  }
  return(0.0e0);
}
