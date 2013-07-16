#include <math.h>
#include <stdlib.h>
 
#ifndef DFPEXTERNALS
//EXTERNAL FUNCTION DEFINITIONS
//Required External Function

double GonthierSCM (
  double *r1,
  double *r2,
  double *v1,
  double *v2,
  double *R1,
  double *R2,
  double *w1,
  double *w2,
  double *leParams,
  double *aOut)
{
  int i;
  int iterMax;
  double height;
  double x;
  double pieVal;
  double muVel;
  double wNormalMag;
  double wTangentMag;
  double velTangentMag;
  double vcn[3];
  double velRelative[3];
  double velNormal[3];
  double velTangent[3];
  double wRelative[3];
  double wTangent[3];
  double wNormal[3];
  double contactForceMag;
  double contactForce[3];
  double linearFrictionForce[3];
  double linearFrictionMoment[3];
  double spinFrictionTorque[3];
  double rollingResistance[3];
  double interpenVolume;
  double interpenJ[2];
  double f;
  double df;
  double xc;
  double pen;
  double muT;
  double muSpin;
  double radGy;
  double sphereCenter[3];
  double floorNormal[3];
  double floorPoint[3];
  double R;
  double k;
  double c;
  double floorTangent1[3];
  double floorTangent2[3];
  double muS;
  double muD;
  double velS;
  double velD;
  double delta;
  double pSphere[3];
  double normal[3];
  double tangent1[3];
  double tangent2[3];
  double pAction[3];
  double tmp1[3];
  double tmp2[3];
  pSphere[0] = 0;
  pSphere[1] = 0;
  pSphere[2] = 0;
  normal[0] = 0;
  normal[1] = 0;
  normal[2] = 0;
  tangent1[0] = 0;
  tangent1[1] = 0;
  tangent1[2] = 0;
  tangent2[0] = 0;
  tangent2[1] = 0;
  tangent2[2] = 0;
  pAction[0] = 0;
  pAction[1] = 0;
  pAction[2] = 0;
  tmp1[0] = 0;
  tmp1[1] = 0;
  tmp1[2] = 0;
  tmp2[0] = 0;
  tmp2[1] = 0;
  tmp2[2] = 0;
  vcn[0] = 0;
  vcn[1] = 0;
  vcn[2] = 0;
  velRelative[0] = 0;
  velRelative[1] = 0;
  velRelative[2] = 0;
  velNormal[0] = 0;
  velNormal[1] = 0;
  velNormal[2] = 0;
  velTangent[0] = 0;
  velTangent[1] = 0;
  velTangent[2] = 0;
  wRelative[0] = 0;
  wRelative[1] = 0;
  wRelative[2] = 0;
  wTangent[0] = 0;
  wTangent[1] = 0;
  wTangent[2] = 0;
  wNormal[0] = 0;
  wNormal[1] = 0;
  wNormal[2] = 0;
  contactForce[0] = 0;
  contactForce[1] = 0;
  contactForce[2] = 0;
  linearFrictionForce[0] = 0;
  linearFrictionForce[1] = 0;
  linearFrictionForce[2] = 0;
  spinFrictionTorque[0] = 0;
  spinFrictionTorque[1] = 0;
  spinFrictionTorque[2] = 0;
  rollingResistance[0] = 0;
  rollingResistance[1] = 0;
  rollingResistance[2] = 0;
  interpenJ[0] = 0;
  interpenJ[1] = 0;
  sphereCenter[0] = 0;
  sphereCenter[1] = 0;
  sphereCenter[2] = 0;
  floorNormal[0] = 0;
  floorNormal[1] = 0;
  floorNormal[2] = 0;
  floorPoint[0] = 0;
  floorPoint[1] = 0;
  floorPoint[2] = 0;
  floorTangent1[0] = 0;
  floorTangent1[1] = 0;
  floorTangent1[2] = 0;
  floorTangent2[0] = 0;
  floorTangent2[1] = 0;
  floorTangent2[2] = 0;
  linearFrictionMoment[0] = 0;
  linearFrictionMoment[1] = 0;
  linearFrictionMoment[2] = 0;
  for (i = 1; i <= 12; i++)
    aOut[i - 1] = 0.0e0;
  iterMax = 10000;
  sphereCenter[0] = leParams[0];
  sphereCenter[1] = leParams[1];
  sphereCenter[2] = leParams[2];
  floorNormal[0] = leParams[3];
  floorNormal[1] = leParams[4];
  floorNormal[2] = leParams[5];
  floorPoint[0] = leParams[6];
  floorPoint[1] = leParams[7];
  floorPoint[2] = leParams[8];
  floorTangent1[0] = leParams[9];
  floorTangent1[1] = leParams[10];
  floorTangent1[2] = leParams[11];
  floorTangent2[0] = leParams[12];
  floorTangent2[1] = leParams[13];
  floorTangent2[2] = leParams[14];
  R = leParams[15];
  k = leParams[16];
  c = leParams[17];
  muS = leParams[18];
  muD = leParams[19];
  velS = leParams[20];
  velD = leParams[21];
  pieVal = 0.31415926535897932384626433832795e1;
  pSphere[0] = sphereCenter[0] * R1[0] + sphereCenter[1] * R1[3] + sphereCenter[2] * R1[6];
  pSphere[1] = sphereCenter[0] * R1[1] + sphereCenter[1] * R1[4] + sphereCenter[2] * R1[7];
  pSphere[2] = sphereCenter[0] * R1[2] + sphereCenter[1] * R1[5] + sphereCenter[2] * R1[8];
  normal[0] = floorNormal[0] * R2[0] + floorNormal[1] * R2[3] + floorNormal[2] * R2[6];
  normal[1] = floorNormal[0] * R2[1] + floorNormal[1] * R2[4] + floorNormal[2] * R2[7];
  normal[2] = floorNormal[0] * R2[2] + floorNormal[1] * R2[5] + floorNormal[2] * R2[8];
  height = floorNormal[0] * floorPoint[0] + floorNormal[1] * floorPoint[1] + floorNormal[2] * floorPoint[2];
  pAction[0] = r1[0] + pSphere[0] - R * normal[0];
  pAction[1] = r1[1] + pSphere[1] - R * normal[1];
  pAction[2] = r1[2] + pSphere[2] - R * normal[2];
  x = (pAction[0] - r2[0]) * normal[0] + (pAction[1] - r2[1]) * normal[1] + (pAction[2] - r2[2]) * normal[2] - height;
  if (x < 0.0e0 && 0.0e0 < R)
  {
    pen = fabs(x);
    if (0.2e1 * R < pen)
      pen = 0.2e1 * R;
    pAction[0] = r1[0] + pSphere[0] - (R - pen) * normal[0];
    pAction[1] = r1[1] + pSphere[1] - (R - pen) * normal[1];
    pAction[2] = r1[2] + pSphere[2] - (R - pen) * normal[2];
    tangent1[0] = floorTangent1[0] * R2[0] + floorTangent1[1] * R2[3] + floorTangent1[2] * R2[6];
    tangent1[1] = floorTangent1[0] * R2[1] + floorTangent1[1] * R2[4] + floorTangent1[2] * R2[7];
    tangent1[2] = floorTangent1[0] * R2[2] + floorTangent1[1] * R2[5] + floorTangent1[2] * R2[8];
    tangent2[0] = floorTangent2[0] * R2[0] + floorTangent2[1] * R2[3] + floorTangent2[2] * R2[6];
    tangent2[1] = floorTangent2[0] * R2[1] + floorTangent2[1] * R2[4] + floorTangent2[2] * R2[7];
    tangent2[2] = floorTangent2[0] * R2[2] + floorTangent2[1] * R2[5] + floorTangent2[2] * R2[8];
    tmp1[0] = w1[1] * (pAction[2] - r1[2]) - w1[2] * (pAction[1] - r1[1]);
    tmp1[1] = w1[2] * (pAction[0] - r1[0]) - w1[0] * (pAction[2] - r1[2]);
    tmp1[2] = w1[0] * (pAction[1] - r1[1]) - w1[1] * (pAction[0] - r1[0]);
    tmp2[0] = w2[1] * (pAction[2] - r2[2]) - w2[2] * (pAction[1] - r2[1]);
    tmp2[1] = w2[2] * (pAction[0] - r2[0]) - w2[0] * (pAction[2] - r2[2]);
    tmp2[2] = w2[0] * (pAction[1] - r2[1]) - w2[1] * (pAction[0] - r2[0]);
    velRelative[0] = v1[0] + tmp1[0] - v2[0] - tmp2[0];
    velRelative[1] = v1[1] + tmp1[1] - v2[1] - tmp2[1];
    velRelative[2] = v1[2] + tmp1[2] - v2[2] - tmp2[2];
    velNormal[0] = velRelative[0] * normal[0];
    velNormal[1] = velRelative[1] * normal[1];
    velNormal[2] = velRelative[2] * normal[2];
    velTangent[0] = velRelative[0] * (tangent1[0] + tangent2[0]);
    velTangent[1] = velRelative[1] * (tangent1[1] + tangent2[1]);
    velTangent[2] = velRelative[2] * (tangent1[2] + tangent2[2]);
    wRelative[0] = w1[0] - w2[0];
    wRelative[1] = w1[1] - w2[1];
    wRelative[2] = w1[2] - w2[2];
    wTangent[0] = wRelative[0] * (tangent1[0] + tangent2[0]);
    wTangent[1] = wRelative[1] * (tangent1[1] + tangent2[1]);
    wTangent[2] = wRelative[2] * (tangent1[2] + tangent2[2]);
    wNormal[0] = wRelative[0] * normal[0];
    wNormal[1] = wRelative[1] * normal[1];
    wNormal[2] = wRelative[2] * normal[2];
    velTangentMag = pow(pow(velTangent[0], 0.2e1) + pow(velTangent[1], 0.2e1) + pow(velTangent[2], 0.2e1), 0.5e0);
    wNormalMag = pow(pow(wNormal[0], 0.2e1) + pow(wNormal[1], 0.2e1) + pow(wNormal[2], 0.2e1), 0.5e0);
    wTangentMag = pow(pow(wTangent[0], 0.2e1) + pow(wTangent[1], 0.2e1) + pow(wTangent[2], 0.2e1), 0.5e0);
    interpenVolume = 0.1e1 / 0.3e1 * pieVal * pen * pen * (0.3e1 * R - pen);
    xc = 0.5e0 * pen;
    df = 0.0e0;
    f = interpenVolume;
    i = 0;
    while (0.1e-3 * interpenVolume < fabs(f) && i < iterMax)
    {
      f = pieVal * xc * xc * (0.3e1 * R - xc) / 0.3e1 - 0.5e0 * interpenVolume;
      df = pieVal * xc * (0.6e1 * R - 0.3e1 * xc) / 0.3e1;
      xc = xc - f / df;
      i = i + 1;
    }
    if (iterMax <= i)
    {
      pen = 0.1e-2 * R;
      xc = 0.5e0 * pen;
    }
    interpenJ[0] = 0.1e1 / 0.30e2 * pieVal * pow(pen, 0.3e1) * (0.3e1 * pen * pen - 0.15e2 * R * pen + 0.20e2 * R * R);
    interpenJ[1] = 0.1e1 / 0.60e2 * pieVal * pen * pen * (-0.9e1 * pow(pen, 0.3e1) + 0.15e2 * R * pen * pen + 0.30e2 * pen * pen * xc - 0.20e2 * pen * xc * xc + 0.20e2 * R * R * pen - 0.80e2 * pen * R * xc + 0.60e2 * xc * xc * R);
    vcn[0] = -velNormal[0];
    vcn[1] = -velNormal[1];
    vcn[2] = -velNormal[2];
    if (c * (vcn[0] + vcn[1] + vcn[2]) < -0.1e1)
    {
      contactForce[0] = 0.0e0;
      contactForce[1] = 0.0e0;
      contactForce[2] = 0.0e0;
    }
    else
    {
      contactForce[0] = k * interpenVolume * (0.1e1 + vcn[0] * c) * normal[0];
      contactForce[1] = k * interpenVolume * (0.1e1 + vcn[1] * c) * normal[1];
      contactForce[2] = k * interpenVolume * (0.1e1 + vcn[2] * c) * normal[2];
    }
    contactForceMag = pow(pow(contactForce[0], 0.2e1) + pow(contactForce[1], 0.2e1) + pow(contactForce[2], 0.2e1), 0.5e0);
    if (0.0e0 < wTangentMag)
    {
      rollingResistance[0] = -k * c * interpenJ[1] * wTangent[0];
      rollingResistance[1] = -k * c * interpenJ[1] * wTangent[1];
      rollingResistance[2] = -k * c * interpenJ[1] * wTangent[2];
    }
    else
    {
      rollingResistance[0] = 0.0e0;
      rollingResistance[1] = 0.0e0;
      rollingResistance[2] = 0.0e0;
    }
    muVel = velTangentMag;
    if (0.0e0 < muVel)
      if (velD < muVel)
        muT = muD;
      else if (velS < muVel && muVel < velD)
      {
        delta = (muVel - velS) / (velD - velS);
        muT = muS + (muD - muS) * delta * delta * (0.30e1 - 0.20e1 * delta);
      }
      else
      {
        delta = (muVel + velS) * 0.5000000000e0 / velS;
        muT = -muS + 0.2e1 * muS * delta * delta * (0.30e1 - 0.20e1 * delta);
      }
    if (0.0e0 < velTangentMag)
    {
      linearFrictionForce[0] = -muT * contactForceMag * velTangent[0];
      linearFrictionForce[1] = -muT * contactForceMag * velTangent[1];
      linearFrictionForce[2] = -muT * contactForceMag * velTangent[2];
    }
    else
    {
      linearFrictionForce[0] = 0.0e0;
      linearFrictionForce[1] = 0.0e0;
      linearFrictionForce[2] = 0.0e0;
    }
    linearFrictionMoment[0] = (linearFrictionForce[1] * normal[2] - linearFrictionForce[2] * normal[1]) * (R - pen);
    linearFrictionMoment[1] = (-linearFrictionForce[0] * normal[2] + linearFrictionForce[2] * normal[0]) * (R - pen);
    linearFrictionMoment[2] = (linearFrictionForce[0] * normal[1] - linearFrictionForce[1] * normal[0]) * (R - pen);
    radGy = sqrt(interpenJ[0] / interpenVolume);
    muVel = radGy * wNormalMag;
    if (0.0e0 < muVel)
      if (velD < muVel)
        muSpin = muD;
      else if (velS < muVel && muVel < velD)
      {
        delta = (muVel - velS) / (velD - velS);
        muSpin = muS + (muD - muS) * delta * delta * (0.30e1 - 0.20e1 * delta);
      }
      else
      {
        delta = (muVel + velS) * 0.5000000000e0 / velS;
        muSpin = -muS + 0.2e1 * muS * delta * delta * (0.30e1 - 0.20e1 * delta);
      }
    if (0.0e0 < wNormalMag)
    {
      spinFrictionTorque[0] = -muSpin * contactForceMag * radGy * wNormal[0];
      spinFrictionTorque[1] = -muSpin * contactForceMag * radGy * wNormal[1];
      spinFrictionTorque[2] = -muSpin * contactForceMag * radGy * wNormal[2];
    }
    else
    {
      spinFrictionTorque[0] = 0.0e0;
      spinFrictionTorque[1] = 0.0e0;
      spinFrictionTorque[2] = 0.0e0;
    }
    aOut[0] = contactForce[0] + linearFrictionForce[0];
    aOut[1] = contactForce[1] + linearFrictionForce[1];
    aOut[2] = contactForce[2] + linearFrictionForce[2];
    aOut[3] = rollingResistance[0] + spinFrictionTorque[0] + linearFrictionMoment[0];
    aOut[4] = rollingResistance[1] + spinFrictionTorque[1] + linearFrictionMoment[1];
    aOut[5] = rollingResistance[2] + spinFrictionTorque[2] + linearFrictionMoment[2];
    aOut[6] = -aOut[0];
    aOut[7] = -aOut[1];
    aOut[8] = -aOut[2];
    aOut[9] = -aOut[3];
    aOut[10] = -aOut[4];
    aOut[11] = -aOut[5];
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
double CParamsContact11[22];
 
CParamsContact11[0] = vParam[6];
CParamsContact11[1] = vParam[7];
CParamsContact11[2] = vParam[8];
CParamsContact11[3] = vParam[13];
CParamsContact11[4] = vParam[14];
CParamsContact11[5] = vParam[15];
CParamsContact11[6] = vParam[16];
CParamsContact11[7] = vParam[17];
CParamsContact11[8] = vParam[18];
CParamsContact11[9] = vParam[19];
CParamsContact11[10] = vParam[21];
CParamsContact11[11] = vParam[23];
CParamsContact11[12] = vParam[20];
CParamsContact11[13] = vParam[22];
CParamsContact11[14] = vParam[24];
CParamsContact11[15] = vParam[4];
CParamsContact11[16] = vParam[9];
CParamsContact11[17] = vParam[5];
CParamsContact11[18] = vParam[12];
CParamsContact11[19] = vParam[11];
CParamsContact11[20] = vParam[26];
CParamsContact11[21] = vParam[25];
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
zz[31] = GonthierSCM(r1Contact11, r2Contact11, v1Contact11, v2Contact11, R1Contact11, R2Contact11, w1Contact11, w2Contact11, CParamsContact11, CFContact11);
zz[32] = CFContact11[0];
zz[33] = CFContact11[1];
zz[34] = CFContact11[2];
zz[35] = CFContact11[3];
zz[36] = CFContact11[4];
zz[37] = CFContact11[5];
zz[54] = vParam[10];
zz[55] = vParam[10];
zz[56] = vParam[10];
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
zz[70] = -0.1e1 * vParam[0] * vParam[10] + zz[33];
zz[71] = zz[34];
zz[13] = zz[3] * w2Contact11[1] - 0.1e1 * zz[30] * w2Contact11[2];
zz[3] = -0.1e1 * zz[3] * w2Contact11[0] + zz[7] * w2Contact11[2];
zz[30] = -0.1e1 * zz[7] * w2Contact11[1] + zz[30] * w2Contact11[0];
zz[7] = zz[8] * w1Contact11[2] + zz[12] * w1Contact11[0] + zz[15] * w1Contact11[1];
zz[16] = zz[10] * w1Contact11[0] + zz[14] * w1Contact11[1] + zz[15] * w1Contact11[2];
zz[17] = zz[9] * w1Contact11[0] + zz[10] * w1Contact11[1] + zz[12] * w1Contact11[2];
zz[9] = -0.1e1 * zz[3] * zz[10] - 0.1e1 * zz[7] * w1Contact11[1] - 0.1e1 * zz[9] * zz[13] - 0.1e1 * zz[12] * zz[30] + zz[16] * w1Contact11[2] + zz[35];
zz[7] = -0.1e1 * zz[3] * zz[14] + zz[7] * w1Contact11[0] - 0.1e1 * zz[10] * zz[13] - 0.1e1 * zz[15] * zz[30] - 0.1e1 * zz[17] * w1Contact11[2] + zz[36];
zz[30] = -0.1e1 * zz[3] * zz[15] - 0.1e1 * zz[8] * zz[30] - 0.1e1 * zz[12] * zz[13] - 0.1e1 * zz[16] * w1Contact11[0] + zz[17] * w1Contact11[1] + zz[37];
zz[72] = zz[1] * zz[7] + zz[4] * zz[30] + zz[9] * zz[22];
zz[73] = zz[2] * zz[7] - 0.1e1 * zz[5] * zz[9] + zz[29] * zz[30];
zz[74] = zz[9] * zz[25] + zz[11] * (-0.1e1 * zz[7] * zz[26] + zz[28] * zz[30]);
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
 
