function GrabnerDiskContactreturn = GrabnerDiskContact(r20, r10, v20, v10, R20, R10, w20, w10, leParams, aOut)
  for i = 1:12
    aOut(i) = 0;
  end
  epsRoot = 0.1490116119384766e-7;
  sx = leParams(1);
  sy = leParams(2);
  sz = leParams(3);
  nx = leParams(4);
  ny = leParams(5);
  nz = leParams(6);
  t1x = leParams(7);
  t1y = leParams(8);
  t1z = leParams(9);
  t2x = leParams(10);
  t2y = leParams(11);
  t2z = leParams(12);
  r = leParams(13);
  a = leParams(14);
  k = leParams(15);
  p = leParams(16);
  c = leParams(17);
  mus = leParams(18);
  mud = leParams(19);
  stVel = leParams(20);
  dyVel = leParams(21);
  nP(1) = R10(1) * nx + R10(2) * ny + R10(3) * nz;
  nP(2) = R10(4) * nx + R10(5) * ny + R10(6) * nz;
  nP(3) = R10(7) * nx + R10(8) * ny + R10(9) * nz;
  tP1(1) = R10(1) * t1x + R10(2) * t1y + R10(3) * t1z;
  tP1(2) = R10(4) * t1x + R10(5) * t1y + R10(6) * t1z;
  tP1(3) = R10(7) * t1x + R10(8) * t1y + R10(9) * t1z;
  tP2(1) = R10(1) * t2x + R10(2) * t2y + R10(3) * t2z;
  tP2(2) = R10(4) * t2x + R10(5) * t2y + R10(6) * t2z;
  tP2(3) = R10(7) * t2x + R10(8) * t2y + R10(9) * t2z;
  rO(1) = R10(1) * sx + R10(2) * sy + R10(3) * sz + r10(1);
  rO(2) = R10(4) * sx + R10(5) * sy + R10(6) * sz + r10(2);
  rO(3) = R10(7) * sx + R10(8) * sy + R10(9) * sz + r10(3);
  rM(1) = r20(1);
  rM(2) = r20(2);
  rM(3) = r20(3);
  rMO(1) = rM(1) - rO(1);
  rMO(2) = rM(2) - rO(2);
  rMO(3) = rM(3) - rO(3);
  nC(1) = R20(7);
  nC(2) = R20(8);
  nC(3) = R20(9);
  tmp1(1) = nC(2) * nP(3) - nC(3) * nP(2);
  tmp1(2) = -nC(1) * nP(3) + nC(3) * nP(1);
  tmp1(3) = nC(1) * nP(2) - nC(2) * nP(1);
  tmp2(1) = nC(2) * tmp1(3) - nC(3) * tmp1(2);
  tmp2(2) = -nC(1) * tmp1(3) + nC(3) * tmp1(1);
  tmp2(3) = nC(1) * tmp1(2) - nC(2) * tmp1(1);
  tmpf1 = tmp2(1) * nP(1) + tmp2(2) * nP(2) + tmp2(3) * nP(3);
  if (0.0e0 < tmpf1)
    tmp2(1) = -0.10e1 * tmp2(1);
    tmp2(2) = -0.10e1 * tmp2(2);
    tmp2(3) = -0.10e1 * tmp2(3);
  end
  tmpf1 = (tmp1(1) * tmp1(1) + tmp1(2) * tmp1(2) + tmp1(3) * tmp1(3)) ^ 0.5e0;
  if (tmpf1 < epsRoot)
    tmpf1 = epsRoot;
  end
  alpha = acos(abs(nP(1) * nC(1) + nP(2) * nC(2) + nP(3) * nC(3)));
  radius = r * (0.1e1 - exp(-a * sin(alpha)));
  rPMuv(1) = tmp2(1) / tmpf1;
  rPMuv(2) = tmp2(2) / tmpf1;
  rPMuv(3) = tmp2(3) / tmpf1;
  rPM(1) = radius * rPMuv(1);
  rPM(2) = radius * rPMuv(2);
  rPM(3) = radius * rPMuv(3);
  rPMO(1) = rMO(1) + rPM(1);
  rPMO(2) = rMO(2) + rPM(2);
  rPMO(3) = rMO(3) + rPM(3);
  gN = rPMO(1) * nP(1) + rPMO(2) * nP(2) + rPMO(3) * nP(3);
  if (gN < 0.0e0)
    vP20(1) = -w20(3) * rPM(2) + w20(2) * rPM(3) + v20(1);
    vP20(2) = w20(3) * rPM(1) - w20(1) * rPM(3) + v20(2);
    vP20(3) = -w20(2) * rPM(1) + w20(1) * rPM(2) + v20(3);
    tmpf1 = rPMO(1) * tP1(1) + rPMO(2) * tP1(2) + rPMO(3) * tP1(3);
    tmpf2 = rPMO(1) * tP2(1) + rPMO(2) * tP2(2) + rPMO(3) * tP2(3);
    rQO(1) = tmpf1 * tP1(1) + tmpf2 * tP2(1);
    rQO(2) = tmpf1 * tP1(2) + tmpf2 * tP2(2);
    rQO(3) = tmpf1 * tP1(3) + tmpf2 * tP2(3);
    vQ10(1) = -w10(3) * rQO(2) + w10(2) * rQO(3) + v10(1);
    vQ10(2) = w10(3) * rQO(1) - w10(1) * rQO(3) + v10(2);
    vQ10(3) = -w10(2) * rQO(1) + w10(1) * rQO(2) + v10(3);
    vPQ(1) = vP20(1) - vQ10(1);
    vPQ(2) = vP20(2) - vQ10(2);
    vPQ(3) = vP20(3) - vQ10(3);
    velNormal = vPQ(1) * nP(1) + vPQ(2) * nP(2) + vPQ(3) * nP(3);
    velTangent(1) = vPQ(1) * (tP1(1) + tP2(1));
    velTangent(2) = vPQ(2) * (tP1(2) + tP2(2));
    velTangent(3) = vPQ(3) * (tP1(3) + tP2(3));
    velTangentMag = (velTangent(1) * velTangent(1) + velTangent(2) * velTangent(2) + velTangent(3) * velTangent(3)) ^ 0.5e0;
    mu = 0.0e0;
    delta = 0.0e0;
    if (0.0e0 < velTangentMag)
      if (dyVel < velTangentMag)
        mu = mud;
      elseif (stVel < velTangentMag & velTangentMag < dyVel)
        delta = (velTangentMag - stVel) / (dyVel - stVel);
        mu = mus + (mud - mus) * delta * delta * (0.30e1 - 0.20e1 * delta);
      else
        delta = (velTangentMag + stVel) / 0.20e1 / stVel;
        mu = -mus + 0.2e1 * mus * delta * delta * (0.30e1 - 0.20e1 * delta);
      end
    end
    forceNormal = k * abs(gN) ^ p * (0.1e1 - c * velNormal);
    if (forceNormal < 0.0e0)
      forceNormal = 0.0e0;
    end
    forceTangentMag = -mu * forceNormal;
    velEps = 0.10000e4 * epsRoot;
    if (velEps < velTangentMag)
      forceTangentUV(1) = velTangent(1) / velTangentMag;
      forceTangentUV(2) = velTangent(2) / velTangentMag;
      forceTangentUV(3) = velTangent(3) / velTangentMag;
    else
      tmpf1 = velTangent(1) / velEps;
      forceTangentUV(1) = tmpf1 * (0.30e1 / 0.20e1 * abs(tmpf1) - 0.10e1 / 0.20e1 * abs(tmpf1 * tmpf1 * tmpf1));
      tmpf1 = velTangent(2) / velEps;
      forceTangentUV(2) = tmpf1 * (0.30e1 / 0.20e1 * abs(tmpf1) - 0.10e1 / 0.20e1 * abs(tmpf1 * tmpf1 * tmpf1));
      tmpf1 = velTangent(3) / velEps;
      forceTangentUV(3) = tmpf1 * (0.30e1 / 0.20e1 * abs(tmpf1) - 0.10e1 / 0.20e1 * abs(tmpf1 * tmpf1 * tmpf1));
    end
    rQM(1) = rQO(1) - rMO(1);
    rQM(2) = rQO(2) - rMO(2);
    rQM(3) = rQO(3) - rMO(3);
    aOut(1) = (forceNormal * nP(1) + forceTangentMag * forceTangentUV(1));
    aOut(2) = (forceNormal * nP(2) + forceTangentMag * forceTangentUV(2));
    aOut(3) = (forceNormal * nP(3) + forceTangentMag * forceTangentUV(3));
    aOut(4) = (rPM(2) * aOut(3) - rPM(3) * aOut(2));
    aOut(5) = (-rPM(1) * aOut(3) + rPM(3) * aOut(1));
    aOut(6) = (rPM(1) * aOut(2) - rPM(2) * aOut(1));
    aOut(7) = -aOut(1);
    aOut(8) = -aOut(2);
    aOut(9) = -aOut(3);
    aOut(10) = (rPMO(2) * aOut(9) - rPMO(3) * aOut(8));
    aOut(11) = (-rPMO(1) * aOut(9) + rPMO(3) * aOut(7));
    aOut(12) = (rPMO(1) * aOut(8) - rPMO(2) * aOut(7));
  end
  GrabnerDiskContactreturn = 0.0e0;
