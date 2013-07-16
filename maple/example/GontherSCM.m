function GonthierSCMreturn = GonthierSCM(r1, r2, v1, v2, R1, R2, w1, w2, leParams, aOut)
  pSphere = [0 0 0];
  normal = [0 0 0];
  tangent1 = [0 0 0];
  tangent2 = [0 0 0];
  pAction = [0 0 0];
  tmp1 = [0 0 0];
  tmp2 = [0 0 0];
  vcn = [0 0 0];
  velRelative = [0 0 0];
  velNormal = [0 0 0];
  velTangent = [0 0 0];
  wRelative = [0 0 0];
  wTangent = [0 0 0];
  wNormal = [0 0 0];
  contactForce = [0 0 0];
  linearFrictionForce = [0 0 0];
  spinFrictionTorque = [0 0 0];
  rollingResistance = [0 0 0];
  interpenJ = [0 0];
  sphereCenter = [0 0 0];
  floorNormal = [0 0 0];
  floorPoint = [0 0 0];
  floorTangent1 = [0 0 0];
  floorTangent2 = [0 0 0];
  linearFrictionMoment = [0 0 0];
  for i = 1:12
    aOut(i) = 0.0e0;
  end
  iterMax = 10000;
  sphereCenter(1) = leParams(1);
  sphereCenter(2) = leParams(2);
  sphereCenter(3) = leParams(3);
  floorNormal(1) = leParams(4);
  floorNormal(2) = leParams(5);
  floorNormal(3) = leParams(6);
  floorPoint(1) = leParams(7);
  floorPoint(2) = leParams(8);
  floorPoint(3) = leParams(9);
  floorTangent1(1) = leParams(10);
  floorTangent1(2) = leParams(11);
  floorTangent1(3) = leParams(12);
  floorTangent2(1) = leParams(13);
  floorTangent2(2) = leParams(14);
  floorTangent2(3) = leParams(15);
  R = leParams(16);
  k = leParams(17);
  c = leParams(18);
  muS = leParams(19);
  muD = leParams(20);
  velS = leParams(21);
  velD = leParams(22);
  pieVal = 0.31415926535897932384626433832795e1;
  pSphere(1) = R1(1) * sphereCenter(1) + R1(4) * sphereCenter(2) + R1(7) * sphereCenter(3);
  pSphere(2) = R1(2) * sphereCenter(1) + R1(5) * sphereCenter(2) + R1(8) * sphereCenter(3);
  pSphere(3) = R1(3) * sphereCenter(1) + R1(6) * sphereCenter(2) + R1(9) * sphereCenter(3);
  normal(1) = R2(1) * floorNormal(1) + R2(4) * floorNormal(2) + R2(7) * floorNormal(3);
  normal(2) = R2(2) * floorNormal(1) + R2(5) * floorNormal(2) + R2(8) * floorNormal(3);
  normal(3) = R2(3) * floorNormal(1) + R2(6) * floorNormal(2) + R2(9) * floorNormal(3);
  height = floorNormal(1) * floorPoint(1) + floorNormal(2) * floorPoint(2) + floorNormal(3) * floorPoint(3);
  pAction(1) = r1(1) + pSphere(1) - R * normal(1);
  pAction(2) = r1(2) + pSphere(2) - R * normal(2);
  pAction(3) = r1(3) + pSphere(3) - R * normal(3);
  x = (pAction(1) - r2(1)) * normal(1) + (pAction(2) - r2(2)) * normal(2) + (pAction(3) - r2(3)) * normal(3) - height;
  if (x < 0.0e0 & 0.0e0 < R)
    pen = abs(x);
    if (0.2e1 * R < pen)
      pen = 0.2e1 * R;
    end
    pAction(1) = r1(1) + pSphere(1) - (R - pen) * normal(1);
    pAction(2) = r1(2) + pSphere(2) - (R - pen) * normal(2);
    pAction(3) = r1(3) + pSphere(3) - (R - pen) * normal(3);
    tangent1(1) = R2(1) * floorTangent1(1) + R2(4) * floorTangent1(2) + R2(7) * floorTangent1(3);
    tangent1(2) = R2(2) * floorTangent1(1) + R2(5) * floorTangent1(2) + R2(8) * floorTangent1(3);
    tangent1(3) = R2(3) * floorTangent1(1) + R2(6) * floorTangent1(2) + R2(9) * floorTangent1(3);
    tangent2(1) = R2(1) * floorTangent2(1) + R2(4) * floorTangent2(2) + R2(7) * floorTangent2(3);
    tangent2(2) = R2(2) * floorTangent2(1) + R2(5) * floorTangent2(2) + R2(8) * floorTangent2(3);
    tangent2(3) = R2(3) * floorTangent2(1) + R2(6) * floorTangent2(2) + R2(9) * floorTangent2(3);
    tmp1(1) = w1(2) * (pAction(3) - r1(3)) - w1(3) * (pAction(2) - r1(2));
    tmp1(2) = w1(3) * (pAction(1) - r1(1)) - w1(1) * (pAction(3) - r1(3));
    tmp1(3) = w1(1) * (pAction(2) - r1(2)) - w1(2) * (pAction(1) - r1(1));
    tmp2(1) = w2(2) * (pAction(3) - r2(3)) - w2(3) * (pAction(2) - r2(2));
    tmp2(2) = w2(3) * (pAction(1) - r2(1)) - w2(1) * (pAction(3) - r2(3));
    tmp2(3) = w2(1) * (pAction(2) - r2(2)) - w2(2) * (pAction(1) - r2(1));
    velRelative(1) = v1(1) + tmp1(1) - v2(1) - tmp2(1);
    velRelative(2) = v1(2) + tmp1(2) - v2(2) - tmp2(2);
    velRelative(3) = v1(3) + tmp1(3) - v2(3) - tmp2(3);
    velNormal(1) = velRelative(1) * normal(1);
    velNormal(2) = velRelative(2) * normal(2);
    velNormal(3) = velRelative(3) * normal(3);
    velTangent(1) = velRelative(1) * (tangent1(1) + tangent2(1));
    velTangent(2) = velRelative(2) * (tangent1(2) + tangent2(2));
    velTangent(3) = velRelative(3) * (tangent1(3) + tangent2(3));
    wRelative(1) = w1(1) - w2(1);
    wRelative(2) = w1(2) - w2(2);
    wRelative(3) = w1(3) - w2(3);
    wTangent(1) = wRelative(1) * (tangent1(1) + tangent2(1));
    wTangent(2) = wRelative(2) * (tangent1(2) + tangent2(2));
    wTangent(3) = wRelative(3) * (tangent1(3) + tangent2(3));
    wNormal(1) = wRelative(1) * normal(1);
    wNormal(2) = wRelative(2) * normal(2);
    wNormal(3) = wRelative(3) * normal(3);
    velTangentMag = (velTangent(1) * velTangent(1) + velTangent(2) * velTangent(2) + velTangent(3) * velTangent(3)) ^ 0.5e0;
    wNormalMag = (wNormal(1) * wNormal(1) + wNormal(2) * wNormal(2) + wNormal(3) * wNormal(3)) ^ 0.5e0;
    wTangentMag = (wTangent(1) * wTangent(1) + wTangent(2) * wTangent(2) + wTangent(3) * wTangent(3)) ^ 0.5e0;
    interpenVolume = 0.1e1 / 0.3e1 * pieVal * pen * pen * (0.3e1 * R - pen);
    xc = 0.5e0 * pen;
    df = 0.0e0;
    f = interpenVolume;
    i = 0;
    while (0.1e-3 * interpenVolume < abs(f) & i < iterMax)
      f = 0.1e1 / 0.3e1 * pieVal * xc * xc * (0.3e1 * R - xc) - 0.5e0 * interpenVolume;
      df = pieVal * xc * (0.6e1 * R - 0.3e1 * xc) / 0.3e1;
      xc = xc - f / df;
      i = i + 1;
    end
    if (iterMax <= i)
      pen = 0.1e-2 * R;
      xc = 0.5e0 * pen;
    end
    interpenJ(1) = 0.1e1 / 0.30e2 * pieVal * pen * pen * pen * (0.3e1 * pen * pen - 0.15e2 * R * pen + 0.20e2 * R * R);
    interpenJ(2) = 0.1e1 / 0.60e2 * pieVal * pen * pen * (-0.9e1 * pen * pen * pen + 0.15e2 * R * pen * pen + 0.30e2 * pen * pen * xc - 0.20e2 * pen * xc * xc + 0.20e2 * R * R * pen - 0.80e2 * pen * R * xc + 0.60e2 * xc * xc * R);
    vcn(1) = -velNormal(1);
    vcn(2) = -velNormal(2);
    vcn(3) = -velNormal(3);
    if (c * (vcn(1) + vcn(2) + vcn(3)) < -0.1e1)
      contactForce(1) = 0.0e0;
      contactForce(2) = 0.0e0;
      contactForce(3) = 0.0e0;
    else
      contactForce(1) = k * interpenVolume * (0.1e1 + vcn(1) * c) * normal(1);
      contactForce(2) = k * interpenVolume * (0.1e1 + vcn(2) * c) * normal(2);
      contactForce(3) = k * interpenVolume * (0.1e1 + vcn(3) * c) * normal(3);
    end
    contactForceMag = (contactForce(1) * contactForce(1) + contactForce(2) * contactForce(2) + contactForce(3) * contactForce(3)) ^ 0.5e0;
    if (0.0e0 < wTangentMag)
      rollingResistance(1) = -k * c * interpenJ(2) * wTangent(1);
      rollingResistance(2) = -k * c * interpenJ(2) * wTangent(2);
      rollingResistance(3) = -k * c * interpenJ(2) * wTangent(3);
    else
      rollingResistance(1) = 0.0e0;
      rollingResistance(2) = 0.0e0;
      rollingResistance(3) = 0.0e0;
    end
    muVel = velTangentMag;
    if (0.0e0 < muVel)
      if (velD < muVel)
        muT = muD;
      elseif (velS < muVel & muVel < velD)
        delta = (muVel - velS) / (velD - velS);
        muT = muS + (muD - muS) * delta * delta * (0.30e1 - 0.20e1 * delta);
      else
        delta = (muVel + velS) / 0.20e1 / velS;
        muT = -muS + 0.2e1 * muS * delta * delta * (0.30e1 - 0.20e1 * delta);
      end
    end
    if (0.0e0 < velTangentMag)
      linearFrictionForce(1) = -muT * contactForceMag * velTangent(1);
      linearFrictionForce(2) = -muT * contactForceMag * velTangent(2);
      linearFrictionForce(3) = -muT * contactForceMag * velTangent(3);
    else
      linearFrictionForce(1) = 0.0e0;
      linearFrictionForce(2) = 0.0e0;
      linearFrictionForce(3) = 0.0e0;
    end
    linearFrictionMoment(1) = (normal(3) * linearFrictionForce(2) - normal(2) * linearFrictionForce(3)) * (R - pen);
    linearFrictionMoment(2) = (-normal(3) * linearFrictionForce(1) + normal(1) * linearFrictionForce(3)) * (R - pen);
    linearFrictionMoment(3) = (normal(2) * linearFrictionForce(1) - normal(1) * linearFrictionForce(2)) * (R - pen);
    radGy = sqrt(interpenJ(1) / interpenVolume);
    muVel = radGy * wNormalMag;
    if (0.0e0 < muVel)
      if (velD < muVel)
        muSpin = muD;
      elseif (velS < muVel & muVel < velD)
        delta = (muVel - velS) / (velD - velS);
        muSpin = muS + (muD - muS) * delta * delta * (0.30e1 - 0.20e1 * delta);
      else
        delta = (muVel + velS) / 0.20e1 / velS;
        muSpin = -muS + 0.2e1 * muS * delta * delta * (0.30e1 - 0.20e1 * delta);
      end
    end
    if (0.0e0 < wNormalMag)
      spinFrictionTorque(1) = -muSpin * contactForceMag * radGy * wNormal(1);
      spinFrictionTorque(2) = -muSpin * contactForceMag * radGy * wNormal(2);
      spinFrictionTorque(3) = -muSpin * contactForceMag * radGy * wNormal(3);
    else
      spinFrictionTorque(1) = 0.0e0;
      spinFrictionTorque(2) = 0.0e0;
      spinFrictionTorque(3) = 0.0e0;
    end
    aOut(1) = contactForce(1) + linearFrictionForce(1);
    aOut(2) = contactForce(2) + linearFrictionForce(2);
    aOut(3) = contactForce(3) + linearFrictionForce(3);
    aOut(4) = rollingResistance(1) + spinFrictionTorque(1) + linearFrictionMoment(1);
    aOut(5) = rollingResistance(2) + spinFrictionTorque(2) + linearFrictionMoment(2);
    aOut(6) = rollingResistance(3) + spinFrictionTorque(3) + linearFrictionMoment(3);
    aOut(7) = -aOut(1);
    aOut(8) = -aOut(2);
    aOut(9) = -aOut(3);
    aOut(10) = -aOut(4);
    aOut(11) = -aOut(5);
    aOut(12) = -aOut(6);
  end
  GonthierSCMreturn = 0.0e0;
