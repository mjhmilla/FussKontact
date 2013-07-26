function mtTq = calcToeTorque(x,vToe)

stateIdx;

th = x(11);
dth= x(4);

kmt = vToe(1);
dmt = vToe(2);
thOff = vToe(3);
ang = th-thOff;
mtTq = -kmt*sign(ang)*(ang)^2 - dmt*dth;

