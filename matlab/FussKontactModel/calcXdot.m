function xdot = calcXdot(t,x,vParams,vToe)

Fx = 0;
Fy = 0;
Fz = 0;
Mx = 0;
My = 0;
Mz = 0;

idx_dth = 4;
idx_th  = 11;

kmt = vToe(1);
dmt = vToe(2);

TK1cK2a = -kmt*x(idx_th)*(1+dmt*x(idx_dth));

vInputs = [Fx,Fy,Fz,Mx,My,Mz,TK1cK2a]';

tmp = xDotMex(t, x, vParams, vInputs);
xdot = tmp(1:1:length(x));