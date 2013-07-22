function err = calcContactForceError(x, ic0, vParams, vToe, expGRF)

Fx = 0;
Fy = 0;
Fz = 0;
Mx = 0;
My = 0;
Mz = 0;

err = 10e8;

%%
%Populate the initial conditions
%%
ic = ic0;
if(length(x) == 1) %height
    ic(10) = ic0(10) + x(1); 
end

if(length(x) == 2) %planar velocity
    ic(1) = ic0(1) + x(1);
    ic(2) = ic0(2) + x(2);
end

%%
% Compute passive toe torque
%%
idx_dth = 4;
idx_th  = 11;
kmt = vToe(1);
dmt = vToe(2);
TK1cK2a = -kmt*ic(idx_th)*(1+dmt*ic(idx_dth));

%%
%Compute the ground reaction force
%%
contactInfo=calcContactForcePosition(0,ic,vParams,[0 0 0 0 0 0 TK1cK2a]');
grfCOP = calcModelGRFCOP(contactInfo);
grf = grfCOP(1:3);
cop = grfCOP(4:6);

if(length(x)==1);
    err = (grf(3)-expGRF(3))^2;
end

if(length(x)==2)
    err = (grf(2)-expGRF(2))^2 + (grf(1)-expGRF(1))^2;    
end