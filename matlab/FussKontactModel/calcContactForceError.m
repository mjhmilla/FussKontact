function err = calcContactForceError(x, ic0, vParams, vToe, expGRF, expCOP)

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

if(length(x) == 1)
    ic(10) = ic0(10) + x(1)./1000; %z
end

if(length(x) == 3)
    ic(10) = ic0(10) + x(1)./1000; %z
    ic(12) = ic0(12) + x(2)./1000; %zeta
    ic(13) = ic0(13) + x(3)./1000; %eta
end


%%
%Compute the ground reaction force
%%
contactInfo=calcContactForcePosition(0,ic,vParams,[0 0 0 0 0 0 0]');
grfCOP = calcModelGRFCOP(contactInfo);
grf = grfCOP(1:3);
cop = grfCOP(4:6);


errGRF = (grf(3)-expGRF(3))^2;
errCOPV = (cop(1:2)' - expCOP(1:2));
errCOP = sum(errCOPV.*errCOPV);
err = errGRF;

if(length(x) == 3)
    err = errGRF + 1000000*errCOP;
end
