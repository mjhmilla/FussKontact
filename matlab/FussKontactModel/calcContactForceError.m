function err = calcContactForceError(var, icIdx, ic0, scaling, vToe, vParams, expX, expGRF, expCOP)

err = 10e8;

Fx = 0;
Fy = 0;
Fz = 0;
Mx = 0;
My = 0;
Mz = 0;
TK1cK2a = 0;

%%
%Populate the initial conditions
%%
ic = ic0;

for i=1:1:length(var)
   ic(icIdx(i)) = ic0(icIdx(i)) + var(i)/scaling; 
end

%%
% Compute passive toe torque
%%

TK1cK2a = calcToeTorque(ic,vToe);



%%
%Compute the GRF and COP error terms
%%
tmp = xDotMex(0, ic, vParams, [0 0 0 0 0 0 TK1cK2a]');
xdot = tmp(1:1:length(ic));
contactInfo = tmp(length(ic)+1:1:length(tmp));

%%
% Compute location and rotation of every relevant frame
%%

%Parameters
buildParamVariableList;

%State Information
x    = ic(8);
y    = ic(9);
z    = ic(10);
th   = ic(11);
zeta = ic(12);
eta  = ic(13);
xi   = ic(14);

%Contact Information
HFx = contactInfo(1);
HFy = contactInfo(2);
HFz = contactInfo(3);
Hcopx  = contactInfo(4);
Hcopy  = contactInfo(5);
Hcopz  = contactInfo(6);

FFx = contactInfo(7);
FFy = contactInfo(8);
FFz = contactInfo(9);
Fcopx  = contactInfo(10);
Fcopy  = contactInfo(11);
Fcopz  = contactInfo(12);
    
%Compute all kinematic positions and frame locations.
calcPosVecsRotMatrices;


%%
% Compute error functions
%%
grfCOP = calcModelGRFCOP(contactInfo);
grf = grfCOP(1:3);
cop = grfCOP(4:6);


errGRFZ = (grf(3)-expGRF(3))^2;
errCOPV = (cop(1:2)' - expCOP(1:2));
errCOP = sum(errCOPV.*errCOPV);

stateIdx;
errddTH = xdot(idx_dth)^2;

errForeFlatV = rForeFoot - [Fcopx Fcopy Fcopz];
errForeFlat  = sum(errForeFlatV.*errForeFlatV);

%%
% Compute the total error 
%%

if(length(var) == 1)
    if(icIdx(1) == 10)
        err = errGRFZ;                
    end
    if(icIdx(1) == 11)
        err = 1e6*errForeFlat; 
    end
end

if(length(var) == 2)    
    copWeight = 0;
    if(icIdx(1) == 12 && icIdx(2) == 13)
       copWeight = 1e7; 
    end
    err = errGRFZ + copWeight*errCOP;    
end





