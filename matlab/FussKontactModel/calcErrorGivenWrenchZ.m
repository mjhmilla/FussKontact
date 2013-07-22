function err = calcErrorGivenWrenchZ(xPose, t, vParams, vToe, expData)

err = inf;

%%
%Get the experimental data at time t
%%
wrenchZ = zeros(6,1);
expX = zeros(14,1);
for i=1:1:6
    wrenchZ(i) = interp1(expData.time, expData.wrenchZ(:,i),t);
end
for i=1:1:14
    expX(i) = interp1(expData.time, expData.mdlState(:,i), t);
end
%%
%Construct the static pose
%%
x = zeros(14,1);
x(8) = expX(8);
x(9) = expX(9);
x(10) = xPose(1);
x(11) = xPose(2);
x(12) = xPose(3);
x(13) = expX(13);%xPose(4);
x(14) = expX(14);%xPose(5);

%%
% Compute passive toe torque
%%
idx_dth = 4;
idx_th  = 11;
kmt = vToe(1);
dmt = vToe(2);
TK1cK2a = -kmt*x(idx_th) - dmt*x(idx_dth);

%%
% Apply the experimental wrench to the foot
%%
FxC = wrenchZ(1);
FyC = wrenchZ(2);
FzC = wrenchZ(3);
MxC = wrenchZ(4);
MyC = wrenchZ(5);
MzC = wrenchZ(6);

vInputs = [FxC,FyC,FzC,MxC,MyC,MzC,TK1cK2a]';

%%
% Compute the accelerations of the foot
%%
tmp = xDotMex(t, x, vParams, vInputs);
xdot = tmp(1:1:length(x));

contactInfo = tmp(length(x)+1:length(tmp));
grfcop = calcModelGRFCOP(contactInfo);

grf = grfcop(1:3);
grfM = sum(grf.*grf)^0.5;
cop = grfcop(4:6);
gapErr = 0;

if(grfM <= eps)
    gapErr = cop(3).*cop(3);
end

%%
% Define the error as the the squared sum of all of the accelerations
% with a term to guide the foot to the ground
%%
xdotSq = sum(xdot(1:7).*xdot(1:7));
xSq = expX(8:14)-x(8:14);
xSq = sum(xSq.*xSq);

err = xdotSq + 100000000*xSq + 1000000*gapErr;
