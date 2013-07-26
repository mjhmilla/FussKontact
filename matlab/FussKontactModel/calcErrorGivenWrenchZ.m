function err = calcErrorGivenWrenchZ(xDelta, t, scaling, xPose0 ,vParams, vToe, expData)

err = inf;

%%
%Get the experimental data at time t
%%
expX = zeros(14,1);
for i=1:1:14
    expX(i) = interp1(expData.time, expData.mdlState(:,i), t);
end
%%
%Construct the static pose
%%

x = zeros(14,1);
for i=1:1:7
   x(7+i) = xPose0(i) + xDelta(i)./scaling; 
end

%%
% Compute passive toe torque
%%
TK1cK2a = calcToeTorque(x,vToe);

%%
% Apply the experimental wrench to the foot
%%

vInputs = [0,0,0,0,0,0,TK1cK2a]';

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
