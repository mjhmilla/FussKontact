clc;
close all;
clear all;


tsim = [0 1];
aniFreq = 100;
%%
% Setup parameters
%%
disp('1. Preprocessing');
getParams; %vParams
getInputs; %vInputs
getIC;     %vIC
vToe = [50, 0.5]; %Stiffness & damping at the nonlinear toe joint
t0 = 0;


%%
%Setup for a simulation
%%

%anonomous function
xdotAFunc = @(t,x) calcXdot(t,x,vParams,vToe);

%%
% Simulate
%%
disp('2. Simulation');
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

ticID = tic;
[t sol] = ode15s(xdotAFunc,tsim,vIC,options);
elapsedTime = toc(ticID);
disp(sprintf('%f s of simulation in %f s',tsim(2), elapsedTime));

contactInfo = zeros(length(t),12);
for i=1:1:length(t)
   contactInfo(i,:) = calcContactForcePosition(t(i),sol(i,:)',vParams,vInputs)'; 
end

inputInfo = zeros(length(t),7);

%%
% Post process
%%
disp('3. Postprocessing');
postprocess(t,sol,contactInfo,inputInfo,aniFreq);
