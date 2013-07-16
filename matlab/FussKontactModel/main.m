clc;
close all;
clear all;


tsim = [0 1];

%%
% Setup parameters
%%
getParams; %vParams
getInputs; %vInputs
getIC;     %vIC
t0 = 0;

disp('Initial Test');
xdot = xDotMex(t0, vIC, vParams, vInputs);
disp(xdot);

%%
%Setup for a simulation
%%

%anonomous function
xdotAFunc = @(t,x) calcXdot(t,x,vParams,vInputs);

%%
% Simulate
%%
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t sol] = ode45(xdotAFunc,tsim,vIC,options);

contactInfo = zeros(length(t),12);
for i=1:1:length(t)
   contactInfo(i,:) = calcContactForcePosition(t(i),sol(i,:)',vParams,vInputs)'; 
end

inputInfo = zeros(length(t),7);

%%
% Post process
%%

