clc;
close all;
clear all;


tsim = [0 0.01];
aniFreq = 100;
expFile = 'data/Walking01.mat';

flag_simulate = 0;


%%
% Setup parameters
%%
disp('1. Preprocessing');
vParams = getParams();
vInputs = getInputs(); 
vIC = getIC(); %    
vToe = [100, 0.5]; %Stiffness & damping at the nonlinear toe joint
t0 = 0;

expData = calcExpDataInModelCoord(expFile);
if(isempty(expData)~=1)
   vIC = expData.mdlState(1,:)'; 
end

%anonomous function
xdotAFunc = @(t,x) calcXdot(t,x,vParams,vToe);

%%
% Simulate/Pose
%%

t = [];
sol = [];
contactInfo = [];
inputInfo = [];

if(flag_simulate == 1)
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
else
    t = expData.time;
    sol = expData.mdlState;
    contactInfo = zeros(length(t),12);
    inputInfo = zeros(length(t),7);    
end


%%
% Post process
%%
disp('3. Postprocessing');
postprocess(t,sol,contactInfo,inputInfo,aniFreq,expData);
