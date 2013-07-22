clc;
close all;
clear all;


tsim = [0 2];
aniPoints = 100*tsim(2);
expFile = 'data/Rotations10.mat';
grfFiltFreq = [];
tmax = 10;

flag_refineIC = 1;
flag_mode = 2;
% 0: go through fitted foot kinematics 
% 1: forward simulation
% 2: forward simulation with controller to match grf
% 3: optimization a la Prof. Kecskemethy's IUTAM paper


vToe = [5, 5]; %Stiffness & damping at the nonlinear toe joint

%%
% Setup the model
%%
disp('1. Preprocessing');
vParams = getParams();
vInputs = getInputs(); 
vIC = getIC(); %    

t0 = 0;

%%
% Compute the initial conditions against experimental data such that 
%  Exact Match: GRF & COP
%  Similiar   : Foot state
%%

expData = calcExpDataInModelCoord(expFile, grfFiltFreq);
if(isempty(expData)~=1)           
   if(flag_mode >= 2)
       t0 = min(expData.time);
       t1 = max(expData.time);
       if(t1 > tmax);
        t1 = tmax;
       end
       
       tsim = [t0 t1];
       
       
       vIC = expData.mdlState(1,:)';          

       if flag_refineIC == 1
           vIC = refineInitialConditions(vIC, expData.grf, expData.cop,...
                                         vParams,vInputs,vToe);

           expData.mdlStateOffset(8:10) = vIC(8:10)'-expData.mdlState(1,8:10);
       else
           
       end
   end   
end




%%
% Simulate/Pose
%%

t = [];
sol = [];
contactInfo = [];
inputInfo = [];


switch(flag_mode)
    case 0
        disp('0: Experimental poses');
        t = expData.time;
        sol = expData.mdlState;
        contactInfo = zeros(length(t),12);
        for i=1:1:length(t)
           contactInfo(i,:) = calcContactForcePosition(...
                            t(i),sol(i,:)',vParams,vInputs)';        
        end        
        inputInfo = zeros(length(t),7);  

    case 1
        disp('1: Forward simulation');  
        xdotAFunc = @(targ,xarg) calcXdot(targ,xarg,vParams,vToe, expData, flag_mode);        
        options = odeset('RelTol',1e-4,'AbsTol',1e-4);
        
        ticID = tic;
        [t sol] = ode45(xdotAFunc,tsim,vIC,options);
        elapsedTime = toc(ticID);
        disp(sprintf('%f s of simulation in %f s',tsim(2), elapsedTime));
        inputInfo = zeros(length(t),7);

    case 2
        disp('2: Forward simulation with controller');        
        
        disp('  Edit later: Major functionality from controller missing');
        xdotAFunc = @(targ,xarg) calcXdot(targ,xarg,vParams,vToe, expData, flag_mode);
        xdotEvent = @(targ,yarg) footEvent(targ,yarg,vParams,expData);
        
        options = odeset('RelTol',1e-4,'AbsTol',1e-5,...
                         'OutputFcn',@getIntegratorOutput);%,'Events',xdotEvent);
        
        ticID = tic;
        [t sol] = ode15s(xdotAFunc,tsim,vIC,options);
        elapsedTime = toc(ticID);
        disp(sprintf('%f s of simulation in %f s',tsim(2), elapsedTime));
        inputInfo = zeros(length(t),7);
                
        
    case 3
        disp('3. Minimize accelerations given applied wrench');
        t = 0;
        solPose = zeros(length(expData.time),7);
        solErr  = zeros(length(expData.time),1);        
        solExitFlag = zeros(length(expData.time),1);
        inputInfo = zeros(length(expData.time),7);
        %Find z, zeta, eta, xi s.t. Fz and its wrench are reproduced
        x0 = [vIC(10), vIC(11), vIC(12)];%, vIC(13), vIC(14)]';
        x0 = double(x0);
        options = optimset('Display','off','Diagnostics','off',...
                           'TolX',1e-12, 'TolFun', 1e-12,'MaxIter',1000);
        for i = 1:1:length(expData.time)
           t = expData.time(i);
           errWrenchZ = @(x) calcErrorGivenWrenchZ(x,t,vParams,vToe,expData);
           
           ticID = tic;
           [sol fval exitflag output grad] = fminunc(errWrenchZ,x0,options);
           elapsedTime = toc(ticID);
           disp(sprintf('Iter %i/%i: %f s',...
               i,length(expData.time),elapsedTime));
           
           x0 = sol;
           err = errWrenchZ(sol);
           
           expX = expData.mdlState(i,:);
           solPose(i,:) = [expX(8:14)];
           solPose(i,3:5)= sol;           
           
           solErr(i) = fval;
           solExitFlag(i) = exitflag;
           
           expWrench = expData.wrenchZ(i,:);
           inputInfo(i,:) = [expWrench, 0];
           
        end
        
        t = expData.time;
        sol = [zeros(size(solPose)) solPose];
        
end


contactInfo = zeros(length(t),12);
for i=1:1:length(t)
   contactInfo(i,:) = calcContactForcePosition(t(i),sol(i,:)',vParams,vInputs)'; 
   if(sum(isnan(contactInfo(i,:)))>0)
       contactInfo(i,:) = calcContactForcePosition(t(i),sol(i,:)',vParams,vInputs)';
   end       
end

if(flag_mode == 2)
   fig = figure;
   mdlGRFCOP = zeros(length(t),6);
   for i=1:1:length(t)
    mdlGRFCOP(i,:) = calcModelGRFCOP(contactInfo(i,:)');   
   end
   
   subplot(1,2,1);
    plot(expData.time,expData.grf);
    hold on;
    plot(t, [mdlGRFCOP(:,1) mdlGRFCOP(:,2) mdlGRFCOP(:,3)],'--');
   subplot(1,2,2);
    plot(expData.cop(:,1), expData.cop(:,2),'b');
   hold on;
    plot(mdlGRFCOP(:,4), mdlGRFCOP(:,5),'r');
   axis equal;
end

%%
% Post process
%%
disp('3. Postprocessing');
mIdx = inf;
for i=1:1:14
   tmp = min(find(isnan(sol(:,i))));
   if(tmp < mIdx-1)
      mIdx = tmp-1;       
   end
end

if(isinf(mIdx) == 1)
   mIdx = length(t); 
else
   disp('NaN in solution, solution trimmed'); 
end


postprocess(t(1:mIdx),sol(1:mIdx,:),contactInfo(1:mIdx,:),inputInfo(1:mIdx,:),aniPoints,expData);
