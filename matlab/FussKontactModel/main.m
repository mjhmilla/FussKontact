clc;
close all;
clear all;


tsim = [0 2];
aniPoints = 100*tsim(2);
expFile = 'Walking01.mat';
grfFiltFreq = [];
tmax = inf;

flag_refineIC = 1;
flag_mode = 2;
% 0: go through fitted foot kinematics 
% 1: forward simulation
% 2: forward simulation with controller to match grf
% 3: optimization a la Prof. Kecskemethy's IUTAM paper



ctrlGAIN = [100, 10, 1];
%Running: [100,10,1];
%Walking: [10, 2, 1];
%Rotations: [10 10 0];
%%
% Setup the model
%%
disp('1. Preprocessing');
[vParams vToe]= getParams();
vInputs = getInputs(); 
vIC = getIC(); %    


t0 = 0;

%%
% Compute the initial conditions against experimental data such that 
%  Exact Match: GRF & COP
%  Similiar   : Foot state
%%

expData = calcExpDataInModelCoord(['data/',expFile], vToe,grfFiltFreq);

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
           vIC = refineInitialConditions(vIC, expData, vParams, vToe);

           expData.mdlStateOffset = vIC'-expData.mdlState(1,:);
       end
   end   
end



save(['modeldata/',expFile], 'expData');


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
        
        xdotAFunc = @(targ,xarg) calcXdot(targ,xarg,vParams,...
                                          vToe, expData, ctrlGAIN, flag_mode);
        xdotEvent = @(targ,yarg) footEvent(targ,yarg,expData);
        
        options = odeset('RelTol',1e-4,'AbsTol',1e-5,...
                         'OutputFcn',@getIntegratorOutput,'Events',xdotEvent);
        
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

        %Find z, th, zeta, eta, xi s.t. Fz        
        x0 = zeros(size(vIC));
        x0(8:14) = vIC(8:14);
        xDelta = zeros(5,1);
        
        
        scaling = 10000;
        options = optimset('Display','off','Diagnostics','off',...
                           'TolX',1e-12, 'TolFun', 1e-12,'MaxIter',1000);
           
        expX0 = zeros(size(x0));
                       
        for i = 1:1:length(expData.time)
           
           %Update the previous solution for this time
           x0(1:7)   = 0;
           x0(8:9)   = expData.mdlState(i,8:9);
           x0(10:14) = x0(10:14) + xDelta(1:5)./scaling;
             
           expX0(8:14) = [x0(8:11); expData.mdlState(i,12:14)'];
           
           %Update the anomomous function 
           errCFcn = @(x)calcContactForceError(x, x0, scaling,vParams,...
                expX0, expData.grf(i,:), expData.cop(i,:));

           %Optimize
           ticID = tic;
           [xDelta err exitflag] = fminunc(errCFcn,zeros(5,1), options);
           elapsedTime = toc(ticID);
           
           %Update the solution given the optimized solution
           x0(1:7)   = 0;
           x0(8:9)   = expData.mdlState(i,8:9);
           x0(10:14) = x0(10:14) + xDelta(1:5)./scaling;
           
           %
           contactInfo=calcContactForcePosition(0,x0,vParams,[0 0 0 0 0 0 0]');
           grfCOP = calcModelGRFCOP(contactInfo);
           grf = grfCOP(1:3);
           cop = grfCOP(4:6);
                                
           %Output
           disp(sprintf('Iter %i/%i: %f s, dFz %f dCOPx %f dCOPy %f',...
                        i,length(expData.time),elapsedTime,...
                        (expData.grf(i,3)-grf(3)),...
                        (expData.cop(i,1)-cop(1)),...
                        (expData.cop(i,2)-cop(2))));
           
           %Store the solution
           solPose(i,:) = [x0(8:14)'];           
           solErr(i) = err;
           solExitFlag(i) = exitflag;

        end
        
        t = expData.time;
        inputInfo = zeros(length(t),7);
        sol = [zeros(size(solPose)) solPose];
        
end


contactInfo = zeros(length(t),12);
for i=1:1:length(t)
   contactInfo(i,:) = calcContactForcePosition(t(i),sol(i,:)',vParams,vInputs)'; 
   if(sum(isnan(contactInfo(i,:)))>0)
       contactInfo(i,:) = calcContactForcePosition(t(i),sol(i,:)',vParams,vInputs)';
   end       
end

if(flag_mode == 2 || flag_mode == 3)
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
