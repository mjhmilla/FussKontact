function err = calcFootError(x, xscaling, optParamIdx, vParam,vToe, ...
                             errorScaling, errorTimeScaling, ...
                             expFileNames, expCTRLGAINS, ...
                             expTimeRange, resultsFolder,...                             
                             flag_animate)
%%
%
% @param x: optimization variables
% @param xscaling: scaling applied to each variable in x
% @param optParamIdx: A nx2 vector. The first column specifies the
%                     parameter list (1 is vParams, 2 is vToe) and the
%                     second specifies the index of the variable being
%                     optimized. The order of indices listed must
%                     correspond to the order of the optimization variables
%                     in x
% @param vParam: the parameter array for this model
% @param vToe: the parameters associated with the toe joint
% @param errorScaling: Scaling applied to the squared position error 
%                      between the ankle on the x,y,z, zeta, eta, xi 
%                      dimensions
% @param errorTimeScaling: Scaling associated with an incomplete simulation
%                          this happens because the kinematic error between
%                          the simulation and the model becomes very large
% @param expFileNames: A cell array containing the names of the processed
%                      experimental data files
% @param expCTRLGAINS: A n x 3 matrix of control gains. Each row
%                      corresponds to the control gains for the i^th 
%                      expFileName
%
% @return err: the kinematic error between the model's ankle frame and the
%              experimental model frame. This error is the squared
%              integrated sum of the difference scaled by the duration of
%              the movement. The movements used to test the foot are
%              walking, running, and a large rotation
%%

err = inf;

%%
% Update vParam and vToe with x
%%

for i=1:1:length(x)
   if(optParamIdx(i,1) == 1)
      pIdx = optParamIdx(i,2);
      vParam(pIdx) = vParam(pIdx) + x(i)/xscaling(i);
   end
   
   if(optParamIdx(i,1) == 2)
       pIdx = optParamIdx(i,2);
       vToe(pIdx) = vToe(pIdx) + x(i)/xscaling(i);
   end    
end

errKIN = zeros(length(expFileNames),length(errorScaling));
errFINISHTIME = zeros(length(expFileNames),1);

for i = 1:1:length(expFileNames)
    load(expFileNames{i});    
    

    expt0 = min(expData.time);       
    expt1 = max(expData.time);
    tDelta = expt1-expt0;
    t0 = expt0 + tDelta*expTimeRange(i,1);
    t1 = expt0 + tDelta*expTimeRange(i,2);
    
    
    tsim = [t0 t1];         
    idx0 = min(find(expData.time >= t0));
    vIC = expData.mdlState(idx0,:)';          


    %2. Refine IC to match the GRF and COP
    vIC = refineInitialConditions(t0, vIC, expData, vParam, vToe);
    expData.mdlStateOffset = vIC'-expData.mdlState(1,:);

    ctrlGAIN = expCTRLGAINS(i,:);
    
    %3. Simulate the model
    flag_mode = 2;
    xdotAFunc = @(targ,xarg) calcXdot(targ,xarg,vParam,...
                                      vToe, expData, ctrlGAIN, flag_mode);
    xdotEvent = @(targ,yarg) footEvent(targ,yarg,expData);

    
    
    options = [];
    if flag_animate == 0
        options = odeset('RelTol',1e-4,'AbsTol',1e-5,'Events',xdotEvent);
    else
        options = odeset('RelTol',1e-4,'AbsTol',1e-5,'Events',xdotEvent,...
                         'OutputFcn', @getIntegratorOutput);
    end
    
    [t sol] = ode15s(xdotAFunc,tsim,vIC,options);

    %4. Resample the solution using the sense of time in expData
    t1 = max(find(expData.time <= max(t)));
    errTime = expData.time(1:1:t1);
    solSampled = zeros(length(errTime), size(sol,2));

    for j=1:1:size(sol,2)
        solSampled(:,j) = interp1(t, sol(:,j),errTime);
    end

    %5. Compute the integral of the squared difference between the simulated
    %   ankle kinematics and the experimental kinematics

    xError = zeros(length(errTime), 6);
    xDist  = zeros(length(errTime),6);
    idxError = [8 9 10 12 13 14]'; %x,y,z,zeta,eta,xi

    for j=1:1:length(idxError)
        idx = idxError(j);
        errV = (solSampled(:,idx)-expData.mdlState(1:t1,idx)).^2;
        xError(:,j) = cumtrapz(errTime, errV)./max(errTime);
        xDist(:,j)  = errV.^0.5;
        errKIN(i,j) = xError(length(errTime),j)*errorScaling(j);
    end

    errFINISHTIME(i) = ( (max(expData.time)-max(errTime))/max(expData.time) )^2;
   
    if(flag_animate == 1)
       contactInfo = zeros(length(t),12);
       inputInfo = zeros(length(t),7); 
       aniPoints = ceil(100*(max(t)-min(t)));
       for j=1:1:length(t)
        contactInfo(j,:) = calcContactForcePosition(t(j),sol(j,:)',...
                                       vParam,[0 0 0 0 0 0 0]')'; 
       end       
            
        fig = figure;
            mdlGRFCOP = zeros(length(t),6);
            for j=1:1:length(t)
                mdlGRFCOP(j,:) = calcModelGRFCOP(contactInfo(j,:)');   
            end
            subplot(1,2,1);
              plot(expData.time,expData.grf);
              hold on;
              plot(t, [mdlGRFCOP(:,1) mdlGRFCOP(:,2) mdlGRFCOP(:,3)],'--');
              title([expFileNames{i}, ' GRF']);
            subplot(1,2,2);
              title([expFileNames{i}, ' COP']);
              plot(expData.cop(:,1), expData.cop(:,2),'b');
              hold on;
              plot(mdlGRFCOP(:,4), mdlGRFCOP(:,5),'r');
              axis equal;

       
       
       
       postprocess(t,sol,contactInfo,inputInfo,aniPoints,expData,...
                   vParam,vToe);
       
       rootDir = pwd;
       
       outDir = rootDir;       
       cd('animation');
           outDir = [outDir,'\',resultsFolder,'\'];
           fname0 = expFileNames{i};
           idx0 = strfind(fname0,'\');
           idx1 = strfind(fname0,'.');
           fname1 = fname0((idx0+1):1:(idx1-1));
           fname2 = ['opt_',fname1,'.wrl'];
       dos(['copy footMdlExp.wrl ', outDir, fname2]);
       cd(rootDir);
              
       output.exp.time = expData.time;
       output.exp.grf  = expData.grf;
       output.exp.cop  = expData.cop;
       output.exp.kin  = expData.mdlState;
       output.mdl.time = t;
       output.mdl.grf  = [mdlGRFCOP(:,1) mdlGRFCOP(:,2) mdlGRFCOP(:,3)];
       output.mdl.cop  = [mdlGRFCOP(:,4) mdlGRFCOP(:,5) mdlGRFCOP(:,6)];
       output.mdl.kin  = sol;
       output.err.time = errTime;
       output.err.x    = xError;
       output.err.xDist     = xDist; 
       output.err.xScaling  = errorScaling;
       output.err.kin       = errKIN(i,:);
       output.err.finishTime= errFINISHTIME(i);
       output.err.finishTimeScaling = errorTimeScaling;
       
       cd(resultsFolder)
        save(['opt_',fname1,'.mat'],'output');
       cd(rootDir);
       
    end
    
end

err = sum(sum(errKIN)) + errorTimeScaling*sum(errFINISHTIME);


                