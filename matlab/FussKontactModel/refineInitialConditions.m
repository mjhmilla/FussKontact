function ic = refineInitialConditions(time0, vIC, expData, vParams, vToe)

ic = vIC;

expIdx0 = min(find(expData.time >= time0));
expGRF = expData.grf;
expCOP = expData.cop;

%Start the foot at rest;
ic(1) = 0;
ic(2) = 0;
ic(3) = 0;
ic(4) = 0;
ic(5) = 0;
ic(6) = 0;
ic(7) = 0;

scaling = 10000;

%%
%Compute z, such that Fz matche experimental data
%%

%%
% Put the foot in the neighborhood of a solution
%%

%1. If both the fore and rear foot are in contact, level the foot out.
ic(11) = 0;
flag_FootFlat = 0;

contactInfo=calcContactForcePosition(0,ic,vParams,[0 0 0 0 0 0 0]');
zHeel = contactInfo(6);
zFore = contactInfo(12);   

idxZ = 10;
idxTh= 11;

if( abs(zHeel-zFore) < 0.02)
    errGapFcn = @(x)calcGapErr(x, ic, vParams);
    [zetaDelta gapErr] = fminunc(errGapFcn,[0]);
    ic(12) = ic(12) + zetaDelta/1000;       
    flag_FootFlat = 1;
end

options = optimset('Display','off','Diagnostics','off',...
                   'TolX',1e-12,'TolFun',1e-12, 'MaxIter', 2000);

if(zHeel - zFore > 0.02 && flag_FootFlat == 0)
    ic(10) = ic(10) - zFore + 0.1; %Put the fore foot in the air
    
    %Align the ankle so the the metarsal joint aligns with the
    %corresponding position on the experimental food
    idxM1D = 0;
    idxM5D = 0;
    for j=1:1:length(expData.mkrNames)
       if( isempty(strfind(expData.mkrNames{j},'M1D'))~=1)
          idxM1D = j; 
       end
       if( isempty(strfind(expData.mkrNames{j},'M5D'))~=1)
          idxM5D = j; 
       end
    end
    
    rM1D = [expData.mkrs(idxM1D,1,expIdx0) ...
        expData.mkrs(idxM1D,2,expIdx0) ...
        expData.mkrs(idxM1D,3,expIdx0)];
    rM5D = [expData.mkrs(idxM5D,1,expIdx0) ...
            expData.mkrs(idxM5D,2,expIdx0) ...
            expData.mkrs(idxM5D,3,expIdx0)];
    rM15D = 0.5.*(rM1D + rM5D);
    
    kin = calcFootFrameKinematics(ic,vParams,vToe);
        
    rMTJ  = kin.Fore.r;
    
    deltaXY = [(rM15D(1)-rMTJ(1)) (rM15D(2)-rMTJ(2)) 0];
    ic(8) = ic(8) + deltaXY(1);
    ic(9) = ic(9) + deltaXY(2);
    
    %Change the toe angle to make the MT pad flat with the ground
    errCFcn = @(x)calcContactForceError(x, [idxTh] ,ic, scaling, ...
                                    vToe,vParams, vIC, ...
                                    expGRF(expIdx0,:), expCOP(expIdx0,:));                                
    [icDelta FzErrSq] = fminunc(errCFcn,[0], options);
    ic(idxTh) = ic(idxTh) + icDelta(1)/scaling; %z        
end

%2. Now make small changes to the height to match Fz
contactInfo=calcContactForcePosition(0,ic,vParams,[0 0 0 0 0 0 0]');
gap = min(contactInfo(6), contactInfo(12));   
ic(10) = ic(10)-gap-0.005;   

errCFcn = @(x)calcContactForceError(x, [idxZ] ,ic, scaling, ...
                                    vToe,vParams, vIC, ...
                                    expGRF(expIdx0,:), expCOP(expIdx0,:));                                
[icDelta FzErrSq] = fminunc(errCFcn,[0], options);
ic(idxZ) = ic(idxZ) + icDelta(1)/scaling; %z


if flag_FootFlat == 1
    %Get the COP correct
    errCFcn = @(x)calcContactForceError(x, [12,13],ic, scaling, ...
                                     vToe,vParams, vIC,...
                                     expGRF(expIdx0,:), expCOP(expIdx0,:));
    [icDelta copErrSq] = fminunc(errCFcn,[0 0], options);

    ic(12) = ic(12) + icDelta(1)./scaling; %zeta
    ic(13) = ic(13) + icDelta(2)./scaling; %eta
end

%%
%Output results
%%
contactInfo=calcContactForcePosition(0,ic,vParams,[0 0 0 0 0 0 0]');
grfCOP = calcModelGRFCOP(contactInfo);
grf = grfCOP(1:3);
cop = grfCOP(4:6);
disp(sprintf('Mdl:Exp FxErr(%f) FyErr(%f) FzErr(%f)',...
    grf(1)-expGRF(expIdx0,1),grf(2)-expGRF(expIdx0,2),...
    grf(3)-expGRF(expIdx0,3)));
    
disp(sprintf('        COPxErr(%f) COPyErr(%f)',...
    cop(1)-expCOP(expIdx0,1),cop(2)-expCOP(expIdx0,2)));
