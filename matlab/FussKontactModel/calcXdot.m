function xdot = calcXdot(t,x,vParams,vToe, expData, flag_mode)



Fx = 0;
Fy = 0;
Fz = 0;
Mx = 0;
My = 0;
Mz = 0;

%%
% Compute passive toe torque
%%
idx_dth = 4;
idx_th  = 11;
kmt = vToe(1);
dmt = vToe(2);
TK1cK2a = -kmt*x(idx_th) - dmt*x(idx_dth);

%%
%Populate input vector
%%

FxC = 0;
FyC = 0;
FzC = 0;
MxC = 0;
MyC = 0;
MzC = 0;   
vInputs = [FxC,FyC,FzC,MxC,MyC,MzC,TK1cK2a]';

if(flag_mode == 2)
    %%
    %Compute the ground reaction force
    %%
    contactInfo = calcContactForcePosition(t,x,vParams,[0 0 0 0 0 0 TK1cK2a]');
    grfCOP = calcModelGRFCOP(contactInfo);
    grf = grfCOP(1:3);
    cop = grfCOP(4:6);

    %%
    % Compute experimental model state at this time
    %%
    expX = zeros(14,1);
    expWRENCH = zeros(6,1);
    expGRF = zeros(3,1);
    expCOP = zeros(3,1);
    for i=1:1:14
        expX(i) = interp1(expData.time, expData.mdlState(:,i), t);
    end
    for i=1:1:3
       expGRF(i) = interp1(expData.time, expData.grf(:,i),t);
       expCOP(i) = interp1(expData.time, expData.cop(:,i),t);
    end
    for i=1:1:3
       expWRENCH(i) =  interp1(expData.time, expData.wrench(:,i),t);
    end

    
    %%
    % Tracking controller
    %%    
        
    
    errPos = expX(8:14) - (x(8:14)-expData.mdlStateOffset(8:14)');
    errVel = expX(1:7)- (x(1:7)-expData.mdlStateOffset(1:7)');
    
    
    %disp( sum(errPos.*errPos));
    
    kPosVel = [1e-1 0  0   0   0 0 0;...
             0 1e-1 0   0   0 0 0;...
             0 0  1e-1  0   0 0 0;...
             0 0  0   1e-3   0 0 0;...
             0 0  0   0   1e-3 0 0;...
             0 0  0   0   0 1e-3 0;...
             0 0  0   0   0 0 1e-1].*10; %In Mz entry: 1e-3
    
%     kVel = [1 0  0   0   0 0 0;...
%              0 1 0   0   0 0 0;...
%              0 0  1  0   0 0 0;...
%              0 0  0  1e-3   0 0 0;...
%              0 0  0   0   1e-1 0 0;...
%              0 0  0   0   0 1e-1 0;...
%              0 0  0   0   0 0 1e-1].*10;

    tmpPOS = kPosVel*errPos;
    tmpVEL = kPosVel*errVel;
    posTRACKING = [tmpPOS(1:3); tmpPOS(5:7)];
    velTRACKING = [tmpVEL(1:3); tmpVEL(5:7)];
    
    posToeTRACKING = [tmpPOS(4)];
    velToeTRACKING = [tmpVEL(4)];
    
    %%
    % Control
    %%
    errV = zeros(6,1);

    expR = expCOP-x(8:10);
    expM = cross(expGRF,expR);
    
    expWRENCH(4) = expM(1);
    expWRENCH(5) = expM(2);
    expWRENCH(6) = expM(3);
    
    simR = cop - x(8:10);    
    simM = cross(grf,simR);

    
    errV(1) = expGRF(1)-grf(1);
    errV(2) = expGRF(2)-grf(2);
    errV(3) = expGRF(3)-grf(3);
    errV(4) = expM(1)-simM(1);
    errV(5) = expM(2)-simM(2);
    errV(6) = expM(3)-simM(3);
    
    
    
    fbWRENCH = zeros(6,1);
    fbWRENCH(1) = -0.1*errV(1);    
    fbWRENCH(2) = -0.1*errV(2);    
    fbWRENCH(3) = -50*errV(3);     
    fbWRENCH(4) = 10*errV(4);
    fbWRENCH(5) = 10*errV(5);
    fbWRENCH(6) = 10*errV(6);
    
    
    %fbWRENCH = zeros(size(fbWRENCH));

    %disp(errV');

    %      vx   vy   vz   wx   wy   wz
    vel = [x(1) x(2) x(3) x(5) x(6) x(7)]';
    damp = [1 0 0 0 0 0;...
            0 1 0 0 0 0;...
            0 0 1 0 0 0;...
            0 0 0 5 0 0;...
            0 0 0 0 5 0;...
            0 0 0 0 0 5];

    dampWRENCH = -damp*vel.*(1);    

    %%
    %Enable or disable components
    %%

    s = (expGRF(3)/expData.grfMAXZ)*0.5;
    

    for i = 1:1:6
        if(i ~= 5)
        posTRACKING(i) = posTRACKING(i).*(1-s);
        velTRACKING(i) = velTRACKING(i).*(1-s);
        end
    end
    
    posToeTRACKING = posToeTRACKING.*(1-s);
    velToeTRACKING = velToeTRACKING.*(1-s);
    
    expWRENCH = expWRENCH.*(s+0.5);
    fbWRENCH = fbWRENCH.*(s+0.5);
    dampWRENCH = dampWRENCH.*(s+0.5);
    
    
    
    ctrl = expWRENCH+fbWRENCH+posTRACKING+velTRACKING+dampWRENCH;
    
    
    FxC = ctrl(1);
    FyC = ctrl(2);
    FzC = ctrl(3);

    MxC = ctrl(4);
    MyC = ctrl(5);
    MzC = ctrl(6);

    Mtoe = posToeTRACKING + velToeTRACKING;
    
    vInputs = [FxC,FyC,FzC,MxC,MyC,MzC,TK1cK2a]';
end


tmp = xDotMex(t, x, vParams, vInputs);
xdot = tmp(1:1:length(x));