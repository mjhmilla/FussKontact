function xdot = calcXdot(t,x,vParams,vToe, expData, ctrlWeights, flag_mode)
%%
% Computes the state derivative of the foot contact model. If flag_mode is 
% 2, then control forces are computed to drive the model foot to have the
% same ground reaction force vector and COP position as was observed in
% the experimental data
% 
% @param t: time
% @param x: the state vector (vx,vy,vz,dth,wx,wy,wz,x,y,z,th,zeta,eta,xi)
% @param vToe: the stiffness and damping (a 1x2 vector) of the toe joint
% @param expData: a structure that contains experimental data, as well
%                 as very specific transformations of that data for this
%                 foot contact model.
% @param ctrlWeights: relative weights applied to the position and force
%                     feedback controller that is used when flag_mode=2.
% @param flag_mode: 0 the foot passively falls. If this parameter is 2
%                   then a control wrench is applied to the foot to
%                   drive it to reproduce the observed ground reaction 
%                   force and cop trajectory
%
%%

FxC = 0;
FyC = 0;
FzC = 0;
MxC = 0;
MyC = 0;
MzC = 0;

%%
% Compute passive toe torque
%%
TK1cK2a = calcToeTorque(x,vToe);

%%
%Populate input vector
%%
vInputs = [FxC,FyC,FzC,MxC,MyC,MzC,TK1cK2a]';

if(flag_mode == 2)
    %%
    % Get the control weights
    %%
    expXGain = ctrlWeights(1);
    dampGain = ctrlWeights(2);
    smax     = ctrlWeights(3);
    
    %%
    %Compute the ground reaction force
    %%
    contactInfo = calcContactForcePosition(t,x,vParams(1:50),[0 0 0 0 0 0 TK1cK2a]');
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
    
    kPosVel = [1e-1 0  0   0    0    0    0;...
               0 1e-1  0   0    0    0    0;...
               0 0  1e-1   0    0    0    0;...
               0 0  0      1e-3 0    0    0;...
               0 0  0      0    1e-3    0    0;...
               0 0  0      0    0    1e-3    0;...
               0 0  0      0    0    0    1e-3].*expXGain; %In Mz entry: 1e-3
    
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
    fbWRENCH(1) = 0*errV(1);    
    fbWRENCH(2) = 0*errV(2);    
    fbWRENCH(3) = -50*errV(3);     
    fbWRENCH(4) = 20*errV(4);
    fbWRENCH(5) = 20*errV(5);
    fbWRENCH(6) = 1*errV(6);
    
    
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

    dampWRENCH = -damp*vel.*(dampGain);    

    %%
    %Enable or disable components
    %%

    s = (expGRF(3)/expData.grfMAXZ)*smax;
    

    for i = 1:1:6
        %if(i ~= 5)
        posTRACKING(i) = posTRACKING(i).*(1-s);
        velTRACKING(i) = velTRACKING(i).*(1-s);
        %end
    end
    
    posToeTRACKING = posToeTRACKING.*(1-s);
    velToeTRACKING = velToeTRACKING.*(1-s);
    
    expWRENCH =   expWRENCH;
    fbWRENCH  =    fbWRENCH;
    dampWRENCH = dampWRENCH.*s;
    
    
    
    
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


tmp = xDotMex(t, x, vParams(1:50), vInputs);
xdot = tmp(1:1:length(x));