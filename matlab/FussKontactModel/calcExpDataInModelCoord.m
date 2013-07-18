function pData = calcExpDataInModelCoord(expFile)
%%
% Will read in processed *.mat files from a previous experiment (record
% 3D foot kinematics and ground reaction forces) and will transform them
% into the model coordinates for this model. 
%
% @param expFile: string to a *.mat file generated from the 5 June 2013
%                 foot contact data recording in Essen
%
%
% @return pData 
%           expData.state: mapping of data state to the model's coordinates
%
%
%
%%
    toeOffsetAngle = -1.972594094578409e-001;

    pData = [];
    if(isempty(expFile)==1)
       return; 
    end
    
    
%%    
%1. Read in file setup output structure. 
%%
   expData = load(expFile);
   
   rows = size(expData.data.heelV,1);   %
   pData.time  = zeros(rows,1);         %
   pData.heel.xyz  = zeros(rows,3);     %1b
   pData.heel.R    = zeros(rows,9);     %1b
   pData.heel.vxyz = zeros(rows,3);     %1c
   pData.heel.vR   = zeros(rows,9);     %1c 
   pData.meta.R    = zeros(rows,9);     %1b
   pData.meta.vR   = zeros(rows,9);     %1c 
   pData.grf   = zeros(rows,3);         %1a
   pData.cop   = zeros(rows,3);         %1a
   pData.mkrs  = zeros(14,3,rows);      %1b
   pData.mkrNames = expData.data.mkrName; %here
   
   pData.mdlState = zeros(rows,14);     %2
   
%%   
%1a. Resample ground reaction forces to have the same sense of time
%%
   cameraRate   = expData.data.metadata.cameraRate;
   grfTime      = expData.data.grf.time;  
   elapsedTime  = max(grfTime)-min(grfTime);
   pData.time   = double([0:1/cameraRate:elapsedTime]');
   
   for i=1:1:3
    pData.grf(:,i)  = interp1(grfTime, expData.data.grf.f(:,i),...
                                    pData.time); 
    pData.cop(:,i)  = interp1(grfTime, expData.data.grf.cop(:,i)./1000.0,...
                                pData.time);                        
   end
%%
%1b. Low pass filter and copy 
%%
    pData.heel.xyz = double(expData.data.heelV);
    pData.heel.R   = double(expData.data.heelR);
    pData.meta.R   = double(expData.data.foreR);
    
    for i=1:1:rows
       for j=1:1:14
          for k=1:1:3
              %Time is the last vector so that columns can be 
              %picked off without running yet another for loop
             pData.mkrs(j,k,i) = expData.data.allVideoData(j,i,k)./1000.0; 
          end
       end
    end
    
%%    
%1c. Take derivatives of the ankle frame and the rotation matrix vector
%    for later use.
%%
    for i=1:1:3
       pData.heel.vxyz(:,i) = dspline(pData.time, pData.heel.xyz(:,i),...
                                      pData.time);
    end
    for i=1:1:9
       pData.heel.vR(:,i) = dspline(pData.time, pData.heel.R(:,i),...
                                    pData.time); 
       pData.meta.vR(:,i) = dspline(pData.time, pData.meta.R(:,i),...
                                    pData.time); 
    end
%% 
%2. Loop through and map the experimental foot state to the 
%   model's state. For now ignore the metatarsal joint
%%

disp('  Metatarsal angle of model approximately matches experimental data');

for i=1:1:length(pData.time)
    vx0     = 0;
    vy0     = 0;
    vz0     = 0;
    dth0    = 0;
    wx0     = 0;
    wy0     = 0;
    wz0     = 0;
    x0      = 0;
    y0      = 0;
    z0      = 0;
    th0     = 0;
    zeta0   = 0;
    eta0    = 0;
    xi0     = 0;
    
    %%
    %Ankle position
    %%
    rAk = expData.data.heelV(i,:);
    x0 = rAk(1);
    y0 = rAk(2);
    z0 = rAk(3);
    
    
    %%
    %Foot orientation in EA123
    % n.b. Rotation matrix from the local frame to the global frame
    %      The transpose is necessary as heelR contains the rotation
    %      matrix from the global to the local frame.
    RMAk= [expData.data.heelR(i,1:3);...
          expData.data.heelR(i,4:6);...
          expData.data.heelR(i,7:9)]';

    %n.b.Again, Dynaflexpro uses an odd sign convention for all of its basic
    %   rotation matrices (convenient for constructing R01 without taking
    %   any transposes).
    zeta0 = atan2( -RMAk(3,2), RMAk(3,3) );   
    eta0  = atan2(RMAk(3,1)*cos(zeta0), RMAk(3,3));   
    Rx = [1 0 0; 0 cos(zeta0) sin(zeta0); 0 -sin(zeta0) cos(zeta0)];
    Ry = [cos(eta0) 0 -sin(eta0); 0 1 0; sin(eta0) 0 cos(eta0)];
   
    %R = Rz*Ry*Rx
    %Rz = R*Rx'*Ry'
    Rz = RMAk * Rx' * Ry';   
    xi0 = atan2(-Rz(2,1),Rz(1,1));

    %%
    %Metatarsal angle
    %  The frame on the metatarsal has been constructed so that the 
    %  x axis on the toe is parallel to the x axis on the heel
    %%
    RMmeta = [expData.data.foreR(i,1:3);...
              expData.data.foreR(i,4:6);...
              expData.data.foreR(i,7:9)]'; 
    RMmetafoot = RMmeta*RMAk';
          
    %Negative sign because the basic rotation matrix convention used by
    %DFP is from the local frame to the global (negative sign in Y col)
    th0 = atan2(-RMmetafoot(3,2),RMmetafoot(2,2)) - toeOffsetAngle;
    
    %%
    %Ankle velocity
    %%
    vx0 = pData.heel.vxyz(i,1);
    vy0 = pData.heel.vxyz(i,2);
    vz0 = pData.heel.vxyz(i,3);
    
    %%
    %Ankle angular velocity
    % rp = r + As'p                      [1]
    % d/dt rp = d/dt r + d/dt A s'p      [2]   
    % n.b.
    %      I = A^T A                     [3]
    %      0 = d/dt A^T A + A^T d/dt A   [4]
    %      
    %      s'p = A^T sp                  [5]
    % 
    % Now 2 reads      
    % d/dt rp = d/dt r + d/dt A A^T sp   [6]
    %  
    % Since in vector form
    %
    % d/dt rp = d/dt r + w x sp          [7]
    %
    % ~w = d/dt A A^T                    [8]
    %
    % ~w = [0  -wz wy]
    %      [wz  0 -wx]
    %      [-wy wx  0]    
    %%
    
    dRMAk = [pData.heel.vR(i,1:3);...
             pData.heel.vR(i,4:6);
             pData.heel.vR(i,7:9)];
    %RMAk is equivalent to A: rotates from the local frame to the global    
    omegaSS = dRMAk*(RMAk); 
    wx0 = omegaSS(3,2);
    wy0 = omegaSS(1,3);
    wz0 = omegaSS(2,1);
    
    %Metatarsal angular velocity: compute the angular velocity of
    %the meta tarsal frame, project it onto the angle X axis
    %take the difference between this and wx of the foot.
    dRMMeta = [pData.meta.vR(i,1:3);...
               pData.meta.vR(i,4:6);
               pData.meta.vR(i,7:9)];
    omegaSSM = dRMMeta*(RMmeta);
    wxM0 = omegaSSM(3,2);
    wyM0 = omegaSSM(1,3);
    wzM0 = omegaSSM(2,1);
    
    dth0 = (wxM0-wx0);
    
    pData.mdlState(i,:) = double([vx0 vy0 vz0 dth0   wx0  wy0 wz0 ...
                            x0  y0  z0  th0 zeta0 eta0 xi0]);
   
   
   
end


