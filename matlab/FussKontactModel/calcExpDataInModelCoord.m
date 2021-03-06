function pData = calcExpDataInModelCoord(expFile, vToe, grfFiltFreq)
%%
% Will read in processed *.mat files from a previous experiment (record
% 3D foot kinematics and ground reaction forces) and will transform them
% into the model coordinates for this model. This function also trims the
% data so that only times where the model is in contact with the ground are
% included.
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
    toeOffsetAngle = vToe(3);

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
   pData.wrench= zeros(rows,6); %Fx,Fy,Fz, Mx, My, Mz
   pData.wrenchZ= zeros(rows,6); %Wrench caused by only the normal force
   pData.grf   = zeros(rows,3);         %1a
   pData.cop   = zeros(rows,3);         %1a
   pData.mkrs  = zeros(14,3,rows);      %1b
   pData.mkrNames = expData.data.mkrName; %here
   
   pData.mdlState = zeros(rows,14);     %2
   pData.mdlStateOffset = zeros(1,14);
%%   
%1a. Resample ground reaction forces to have the same sense of time
%%
   cameraRate   = expData.data.metadata.cameraRate;
   grfTime      = expData.data.grf.time;     
   elapsedTime  = max(grfTime)-min(grfTime);
   pData.time   = double([0:1/cameraRate:elapsedTime]');
   grfRate = cameraRate*(length(grfTime)/length(pData.time));
   
   if(isempty(grfFiltFreq) == 0)
       disp('  Filtering the experimental GRF and COP');
       [b a] = butter(2, double(grfFiltFreq*2.0/grfRate));
   end
   
   for i=1:1:3
    grfTmp = [];
    copTmp = [];
    
    if(isempty(grfFiltFreq) == 0)
      grfTmp = filtfilt(b,a,   double(expData.data.grf.f(:,i)));
      copTmp = filtfilt(b,a,   double(expData.data.grf.cop(:,i)./1000.0));
    else
       grfTmp =  double(expData.data.grf.f(:,i));
       copTmp =  expData.data.grf.cop(:,i)./1000.0;
    end
    
    pData.grf(:,i)  = interp1(grfTime, grfTmp,...
                                    pData.time); 
    pData.cop(:,i)  = interp1(grfTime, copTmp,...
                                pData.time);
   end

    idx0 = 1;
    while(expData.data.grf.f(idx0,3) <= 0)
        idx0 = idx0+1;
    end
    
    idx1 = length(expData.data.grf.f(:,3));
    while(expData.data.grf.f(idx1,3) <= 0)
        idx1 = idx1-1;
    end
    
    idx = 1;
    while(pData.time(idx) < grfTime(idx0))
       pData.grf(idx,:) = [0,0,0]; 
       pData.cop(idx,:) = [0,0,0];
       idx = idx+1;
    end

    idx = length(pData.time);
    while(pData.time(idx) > grfTime(idx1))
       pData.grf(idx,:) = [0,0,0]; 
       pData.cop(idx,:) = [0,0,0];
       idx = idx-1;
    end
   
   grfM = (pData.grf(:,1).*pData.grf(:,1)...
         +pData.grf(:,2).*pData.grf(:,2)...
         +pData.grf(:,3).*pData.grf(:,3)).^0.5;
   
    %Get the time of contact onset and break. 
    %All data will be trimmed to these values
    t0 = 1;
    t1 = length(pData.time);
    while(grfM(t0) < 1)
        t0 = t0+1;
    end
    
    while(grfM(t1) < 1)
       t1 = t1-1; 
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

%%
%Compute the experimental ankle wrench.
%%
xIdx = 8;
yIdx = 9;
zIdx = 10;


for i=1:1:length(pData.time)
    pData.wrench(i,1) = -pData.grf(i,1);
    pData.wrench(i,2) = -pData.grf(i,2);
    pData.wrench(i,3) = -pData.grf(i,3);
    
    pData.wrenchZ(i,1) = 0;
    pData.wrenchZ(i,2) = 0;
    pData.wrenchZ(i,3) = -pData.grf(i,3);    
    
    grfM = sum(pData.grf(i,:).*pData.grf(i,:)).^0.5;
    if(grfM > 0.0)
        %Compute the moment the grf generates about the ankle
        r = [0;0;0];
        r(1) = pData.cop(i,1)-pData.mdlState(i,xIdx);
        r(2) = pData.cop(i,2)-pData.mdlState(i,yIdx);
        r(3) = pData.cop(i,3)-pData.mdlState(i,zIdx);
        tq = cross(pData.grf(i,:)',r);

        pData.wrench(i,4) = tq(1);
        pData.wrench(i,5) = tq(2);
        pData.wrench(i,6) = tq(3);
        
        grfZ = [0;0;0];
        grfZ(3) = pData.grf(i,3);
        tq = cross(grfZ, r);
        
        pData.wrenchZ(i,4) = tq(1);
        pData.wrenchZ(i,5) = tq(2);
        pData.wrenchZ(i,6) = tq(3);
        
    end
end

%%
%Trim all data so that the foot is always in contact
%%

%t0 = min(find(pData.time > 10));
%t1 = length(pData.time);

pData.time      = pData.time(t0:1:t1)-pData.time(t0);
pData.wrench    = pData.wrench(t0:1:t1,:);
pData.wrenchZ   = pData.wrenchZ(t0:1:t1,:);
pData.heel.xyz  = pData.heel.xyz(t0:1:t1,:); 
pData.heel.R    = pData.heel.R(t0:1:t1,:); 
pData.heel.vxyz = pData.heel.vxyz(t0:1:t1,:); 
pData.heel.vR   = pData.heel.vR(t0:1:t1,:); 
pData.meta.R    = pData.meta.R(t0:1:t1,:); 
pData.meta.vR   = pData.meta.vR(t0:1:t1,:); 
pData.grf       = pData.grf(t0:1:t1,:); 
pData.cop       = pData.cop(t0:1:t1,:);
pData.mkrs      = pData.mkrs(:,:,t0:1:t1);
pData.mdlState  = pData.mdlState(t0:1:t1,:);
pData.grfMAXZ = max(pData.grf(:,3));

