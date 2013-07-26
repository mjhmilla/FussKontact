function success = postprocess(tsol, xsol, contactInfo, vInput, ...
                   aniPoints, expData, vParams, vToe)

success = 0;

buildParamVariableList; %requires vParams and vToe

t0 = tsol(1);
t1 = tsol(length(tsol));
tani = [t0:(t1-t0)/(aniPoints-1):t1]';

xAni            = zeros(length(tani),size(xsol,2));
contactInfoAni  = zeros(length(tani),size(contactInfo,2));
vInputAni       = zeros(length(tani),size(vInput,2)); 

for i=1:1:size(xsol,2)
    xAni(:,i) = interp1(tsol,xsol(:,i),tani);
end

for i=1:1:size(contactInfo,2)
    contactInfoAni(:,i) = interp1(tsol,contactInfo(:,i),tani);
end

for i=1:1:size(vInput,2)
    vInputAni(:,i) = interp1(tsol,vInput(:,i),tani);
end


anklePosOri     = zeros(length(tani),13);
k1PosOri        = zeros(length(tani),13); 
k2PosOri        = zeros(length(tani),13);
metaJointPosOri = zeros(length(tani),13);
heelDiskPosOri  = zeros(length(tani),13);
foreDiskPosOri  = zeros(length(tani),13);

heelForceTorque = zeros(length(tani),10);
foreForceTorque = zeros(length(tani),10);
footForceTorque = zeros(length(tani),10);

for i = 1:1:length(tani)

%%
%turn the arrays into named variables
%%    
    %State Information
    x    = xAni(i,8);
    y    = xAni(i,9);
    z    = xAni(i,10);
    th   = xAni(i,11);
    zeta = xAni(i,12);
    eta   = xAni(i,13);
    xi  = xAni(i,14);
    
    %Input information    
    Fx = vInputAni(i,1);
    Fy = vInputAni(i,2);
    Fz = vInputAni(i,3);
    Mx = vInputAni(i,4);
    My = vInputAni(i,5);
    Mz = vInputAni(i,6);
    TK2cK2a = vInputAni(i,7);
    
    %Contact Information
    HFx = contactInfoAni(i,1);
    HFy = contactInfoAni(i,2);
    HFz = contactInfoAni(i,3);
    Hcopx  = contactInfoAni(i,4);
    Hcopy  = contactInfoAni(i,5);
    Hcopz  = contactInfoAni(i,6);
    
    FFx = contactInfoAni(i,7);
    FFy = contactInfoAni(i,8);
    FFz = contactInfoAni(i,9);
    Fcopx  = contactInfoAni(i,10);
    Fcopy  = contactInfoAni(i,11);
    Fcopz  = contactInfoAni(i,12);
    

%%
%call the script that calculates the position and orientation of every
%body of interest
%%
    calcPosVecsRotMatrices;
    %         rAnkle, RMAnkle
    %         rK1com, RMK1com
    %         rK2com, RMK2com
    %         rK1c, RMK1c
    %         rHeel, RMHeel
    %         rForeFoot, RMForeFoot
    
%%
%build the output matricies of time, position and rotation
%%

    anklePosOri(i,:)     = [tani(i) rAnkle matrix2rowVector(RMAnkle')];
    k1PosOri(i,:)        = [tani(i) rK1com matrix2rowVector(RMK1com')];        
    k2PosOri(i,:)        = [tani(i) rK2com matrix2rowVector(RMK2com')];        
    metaJointPosOri(i,:) = [tani(i) rK1c   matrix2rowVector(RMK1c')];
    heelDiskPosOri(i,:)  = [tani(i) rHeel  matrix2rowVector(RMHeel')];
    foreDiskPosOri(i,:)  = [tani(i) rForeFoot matrix2rowVector(RMForeFoot')];

%%
%build the force output matricies
%%
    if(HFx*HFx + HFy*HFy + HFz*HFz < 0.001)
        Hcopx = rHeel(1);
        Hcopy = rHeel(2);
        Hcopz = rHeel(3);
    end

    if(FFx*FFx + FFy*FFy + FFz*FFz < 0.001)
        Fcopx = rForeFoot(1);
        Fcopy = rForeFoot(2);
        Fcopz = rForeFoot(3);
    end
    

    heelForceTorque(i,:) = [tani(i) Hcopx Hcopy Hcopz HFx HFy HFz 0 0 0];
    foreForceTorque(i,:) = [tani(i) Fcopx Fcopy Fcopz FFx FFy FFz 0 0 0];

    FootF   = [(HFx+FFx) (HFy+FFy) (HFz+FFz)];
    FootCOP = [0 0 0];
    if(Hcopz < Fcopz)
       FootCOP = [Hcopx, Hcopy, Hcopz]; 
    else
       FootCOP = [Fcopx, Fcopy, Fcopz]; 
    end
    
    if(HFz + FFz > 0)
       FootCOP = ([Hcopx, Hcopy, Hcopz].*HFz + ...
                  [Fcopx, Fcopy, Fcopz].*FFz)./(HFz+FFz); 
       
    end
    footForceTorque(i,:) = [tani(i) FootCOP FootF 0 0 0];
    
    
end

%%
%write it all to file
%%
dlmwrite('animation/anklePosOri.dat',    anklePosOri,'\t');
dlmwrite('animation/k1PosOri.dat',       k1PosOri,'\t');
dlmwrite('animation/k2PosOri.dat',       k2PosOri,'\t');
dlmwrite('animation/metaJointPosOri.dat',metaJointPosOri,'\t');
dlmwrite('animation/heelDiskPosOri.dat', heelDiskPosOri,'\t');
dlmwrite('animation/foreDiskPosOri.dat', foreDiskPosOri,'\t');

dlmwrite('animation/heelForceTorque.dat', heelForceTorque,'\t');
dlmwrite('animation/foreForceTorque.dat', foreForceTorque,'\t');
dlmwrite('animation/footForceTorque.dat', footForceTorque,'\t');

dlmwrite('animation/heelDisk.fr',[rH rH rH],' ');
dlmwrite('animation/foreDisk.fr',[rF rF rF],' ');

v1 = ones(size(tani));
v0 = zeros(size(tani));

dlmwrite('animation/origin.dat',...
    [tani v0 v0 v0 v1 v0 v0 v0 v1 v0 v0 v0 v1],...
    '\t');

thX = pi/2;
thY = -pi/2;
thZ = -pi/2;

rX = [1 0 0; 0 cos(thX) -sin(thX); 0 sin(thX) cos(thX)];
rY = [cos(thY) 0 sin(thY); 0 1 0; -sin(thY) 0 cos(thY)];
rZ = [cos(thZ) -sin(thZ) 0; sin(thZ) cos(thZ) 0; 0 0 1];

RM = rY*rZ;

camera = [tani -v1.*(0.5) anklePosOri(1,3).*v1 v1.*0.05...
               RM(1,1).*v1 RM(1,2).*v1 RM(1,3).*v1...
               RM(2,1).*v1 RM(2,2).*v1 RM(2,3).*v1...
               RM(3,1).*v1 RM(3,2).*v1 RM(3,3).*v1];
dlmwrite('animation/camera.dat',camera,'\t');

%%
% If the experimental data structure is not empty, then
%   1. Animate experimental ankle frame
%   2. Positions of the markers
%   3. Force vector
%
%  But, all in grey
%
%%
if(isempty(expData) ~= 1)
    %resample data
    expAk = zeros(length(tani),3);
    expAkR= zeros(length(tani),9);
    
    expGRF = zeros(length(tani),3);
    expCOP = zeros(length(tani),3);    
    expMkr = zeros(14,3,length(tani));

    
    v1 = ones(length(tani),1);
    v0 = zeros(length(tani),1);
    
    RI = [v1 v0 v0 v0 v1 v0 v0 v0 v1];
    
    for i=1:1:3
        expAk(:,i) = interp1(expData.time, expData.heel.xyz(:,i), tani);
        expGRF(:,i)= interp1(expData.time, expData.grf(:,i), tani);
        expCOP(:,i)= interp1(expData.time, expData.cop(:,i), tani);
        for j=1:1:14
           tmp = reshape(expData.mkrs(j,i,:),length(expData.time),1);
           expMkr(j,i,:) = interp1(expData.time,tmp,tani);
        end
    end
    
    for i=1:1:9
       expAkR(:,i) = interp1(expData.time, expData.heel.R(:,i), tani); 
    end
    
    dlmwrite('animation/expAnkle.dat', [tani expAk expAkR],'\t');
    dlmwrite('animation/expGRF.dat',[tani expCOP expGRF v0 v0 v0],'\t');
    
    for i=1:1:14
       xtmp = reshape(expMkr(i,1,:), length(tani),1);
       ytmp = reshape(expMkr(i,2,:), length(tani),1);
       ztmp = reshape(expMkr(i,3,:), length(tani),1);
       
       dlmwrite(['animation/mkr',num2str(i),'.dat'],...
                [tani xtmp ytmp ztmp RI], '\t');
    end
    
end
%%
%compile it
%%
cd('animation');
if(isempty(expData) ~= 1)
    dos('compileFootMdlExp.bat');
    dos('compileFootMdl.bat');
else
    dos('compileFootMdl.bat');
end

cd ..

success = 1;