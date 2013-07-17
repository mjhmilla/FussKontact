function success = postprocess(tsol, xsol, contactInfo, vInput, aniFreq)

success = 0;

getParams;

t0 = tsol(1);
t1 = tsol(length(tsol));
tani = [t0:1/aniFreq:t1]';

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
    if(HFx*HFx + HFy*HFy + HFz*HFz < 0.01)
        Hcopx = rHeel(1);
        Hcopy = rHeel(2);
        Hcopz = rHeel(3);
    end

    if(FFx*FFx + FFy*FFy + FFz*FFz < 0.01)
        Fcopx = rForeFoot(1);
        Fcopy = rForeFoot(2);
        Fcopz = rForeFoot(3);
    end
    

    heelForceTorque(i,:) = [tani(i) Hcopx Hcopy Hcopz HFx HFy HFz 0 0 0];
    foreForceTorque(i,:) = [tani(i) Fcopx Fcopy Fcopz FFx FFy FFz 0 0 0];

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

dlmwrite('animation/heelDisk.fr',[rH rH rH],' ');
dlmwrite('animation/foreDisk.fr',[rF rF rF],' ');

v1 = ones(size(tani));
v0 = zeros(size(tani));

thX = pi/2;
thY = pi/2;
thZ = pi/2;

rX = [1 0 0; 0 cos(thX) -sin(thX); 0 sin(thX) cos(thX)];
rY = [cos(thY) 0 sin(thY); 0 1 0; -sin(thY) 0 cos(thY)];
rZ = [cos(thZ) -sin(thZ) 0; sin(thZ) cos(thZ) 0; 0 0 1];

RM = rY*rZ;

camera = [tani v1.*(0.5) v0 v1.*0.05...
               RM(1,1).*v1 RM(1,2).*v1 RM(1,3).*v1...
               RM(2,1).*v1 RM(2,2).*v1 RM(2,3).*v1...
               RM(3,1).*v1 RM(3,2).*v1 RM(3,3).*v1];
dlmwrite('animation/camera.dat',camera,'\t');

%%
%compile it
%%
%cd('animation');
%dos('compileWRL.bat');
%cd ..

success = 1;