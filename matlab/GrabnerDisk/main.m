clc; 
close all;
clear all;

%%
%Script Configuration
%%
aniFreq = 1000;
tsim = [0 1];

rootDir = pwd;

%%
% Set up initial conditions
%%
vx = 0;
vy = 0;
vz = 0;
wx = 0;
wy = 0;
wz = 0;
x = 0;
y = 0;
z = 0.707106271093168148;
alpha = 0;
beta = pi/4;
zeta = 0;
Wvx = 0;
Wvy = 0;
Wvz = 0;
Wwx = 0;
Wwy = 0;
Wwz = 0;

calcXdot = @(t,x) xdot_disk(t,x, 0);
calcContactForce = @(t,x) xdot_disk(t,x, 1);


ic = [vx,vy,vz, wx,wy,wz,...
      x,y,z, alpha, beta, zeta,...
      Wvx,Wvy,Wvz, Wwx,Wwy,Wwz];

disp(sprintf('x0: vx(%f),vy(%f),vz(%f)',vx,vy,vz));
disp(sprintf('    wx(%f),wy(%f),wz(%f)',wx,wy,wz));
disp(sprintf('     x(%f), y(%f), z(%f)',x,y,z));
disp(sprintf('     a(%f), b(%f),ze(%f)\n\n',alpha,beta,zeta));

xdot = calcXdot(0,ic);  

disp(sprintf('dx0:dvx(%f), dvy(%f), dvz(%f)',xdot(1),xdot(2),xdot(3)));
disp(sprintf('    dwx(%f), dwy(%f), dwz(%f)',xdot(4),xdot(5),xdot(6)));
disp(sprintf('     dx(%f),  dy(%f),  dz(%f)',xdot(7),xdot(8),xdot(9)));
disp(sprintf('     da(%f),  db(%f), dze(%f)\n\n',xdot(10),xdot(11),xdot(12)));

pos = calcPosition(ic);
rm  = calcRotationMatrix(ic);
disp(sprintf('r: x(%f), y(%f), z(%f)', pos(1),pos(2),pos(3)));
disp(sprintf('RM  %f , %f , %f ', rm(1),rm(2), rm(3)));
disp(sprintf('    %f , %f , %f ', rm(4),rm(5), rm(6)));
disp(sprintf('    %f , %f , %f ', rm(7),rm(8), rm(9)));
disp('Maple generates the rotation matrix R02');


%%
% Simulate
%%
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t sol] = ode45(calcXdot,tsim,ic,options);

%%
% Post process to compute KE+PE-W
%%
m = 1;
J = 0.1; %Jxx = Jyy = Jzz = 0.1. All off diagional terms are zero
G = 9.81;

ke = (0.5*m).*( sol(:,1).*sol(:,1) ...
              + sol(:,2).*sol(:,2) ...
              + sol(:,3).*sol(:,3)) ...
   + (0.5.*J).*( sol(:,4).*sol(:,4) ...
               + sol(:,5).*sol(:,5) ...
               + sol(:,6).*sol(:,6));

pe = m*G*sol(:,9);           
           
w = sol(:,13) + sol(:,14) + sol(:,15)...
   +sol(:,16) + sol(:,17) + sol(:,18);

kepew = ke+pe-w;

fig1 = figure;
subplot(2,2,1);
    plot(t,kepew,'b');
    title('Ke+Pe-W');
    xlabel('Time (s)');
    ylabel('Energy (J)')
subplot(2,2,2);
    plot(t,ke,'b');
    hold on;
    plot(t,pe,'r');    
    title('Ke(b), Pe(r)');
    xlabel('Time (s)');
    ylabel('Energy (J)')
subplot(2,2,3);
    plot(t,sol(:,13),'r');
    hold on;
    plot(t,sol(:,14),'g');
    hold on;
    plot(t,sol(:,15),'b');        
    title('Linear Work');
    xlabel('Time (s)');
    ylabel('Energy (J)');
subplot(2,2,4);
    plot(t,sol(:,16),'r');
    hold on;
    plot(t,sol(:,17),'g');
    hold on;
    plot(t,sol(:,18),'b');        
    title('Angular Work');
    xlabel('Time (s)');
    ylabel('Energy (J)');
    
%%
% Post process to generate files to animate results
%%
%resample
tspan = tsim(2)-tsim(1);
npts = length(t);
if(aniFreq > 0)
    npts = floor(tspan*aniFreq);
end
tani = [tsim(1):tspan/(npts-1):tsim(2)]';
solAni = zeros(npts,12);

for i=1:1:12
    solAni(:,i) = interp1(t,sol(:,i),tani);
end



%Get contact forces
contactForce = zeros(npts,9);
for i=1:1:npts
    contactForce(i,1:6) = calcContactForce(tani(i),solAni(i,:))';
end

%Export animation files
rcom = calcPosition(solAni);
RMcom = calcRotationMatrix(solAni);
RMcomt = [RMcom(:,1) RMcom(:,4) RMcom(:,7) ...
          RMcom(:,2) RMcom(:,5) RMcom(:,8) ...
          RMcom(:,3) RMcom(:,6) RMcom(:,9)];
dlmwrite('animation/disk.dat',[tani,rcom,RMcomt],'\t');

dlmwrite('animation/forces.dat',[tani,contactForce],'\t');

v1 = ones(size(tani));
v0 = zeros(size(tani));
camera = [tani v0 v1.*(-5) v1.*(1)  v1 v0 v0 v0 v0 -v1 v0 v1 v0];
dlmwrite('animation/camera.dat',camera,'\t');

%cd('animation');
%dos('compileWRL.bat');
%cd ..

