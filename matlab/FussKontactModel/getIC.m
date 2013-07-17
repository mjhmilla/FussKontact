function vIC = getIC(expFile)

vx0   = 0;
vy0   = 0;
vz0   = 0;
dth0  = 0;
wx0   = 0;
wy0   = 0;
wz0   = 0;
x0    = 0;
y0    = 0;
z0    = 0.1;
th0   = 0;
zeta0 = 0;
eta0   = 0;
xi0  = 0;

if(isempty(expFile) ~= 1)
   %1. Read in file
   expData = load(expFile);
   %2. Get the ankle position and orientation at time 0
   rAk = expData.data.heelV(1,:);
   RMAk= [expData.data.heelR(1,1:3);...
          expData.data.heelR(1,4:6);...
          expData.data.heelR(1,7:9)];
   %3. Get the rotation matrix of the foot,
   %   decompose it into zeta eta xi
   
   
   %4. Get the ankle velocity at time t0
   %5. Get the ankle angular velocity at time t0
end



vIC = [vx0, vy0, vz0, dth0, wx0, wy0, wz0, ...
       x0, y0, z0, th0, zeta0, eta0, xi0]';
   