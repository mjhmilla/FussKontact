function vIC = getIC(expFile)

vx0   = 0.0;
vy0   = 0.0;
vz0   = 0.0;
dth0  = 0.0;
wx0   = 0.0;
wy0   = 0.0;
wz0   = 0.0;
x0    = 0.0;
y0    = 0.0;
z0    = 0.15;
th0   = 0.0;
zeta0 = -0.2;
eta0   = 0.0;
xi0  = 0.0;




vIC = double([vx0, vy0, vz0, dth0, wx0, wy0, wz0, ...
       x0, y0, z0, th0, zeta0, eta0, xi0]');
   