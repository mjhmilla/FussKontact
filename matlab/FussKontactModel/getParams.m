G = 9.81;
Jxx1 = .18e-3;
Jyy1 = .18e-3;
Jzz1 = .18e-3;
Jxx2 = .18e-3;
Jyy2 = .18e-3;
Jzz2 = .18e-3;
K1ax = 0;
K1ay = -.36e-1;
K1az = .41e-1;
K1brx = 0;
K1bry = 0;
K1bx = 0;
K1by = -.36e-1;
K1bz = -.41e-1;
K1cx = 0;
K1cy = .72e-1;
K1cz = -.41e-1;
K1crz =0;
K2ax = 0;
K2ay = 0;
K2az = 0;
aF = 20;
aH = 20;
cF = 1.0;
cH = 1.0;
dyVel = .1;
kF = 1000;
kH = 1000;
m1 = .75;
m2 = .25;
mud = .7;
mus = 1;
nx = 0;
ny = 0;
nz = 1;
pF = 1;
pH = 1;
rF = .5e-1;
rH = .4e-1;
stVel = .5e-1;
sx = 0;
sy = 0;
sz = 0;
t1x = 1;
t1y = 0;
t1z = 0;
t2x = 0;
t2y = 1;
t2z = 0;


vParams = [G,Jxx1,Jxx2,...
             Jyy1,Jyy2,...
             Jzz1,Jzz2,...
             K1ax,K1ay,K1az,...
             K1brx,K1bry,K1bx,K1by,K1bz,...
             K1crz,K1cx,K1cy,K1cz,...
             K2ax,K2ay,K2az,...
             aF,aH,...
             cF,cH,...
             dyVel,...
             kF,kH,...
             m1,m2,...
             mud,mus,...
             nx,ny,nz,...
             pF,pH,...
             rF,rH,...
             stVel,...
             sx,sy,sz,...
             t1x,t1y,t1z,...
             t2x,t2y,t2z]';
