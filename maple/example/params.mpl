leParams:= []

G, Jzz_ffcom, R_mt, c_mt, cx, cy, cz, fFMTcx, fFMTcy, fFMTcz, fFMTx, fFMTy, fFMTz, hFAx, hFAy, hFAz, hFMTx, hFMTy, hFMTz, k_mt, m_ffcom, mt_d, mt_k, mt_theta0, mu, mud, nx, ny, nz, sx, sy, sz, tx1, tx2, ty1, ty2, tz1, tz2, velD, velS

G         = 9.81, 
Jzz_ffcom  = 0.0055, 
Jzz_hfcom  = 0.0055,

R_hl      = 0.04, 
R_mt      = 0.015, 

c_hl      = 1.73, 
c_mt      = 1.73, 

cbf_hl    = 100, 
cbf_mt    = 100, 

cx        = 0, 
cy        = 0, 
cz        = 0, 

hFAx      = -0.0575, 
hFAy      = 0.015, 
hFAz      = 0, 

hFHcx     = -0.0775, 
hFHcy     = -0.053,  
hFHcz     = 0, 

fFMTcx    = 0.0575, 
fFMTcy    = -0.015, 
fFMTcz    = 0, 

hFMTx = 0,
hFMTy = 0,
hFMTz = 0,

fFMTx = -0.08,
fFMTy = 0.01,
fFMTz = 0,

mt_k = 1,
mt_d = 0.01,
mt_theta0 = 0,


k_hl      = 25029995.73, 
k_mt      = 268200000,  

kbf_hl    = 1000, 
kbf_mt    = 1000, 

m_ffcom    = 0.5, 
m_hfcom    = 0.5,

mubf      = 0.5, 
musbf     = 0.7, 

nx        = 0, 
ny        = 1, 
nz        = 0, 

pdthetabf_hl      = 0, 
pdthetabf_mt      = 0, 

 
resetbf_hl        = 2, 
resetbf_mt        = 2, 


svbf              = 0.05, 

sx                = 0, 
sy                = 0, 
sz                = 0, 

thetaDistbf_hl    = 0, 
thetaDistbf_mt    = 0, 

thetabf_hl        = 0, 
thetabf_mt        = 0, 

tx1       = 1, 
tx2       = 0, 
ty1       = 0, 
ty2       = 0, 
tz1       = 0, 
tz2       = 1, 

xbf_hl    = 0, 
xbf_mt    = 0, 

ybf_hl    = 0, 
ybf_mt    = 0, 
 
zbf_hl    = 0, 
zbf_mt    = 0,  

MT_T(t) = 0,

X_FM(t)=0,
Y_FM(t)=0,
theta_akFM(t)=0,

diff(X_FM(t),t)=0,
diff(Y_FM(t),t)=0,
diff(theta_akFM(t),t)=0,

diff(diff(X_FM(t),t),t)=0,
diff(diff(Y_FM(t),t),t)=0,
diff(diff(theta_akFM(t),t),t)=0




