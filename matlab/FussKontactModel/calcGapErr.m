function errGAP = calcGapErr(x, vIC, vParams)

vIC(12) = vIC(12) + x(1)/1000; %Zeta


tmp = xDotMex(0, vIC, vParams, [0 0 0 0 0 0 0]');
xdot = tmp(1:1:length(vIC));
contactInfo = tmp(length(vIC)+1:1:length(tmp));
grfCOP = calcModelGRFCOP(contactInfo);

gapFore = contactInfo(12);
gapHeel = contactInfo(6);

errGAP = 10e5*(gapFore-gapHeel)^2;