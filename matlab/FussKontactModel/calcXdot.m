function xdot = calcXdot(t,x,vParams,vInputs)

tmp = xDotMex(t, x, vParams, vInputs);
xdot = tmp(1:1:length(x));