clc;
close all;
clear all;

getParams; %vParams
getInputs; %vInputs
getIC;     %vIC
t0 = 0;

xdot = xDotMex(t0, vIC, vParams, vInputs);

%disp(xdot);
