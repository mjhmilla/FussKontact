clc;
clear all;
close all;

flag_mode = 2;
% 0: preprocess experimental data
% 1: parameter optimize & animate
% 2: just animate
optMaxIter = 50;
paramsFile = [];%'modeldata\optvParamsvToe20130724.mat'; %has positive toe angle
resultsFolder = 'modeldata\lima20130726';

preExpNames = {};

%%
% Fetch Model Parameters (required for all tasks)
%%
vParam = [];
vToe = [];

if(isempty(paramsFile) == 1)
    [vParam vToe] = getParams();
    disp('Using initial parameters');
else
    load(paramsFile);
    vParam = optOutput.vParamOpt;
    vToe   = optOutput.vToe;
    disp(['Using parameters from: ', paramsFile]);
end
paramIdx;

%%
%Preprocessing
%%
if(flag_mode == 0)    
    for i=1:1:length(preExpNames)
        expData=calcExpDataInModelCoord(['data/',preExpNames{i}],vToe,[]);
        save(['modeldata/',preExpNames{i}], 'expData'); 
    end
end

%%
% Configure the parameter optimization run
%%


optParamIdx = [2 kmt_idx];
x0 = zeros(size(optParamIdx,1),1);
xscaling = [100];

distScaling  = 100^2; %Will make an error of 1cm yield an error of 1
angleScaling = (1/(10*pi/180))^2; %Will make an error of 10 degrees
                                  %yield an error of 1

errorScaling = [distScaling;  distScaling;  distScaling; ...
                angleScaling; angleScaling; angleScaling];

errorTimeScaling = 1e1; %Will make an error of a simulation cut short by
                       %36% yield an error of 1
                       
expFileNames = {'modeldata\Rotations12.mat'};
%'modeldata\Walking01.mat', ...
% 'modeldata\Jogging02.mat', ... 
% 10,2,1; 100,10,1;
% 0 1; 0 1; 
%
% 'modeldata\Rotations10.mat',...
% 10 10 0.5;
% 0 1;
%
%'modeldata\Rotations02.mat',...
% 10 10 0.5; 
% 0 1;

expCtrlGains = [10 10 0.5];            
            
expTimeRange = [0 0.6];

errFcn = @(x)calcFootError(x, xscaling, optParamIdx, vParam,vToe, ...
                             errorScaling, errorTimeScaling,...
                             expFileNames, expCtrlGains, ...
                             expTimeRange, resultsFolder, 0);
                         
                         
outFcn = @(x,optimValues,state)optimizeFootOutputFunction(x,...
                 optimValues,state, xscaling, optParamIdx, vParam, vToe,...
                 resultsFolder);                         
                         
%With all three motions 50 iterations is roughly 12 hours of optimization             
options = optimset('Display','iter','MaxIter',optMaxIter,...
                   'PlotFcns',@optimplotfval, ...
                   'OutputFcn',outFcn);
                                    
%%
% Optimize!
%%
x = [];
if(flag_mode == 1)               
    [x, fval, exitflag] = fminsearch(errFcn, x0, options);    
else
   x = x0; 
end

%%
%Animate the results
%%
if(flag_mode >= 1)
err = calcFootError(x, xscaling, optParamIdx, vParam,vToe, ...
                             errorScaling, errorTimeScaling,...
                             expFileNames, expCtrlGains, ...
                             expTimeRange, resultsFolder, 1);                             
end                    