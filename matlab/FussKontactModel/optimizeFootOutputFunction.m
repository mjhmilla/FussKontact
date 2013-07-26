function stop = optimizeFootOutputFunction(x, optimValues,  state, ...
                                    xscaling, optParamIdx, vParam,vToe,...
                                    resultsFolder)
                                
stop = 0;

for i=1:1:length(x)
   if(optParamIdx(i,1) == 1)
      pIdx = optParamIdx(i,2);
      vParam(pIdx) = vParam(pIdx) + x(i)/xscaling(i);
   end
   
   if(optParamIdx(i,1) == 2)
       pIdx = optParamIdx(i,2);
       vToe(pIdx) = vToe(pIdx) + x(i)/xscaling(i);
   end    
end


optOutput.fval        = optimValues.fval; 
optOutput.optParamIdx = optParamIdx;
optOutput.x           = x;
optOutput.xscaling    = xscaling;
optOutput.vParamOpt   = vParam;
optOutput.vToe        = vToe;

iter = optimValues.iteration;

save([resultsFolder,'\optvParamsvToe_',num2str(iter),'.mat'], 'optOutput');