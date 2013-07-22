function status = getIntegratorOutput(t,y,flag)
status = 0;
if(isempty(t)==0)
    if(t(1) > 0)
        count = length(sprintf('Time: %f \r',t(1)));
        fprintf(1, repmat('\b',1,count));
    end
    
    fprintf('Time: %f \r',t(1));
end