function [value, isterminal, direction] = footEvent(t,x,expData)

expX = zeros(size(x));
for i=1:1:14
    expX(i) = interp1(expData.time, expData.mdlState(:,i), t);
end
errX = x - expX;

rErr   = sum( errX(8:10).^2 )^0.5;
angErr = max(abs( errX(12:14) ));

value = 0;
if(rErr > 0.15 || angErr > pi/3)
 value = 1;
end

isterminal = 1;
direction = 0;

