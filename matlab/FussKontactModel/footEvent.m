function [value, isterminal, direction] = footEvent(t,x,vParams,expData)

contactInfo = calcContactForcePosition(t,x,vParams,[0,0,0,0,0,0,0]');

heelF = contactInfo(1:3);
foreF = contactInfo(7:9);
heelCOP=contactInfo(4:6);
foreCOP=contactInfo(10:12);

grf = heelF+foreF;



expGRF = zeros(3,1);
expCOP = zeros(3,1);
for i=1:1:3
   expGRF(i) = interp1(expData.time, expData.grf(:,i),t);
   expCOP(i) = interp1(expData.time, expData.cop(:,i),t);
end

grfErr = expGRF-grf;
value = 0;
if(sum(grfErr.*grfErr) > 100)
   value = 1; 
   disp('Simulation Terminated: Error is too large');
end

isterminal = 1;
direction = 0;

