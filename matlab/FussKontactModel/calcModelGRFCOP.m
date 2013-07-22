function grfCOP = calcModelGRFCOP(contactInfo)

heelF = contactInfo(1:3);
foreF = contactInfo(7:9);
heelCOP=contactInfo(4:6);
foreCOP=contactInfo(10:12);

grf = heelF+foreF;
grfM = sum(grf.*grf)^0.5;

cop = [];
%Set COP to the lowest patch's cop if the foot is in the air
if(heelCOP(3) < foreCOP(3))
   cop = heelCOP;
else
   cop = foreCOP; 
end

%If there's any ground reaction force, take the actual COP
if grfM > 0.001
cop = [ (heelF(3)*heelCOP(1)+foreF(3)*foreCOP(1))/grf(3); ...
        (heelF(3)*heelCOP(2)+foreF(3)*foreCOP(2))/grf(3); ...
        0];
end

grfCOP = [grf;cop];    