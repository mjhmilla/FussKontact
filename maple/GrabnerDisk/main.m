mus = 1.0;
mud = 0.7;
stVel = 0.001;
dyVel = 0.01;
vel = [-dyVel*2: (dyVel*4)/999 :dyVel*2];

mu = zeros(size(vel));

for i=1:1:length(vel)
   mu(i) = calculateFrictionCoefficient(vel(i), mus,mud,stVel,dyVel); 
end

fig = figure;
plot(vel,mu);