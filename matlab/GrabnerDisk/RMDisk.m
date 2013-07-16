function RM = calcRotationMatrix(xdot)

alpha = xdot(:,10);
beta = xdot(:,11);
zeta = xdot(:,12);


t1 = cos(zeta);
t2 = cos(beta);
t3 = sin(beta);
t4 = sin(alpha);
t5 = sin(zeta);
t6 = cos(alpha);
t7 = t1 .* t4;
t8 = t5 .* t6;
t9 = t1 .* t6;
t10 = t5 .* t4;
RM = [t1 * t2, t7 * t3 + t8, -t9 * t3 + t10,...
      -t5 * t2, -t10 * t3 + t9, t8 * t3 + t7,...
            t3, -t2 * t4, t2 * t6];
        
        
