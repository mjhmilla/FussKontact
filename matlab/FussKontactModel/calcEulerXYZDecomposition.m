function ZetaEtaXi = calcEulerXYZDecomposition(R)

% A    = [rM(i,1) rM(i,2) rM(i,3); ...
%         rM(i,4) rM(i,5) rM(i,6); ...
%         rM(i,7) rM(i,8) rM(i,9)];
% 
% %Extract angle 1
% ea313(i,1) =  atan2(A(3,1), A(3,2));
% 
% %Extract angle 2
% sinBeta = (A(1,3)*A(1,3) + A(2,3)*A(2,3))^0.5;
% ea313(i,2) = atan2(sinBeta, A(3,3));
% 
% %Form component rotation matricies 1 and 2, and compute rotation
% %matrix 3
% s1 = sin(ea313(i,1));
% c1 = cos(ea313(i,1));    
% r1 = [c1 -s1 0; s1 c1 0; 0 0 1];
% 
% s2 = sin(ea313(i,2));
% c2 = cos(ea313(i,2));
% r2 = [1 0 0; 0 c2 -s2; 0 s2 c2];
% 
% r3 = A*r1'*r2';
% ea313(i,3) = atan2( r3(2,1), r3(1,1));
% 
% 
% %%
% %Check
% %%
% eM = A - r3*r2*r1;
% err = norm(eM);
% if (err > rootEPS)
%    disp('EA313 decomposition went bad!!!') 
% end

