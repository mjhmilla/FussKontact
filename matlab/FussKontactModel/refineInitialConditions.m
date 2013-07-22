function ic = refineInitialConditions(ic0, expGRF, expCOP, vParams, vInput, vToe)

ic = ic0;

%Start the foot at rest;
ic(1) = 0;
ic(2) = 0;
ic(3) = 0;
ic(4) = 0;
ic(5) = 0;
ic(6) = 0;
ic(7) = 0;

%%
%Compute z, such that Fz matche experimental data
%%

%set the penetration depth to 1mm
contactInfo=calcContactForcePosition(0,ic,vParams,[0 0 0 0 0 0 0]');
gap = min([contactInfo(6),contactInfo(12)]);
ic(10) = ic(10)-gap-0.001;


options = optimset('Display','off','Diagnostics','off',...
                   'TolX',1e-12,'TolFun',1e-12, 'MaxIter', 1000);

errCFcn = @(x)calcContactForceError(x, ic, vParams, vToe, expGRF(1,:), expCOP(1,:));
[icDelta copErrSq] = fminunc(errCFcn,[0 0 0], options);
errf = errCFcn(icDelta);

ic(10) = ic(10) + icDelta(1)./1000; %z
ic(12) = ic(12) + icDelta(2)./1000; %zeta
ic(13) = ic(13) + icDelta(3)./1000; %eta

%%
%Compute vx,vy such that Fx and Fy match experimental data
%%
%errFFcn = @(x)calcContactForceError(x, ic, vParams, vToe, expGRF(1,:));
%[vxvyDelta fxfyErrSq] = fminunc(errFFcn,[0,0], options);
%ic(1) = ic(1) + vxvyDelta(1);
%ic(2) = ic(2) + vxvyDelta(2);

%%
%Algebra: Compute x, y position so that the COP matches: 
%         straight addition sub
%%
% contactInfo=calcContactForcePosition(0,ic,vParams,[0 0 0 0 0 0 0]');
% grfCOP = calcModelGRFCOP(contactInfo);
% grf = grfCOP(1:3);
% cop = grfCOP(4:6);
% 
% disp('refineInitialConditions: not adjusting COP anymore');
% %ic(8) = ic(8) + expCOP(1,1)-cop(1);
% %ic(9) = ic(9) + expCOP(1,2)-cop(2);
% 
% %Put the toe in the neural position if it is not under load
% foreF = contactInfo(7:9);
% if(sum(foreF.*foreF) < 0.1)
%    ic(4) = 0;
%    ic(11)= 0;
% end

%%
%Output results
%%
contactInfo=calcContactForcePosition(0,ic,vParams,[0 0 0 0 0 0 0]');
grfCOP = calcModelGRFCOP(contactInfo);
grf = grfCOP(1:3);
cop = grfCOP(4:6);
disp(sprintf('Mdl:Exp FxErr(%f) FyErr(%f) FzErr(%f)',...
    grf(1)-expGRF(1,1),grf(2)-expGRF(1,2),grf(3)-expGRF(1,3)));
    
disp(sprintf('        COPxErr(%f) COPyErr(%f)',...
    cop(1)-expCOP(1,1),cop(2)-expCOP(1,2)));
