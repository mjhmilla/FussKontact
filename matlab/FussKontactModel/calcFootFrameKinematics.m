function kin = calcFootFrameKinematics(state, vParams, vToe)

%Parameters
buildParamVariableList;

%State Information
x    = state(8);
y    = state(9);
z    = state(10);
th   = state(11);
zeta = state(12);
eta  = state(13);
xi   = state(14);


%Compute all kinematic positions and frame locations.
calcPosVecsRotMatrices;

kin.ankle.r = rAnkle;
kin.ankle.RM= RMAnkle;
kin.K1com.r = rK1com;
kin.K1com.RM= RMK1com;
kin.K1c.r   = rK1c;
kin.K1c.RM  = RMK1c;
kin.K2com.r = rK2com;
kin.K2com.RM= RMK2com;
kin.Heel.r = rHeel;
kin.Heel.RM=RMHeel;
kin.Fore.r = rForeFoot;
kin.Fore.RM=RMForeFoot;
