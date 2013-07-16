unprotect(MapleSim:-Multibody);
sRootDir := "C:/mjhmilla/library/Maple/MapleSim_DFP/modules/":

#Construct and initialise the various modules
kernelopts(opaquemodules=false):

#Construct and initialise the library manager
MapleSim:-Multibody:-mLibMan := MapleSim:-Multibody:-mConstructors:-CmLibMan():

rsTemplateList := table([]):

rsTemplateList["m3DContact"] := cat(sRootDir, "C:/mjhmilla/library/Maple/MapleSim_DFP/library/mechanical/other/m3DContact.mpl"):

#rsTemplateList["mTire"] := cat(sRootDir,"../library/mechanical/other/Tire_Modified_NoExpMan.mpl"):
#rsTemplateList["mTire"] := cat(sRootDir,"../library/mechanical/other/Tire_Road.mpl"):
rsTemplateList["mTireKinematics"] := cat(sRootDir,"../library/mechanical/tires/TireKinematics.mpl"):
rsTemplateList["mTireForces"] := cat(sRootDir,"../library/mechanical/tires/TireForces.mpl"):
rsTemplateList["mTire"] := cat(sRootDir,"../library/mechanical/tires/Tire.mpl"):
rsTemplateList["mBushing"] := cat(sRootDir,"../library/mechanical/other/Bushing_Modified.mpl"):
rsTemplateList["mGenMotDrv"] := cat(sRootDir,"../library/mechanical/drivers/MotionDriver.mpl"):
rsTemplateList["mTSDA"] := cat(sRootDir,"../library/mechanical/drivers/TSDA.mpl"):
rsTemplateList["mRSDA"] := cat(sRootDir,"../library/mechanical/drivers/RSDA.mpl"):
rsTemplateList["mForce"] := cat(sRootDir,"../library/mechanical/drivers/Force.mpl"):
rsTemplateList["mMoment"] := cat(sRootDir,"../library/mechanical/drivers/Moment.mpl"):
rsTemplateList["mRigidBody"] := cat(sRootDir,"../library/mechanical/bodies/RigidBody.mpl"):
rsTemplateList["mRigidBodyFrame"] := cat(sRootDir,"../library/mechanical/bodies/RigidArm.mpl"):
rsTemplateList["mFlexBeam"] := cat(sRootDir,"../library/mechanical/bodies/FlexBeam.mpl"):
rsTemplateList["mFlexBeamFrame"] := cat(sRootDir,"../library/mechanical/bodies/FlexArm.mpl"):
rsTemplateList["mWeldJt"] := cat(sRootDir,"../library/mechanical/joints/Weld.mpl"):
rsTemplateList["mRevJt"] := cat(sRootDir,"../library/mechanical/joints/Revolute.mpl"):
rsTemplateList["mUnivJt"] := cat(sRootDir,"../library/mechanical/joints/Universal.mpl"):
rsTemplateList["mSphJt"] := cat(sRootDir,"../library/mechanical/joints/Spherical.mpl"):
rsTemplateList["mPrisJt"] := cat(sRootDir,"../library/mechanical/joints/Prismatic.mpl"):
rsTemplateList["mPlanJt"] := cat(sRootDir,"../library/mechanical/joints/Planar.mpl"):
rsTemplateList["mXYZTranJt"] := cat(sRootDir,"../library/mechanical/joints/XYZTran.mpl"):
rsTemplateList["mCylJt"] := cat(sRootDir,"../library/mechanical/joints/Cylindrical.mpl"):
rsTemplateList["mFreeJt"] := cat(sRootDir,"../library/mechanical/joints/FreeJoint.mpl"):


#Load general template
read cat(sRootDir, "/DFPEngine/GeneralTemplateFuncs.mpl"):

for iTemplate from 1 to nops([indices(rsTemplateList)]) do
 MapleSim:-Multibody:-mLibMan:-LoadTemplate(op(indices(rsTemplateList)[iTemplate]),
 							convert(rsTemplateList,list)[iTemplate], "FILE"):
od:

MapleSim:-Multibody:-mLibMan:-rpTemplateFunctions := copy(MapleSim:-Multibody:-mConstants:-rSubFunc):


MapleSim:-Multibody:-mMSimInterface := MapleSim:-Multibody:-mConstructors:-CmDMInterface():
MapleSim:-Multibody:-mCodeGen := MapleSim:-Multibody:-mConstructors:-CmCodeGen(MapleSim:-Multibody:-mLibMan):
MapleSim:-Multibody:-mStoreFileMan := MapleSim:-Multibody:-mConstructors:-CmStoreFileMan():
MapleSim:-Multibody:-mProjectBuilder := MapleSim:-Multibody:-mConstructors:-CmProjectBuilder():
MapleSim:-Multibody:-BuildEQs := MapleSim:-Multibody:-mConstructors:-CmEQBuilder():
MapleSim:-Multibody:-BuildSimulation := MapleSim:-Multibody:-mConstructors:-CmFwdDynSimBuilder():
MapleSim:-Multibody:-BuildExpression := MapleSim:-Multibody:-mConstructors:-CmExpressionBuilder():
MapleSim:-Multibody:-BuildSimCode := MapleSim:-Multibody:-mConstructors:-CmSimCodeBuilder(MapleSim:-Multibody:-mCodeGen):
MapleSim:-Multibody:-BuildPlot := MapleSim:-Multibody:-mConstructors:-CmPlotBuilder():
kernelopts(opaquemodules=true):
