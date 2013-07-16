GonthierSCM := proc(r1::array, r2::array, v1::array, v2::array, R1::array, R2::array, 
				w1::array, w2::array, leParams::array, aOut::array(float))
			
######################################################
#
#
#	Gonthier, McPhee, Lange, Piedboeuf Volumetric Sphere & Plane Contact Model
#		with viscous Coulomb friction (no velocity term)
#		- Matt Millard Jan 18 2007 -
#####################################################
#
#
#
#	[ 1, 2, 3,  4, 5, 6,  7, 8, 9,  10,  11,  12,  13,  14,  15, 16,17,18,  19,  20,   21,   22] 
#	[cx,cy,cz, nx,ny,nz, sx,sy,sz, tx1, ty1, tz1, tx2, ty2, tz2,  R, k, c, mus, mud, velS, velD]
#	
# Parameters:
# Geo:Location of Pad
# 	1-0. cx 	Centre of sphere, expressed in body-fixed frame 
# 	2-1. cy		Centre of sphere, expressed in body-fixed frame 
# 	3-2. cz	Centre of sphere, expressed in body-fixed frame 
# Geo:Normal of floor surface
# 	4-3. nx	Normal vector to surface, expressed in body-fixed frame 
# 	5-4. ny	Normal vector to surface, expressed in body-fixed frame 
# 	6-5. nz	Normal vector to surface, expressed in body-fixed frame 
# Geo: Location of Surface 
# 	7-6. sx	Point on surface, expressed in body-fixed frame 2 - x
# 	8-7. sy	Point on surface, expressed in body-fixed frame 2 - y
# 	9-8. sz	Point on surface, expressed in body-fixed frame 2 - z
## Geo: Tangent of floor surface (passing it in for computation purposes
#	10-9. tx1 Tangent vector to surface, expressed in a body-fixed frame
#	11-10. ty1 Tangent vector to surface, expressed in a body-fixed frame
#	12-11. tz1 Tangent vector to surface, expressed in a body-fixed frame
# Geo: Tangent of floor surface (passing it in for computation purposes
#	13-12. tx2 Tangent vector to surface, expressed in a body-fixed frame
#	14-13. ty2 Tangent vector to surface, expressed in a body-fixed frame
#	15-14. tz2 Tangent vector to surface, expressed in a body-fixed frame
# Geometry & Material Properties 
# 	16-15. R	Radius of sphere								
#  	17-16. k 	Average spring constant - use (1.2*body weight)/d for walking
# 	18-17. c	Average damping factor - tune to fit
# Friction Coefficients
# 	19-18. mus	Stiction friction value
# 	20-19. mud 	Dynamic friction value
# 	21-20. velS	Stiction Velocity
# 	22-21. velD	Dynamic Friction Velocity
				
#############################################
#
#		Create necessary Variables:
#
#############################################								
	local i::integer,	#for iteration
		iterMax::integer,
		height::float, 				#height of sphere above plane (in normal direction)
		x::float, 					#Penetration Depth
		pieVal::float,				#Pi
		muVel::float,				#Velocity used to solve for mu
		wNormalMag::float,			#Magnitude of the relative w of the bodies in the normal direction
		wTangentMag::float,			#Magnitude of the relative w of the bodies in the tangential directions
		velTangentMag::float,			#Magnitude of the tangential velocity;
		vcn::array(1..3),			#Normal velocity of the centroid of the interpenation volume.
		velRelative::array(1..3),		#Relative velocity of the two surfaces
		velNormal::array(1..3), 			#Normal velocity
		velTangent::array(1..3),		#Tangential velocity
		wRelative::array(1..3),			#Relative angular velocity of the two surfaces
		wTangent::array(1..3),			#Tangential angular velocity of the sphere
		wNormal::array(1..3),			#Normal angular velocity of the sphere.
		contactForceMag::float,
		contactForce::array(1..3),		#Contact force
		linearFrictionForce::array(1..3),	#Linear friction force magnitude
		linearFrictionMoment::array(1..3),	#Moment created by the tangential friction forces.
		spinFrictionTorque::array(1..3),	#Spinning friction torque magnitude
		rollingResistance::array(1..3),		#Rolling resistance
		interpenVolume::float,			#Interpenation volume
		interpenJ::array(1..2),			#Contact surface force weighted second moment of inertia 
								#  just xx, yy, zz no cross products due to symmetry.
								# 1:  Normal
								# 2: Tangential
		f::float,				#Centroid: F in Newton's method used to find centroid of the interpen.Vol
		df::float,				#Centroid: dF in Newton's method used to find the centroid of the interpen.Vol
		xc::float,				#Centroid: Used in the calculation of J
		pen::float,				#General = abs(x). Just so I don't forget to take the absolute value.
		
		muT::float,				#Tangential friction coefficient;
		muSpin::float,				#Spinning friction coefficient;
		radGy::float,				#Radius of gyration, for spinning friction.

		#For the sake of ease of writing, taking all leParams and assigning to 
		#appropriately named variables.
		
		sphereCenter::array(1..3),	#Sphere center 			
		floorNormal::array(1..3),	#Normal vector of the floor	
		floorPoint::array(1..3),	#A point on the floor.		
		R::float,			#Sphere radius		
		k::float,			#Contact stiffness coefficient	
		c::float,			#Damping coefficient		

		floorTangent1::array(1..3),	#Floor tangent 1		
		floorTangent2::array(1..3),	#Floor tangent 2		
		muS::float,			#Mu Stiction				
		muD::float,			#Mu Dynamic Friction		
		velS::float,			#Stiction transition velocity	
		velD::float,			#Dynamic friction transition vel
		delta::float,			#For the cubic step 
		pSphere::array(1..3), 		#Loc of center of sphere in inertial frame
		normal::array(1..3),		#Normal of contacting surfaces
		tangent1::array(1..3),		#Tangent1 of contact surfaces
		tangent2::array(1..3),		#Tangent2 of contacting surfaces
		pAction::array(1..3),		#Loc of contact point
		tmp1::array(1..3), 		#X product of contact point of body 1 (sphere) w.r.t. inertial frame
		tmp2::array(1..3): 		#X product of contact point of body 2 (floor) w.r.t. inertial frame

		#initialise arrays
		pSphere 	:= array(1..3):
		normal 		:= array(1..3):
		tangent1	:= array(1..3):
		tangent2	:= array(1..3):
		pAction 	:= array(1..3):
		tmp1 	:= array(1..3):
		tmp2 	:= array(1..3):	
		
		vcn			:= array(1..3):
		velRelative		:= array(1..3):
		velNormal 		:= array(1..3):		
		velTangent 		:= array(1..3):
		wRelative		:= array(1..3):
		wTangent 		:= array(1..3):			
		wNormal 		:= array(1..3):			
		contactForce 		:= array(1..3):		
		linearFrictionForce 	:= array(1..3):	
		spinFrictionTorque 	:= array(1..3):	
		rollingResistance 	:= array(1..3):	
		interpenJ 		:= array(1..2):		

		sphereCenter 	:= array(1..3):	
		floorNormal 	:= array(1..3):	
		floorPoint 		:= array(1..3):		
		floorTangent1 	:= array(1..3):	
		floorTangent2 	:= array(1..3):

		linearFrictionMoment:=array(1..3):
	
	for i from 1 to 12 do
		aOut[i] := 0:
	od:
#############################################
#
#		Set parameter lists up:
#
#############################################
	iterMax 	:=		10000: #Maximum number of iterations Newton's method is allowed to go for
	sphereCenter[1]:=		leParams[1]: 	#Sphere center 				leParams 1
	sphereCenter[2]:=		leParams[2]: 	#Sphere center 				leParams 2
	sphereCenter[3]:=		leParams[3]: 	#Sphere center 				leParams 3
	
	floorNormal[1] :=		leParams[4]:	#Normal vector of the floor		leParams 4
	floorNormal[2] :=		leParams[5]:	#Normal vector of the floor		leParams 5
	floorNormal[3] :=		leParams[6]:	#Normal vector of the floor		leParams 6
	
	floorPoint[1] :=		leParams[7]:		#A point on the floor.		leParams 7
	floorPoint[2] :=		leParams[8]:		#A point on the floor.		leParams 8
	floorPoint[3] :=		leParams[9]:		#A point on the floor.		leParams 9
							
	floorTangent1[1] := 		leParams[10]:	#Floor tangent 1			leParams 10
	floorTangent1[2] :=		leParams[11]:	#Floor tangent 1			leParams 11
	floorTangent1[3] :=		leParams[12]:	#Floor tangent 1			leParams 12
	
	floorTangent2[1] :=		leParams[13]:	#Floor tangnet 2			leParams 13
	floorTangent2[2] :=		leParams[14]:	#Floor tangnet 2			leParams 14
	floorTangent2[3] :=		leParams[15]:	#Floor tangnet 2			leParams 15
	
			R:=		leParams[16]:	#Sphere radius				leParams 16
			k:=		leParams[17]:	#Volumetic Contact stiffness 		leParams 17
			c:=		leParams[18]:	#Damping coefficient			leParams 18
	
			muS:=		leParams[19]:	#Mu Stiction				leParams 19
			muD:=		leParams[20]:	#Mu Dynamic Friction			leParams 20
			velS:=		leParams[21]: 	#Stiction transition velocity		leParams 21
			velD:=		leParams[22]:	#Dynamic friction transition vel	leParams 22
	pieVal := 3.1415926535897932384626433832795:
#############################################
#
#		Calculate penetration depth: Contact or No contact?
#
#############################################		
	
	# Convert to global reference frame using rotation transformation matrices R1, R2
	pSphere[1] := R1[1]*sphereCenter[1] + R1[4]*sphereCenter[2] + R1[7]*sphereCenter[3]:
	pSphere[2] := R1[2]*sphereCenter[1] + R1[5]*sphereCenter[2] + R1[8]*sphereCenter[3]:
	pSphere[3] := R1[3]*sphereCenter[1] + R1[6]*sphereCenter[2] + R1[9]*sphereCenter[3]:
		
	#Normal of the floor in an inertial frame reference. Included here to make it easy to expand this code to
	#other things.
	normal[1] := R2[1]*floorNormal[1] + R2[4]*floorNormal[2] + R2[7]*floorNormal[3]:
	normal[2] := R2[2]*floorNormal[1] + R2[5]*floorNormal[2] + R2[8]*floorNormal[3]:
	normal[3] := R2[3]*floorNormal[1] + R2[6]*floorNormal[2] + R2[9]*floorNormal[3]:
	
	
	
	# height of surface from origin of body-fixed frame 2
	# project point on surface onto normal vector
	# reference frame irrelevent; as long as consistent
	height := floorNormal[1]*floorPoint[1] + floorNormal[2]*floorPoint[2] + floorNormal[3]*floorPoint[3]:
	
	# Point of action
	
	pAction[1] := r1[1] + pSphere[1] - R*normal[1]:
	pAction[2] := r1[2] + pSphere[2] - R*normal[2]:
	pAction[3] := r1[3] + pSphere[3] - R*normal[3]:
	
	# Depth of penetration
	x := (pAction[1]-r2[1])*normal[1] + (pAction[2]-r2[2])*normal[2] + (pAction[3]-r2[3])*normal[3] - height:
	
#############################################
#
#		If surfaces are in contact:
#
#############################################			
	if (x < 0 and R > 0) then
		pen:=abs(x):
		
		if (pen > 2*R) then
			pen := 2*R:
		fi:	
		#Recalculate pAction for calculations involving the tangential velocity of
		#the contact surface. pAction is now the centroid of the contact surface.
		#rather than a point on the contact sphere's surface.
		pAction[1] := r1[1] + pSphere[1] - (R-pen)*normal[1]:
		pAction[2] := r1[2] + pSphere[2] - (R-pen)*normal[2]:
		pAction[3] := r1[3] + pSphere[3] - (R-pen)*normal[3]:
	#############################################
	#
	#		Calculate tangential, normal linear and rotational velocities, 
	#		volume of penetration, moment of inertia of the contact surface
	#
	#############################################	
		tangent1[1] := R2[1]*floorTangent1[1] + R2[4]*floorTangent1[2] + R2[7]*floorTangent1[3]:
		tangent1[2] := R2[2]*floorTangent1[1] + R2[5]*floorTangent1[2] + R2[8]*floorTangent1[3]:
		tangent1[3] := R2[3]*floorTangent1[1] + R2[6]*floorTangent1[2] + R2[9]*floorTangent1[3]:
		
		tangent2[1] := R2[1]*floorTangent2[1] + R2[4]*floorTangent2[2] + R2[7]*floorTangent2[3]:
		tangent2[2] := R2[2]*floorTangent2[1] + R2[5]*floorTangent2[2] + R2[8]*floorTangent2[3]:
		tangent2[3] := R2[3]*floorTangent2[1] + R2[6]*floorTangent2[2] + R2[9]*floorTangent2[3]:
	
		# Cross product of angular velocity of sphere with radial arm of point of action
		tmp1[1] := w1[2]*(pAction[3] - r1[3]) - w1[3]*(pAction[2] - r1[2]):
		tmp1[2] := w1[3]*(pAction[1] - r1[1]) - w1[1]*(pAction[3] - r1[3]):
		tmp1[3] := w1[1]*(pAction[2] - r1[2]) - w1[2]*(pAction[1] - r1[1]):
	
		# Cross product of angular velocity of surface with radial arm of point of action
		tmp2[1] := w2[2]*(pAction[3] - r2[3]) - w2[3]*(pAction[2] - r2[2]):
		tmp2[2] := w2[3]*(pAction[1] - r2[1]) - w2[1]*(pAction[3] - r2[3]):
		tmp2[3] := w2[1]*(pAction[2] - r2[2]) - w2[2]*(pAction[1] - r2[1]):
	
		# RELATIVE VELOCITY at POINT OF ACTION, projected onto the normal
		velRelative[1] := (v1[1] + tmp1[1] - v2[1] - tmp2[1]):
		velRelative[2] := (v1[2] + tmp1[2] - v2[2] - tmp2[2]):
		velRelative[3] := (v1[3] + tmp1[3] - v2[3] - tmp2[3]):
		
		velNormal[1] := velRelative[1]*normal[1]:
		velNormal[2] := velRelative[2]*normal[2]:
		velNormal[3] := velRelative[3]*normal[3]:
		
		velTangent[1] := velRelative[1]*(tangent1[1] + tangent2[1]):
		velTangent[2] := velRelative[2]*(tangent1[2] + tangent2[2]):
		velTangent[3] := velRelative[3]*(tangent1[3] + tangent2[3]):
	
		wRelative[1] := w1[1]-w2[1]:
		wRelative[2] := w1[2]-w2[2]:
		wRelative[3] := w1[3]-w2[3]:
		
		wTangent[1] := wRelative[1]*(tangent1[1] + tangent2[1]):
		wTangent[2] := wRelative[2]*(tangent1[2] + tangent2[2]):
		wTangent[3] := wRelative[3]*(tangent1[3] + tangent2[3]):
		
		wNormal[1] := wRelative[1]*normal[1]:
		wNormal[2] := wRelative[2]*normal[2]:
		wNormal[3] := wRelative[3]*normal[3]:
		
		velTangentMag 	:= (velTangent[1]*velTangent[1]+velTangent[2]*velTangent[2]+velTangent[3]*velTangent[3])^0.5:
		wNormalMag 	:= (wNormal[1]*wNormal[1] + wNormal[2]*wNormal[2] + wNormal[3]*wNormal[3])^0.5:
		wTangentMag 	:= (wTangent[1]*wTangent[1] + wTangent[2]*wTangent[2] + wTangent[3]*wTangent[3])^0.5:
		#Calculate the interpenation volume:
			interpenVolume := (1/3)*pieVal*pen*pen*(3*R-pen):
		#Calculate the centroid of the interpenation volume
			#I found this converges so quickly that Newton's method is
			#actually as fast as a good spline fit. So I use Newton's method

			xc	:=0.5*pen:
			df	:=0.0:
			f	:=interpenVolume:
			i := 0:
			while (abs(f) > 0.0001*interpenVolume and i < iterMax) do
				f := (1/3)*pieVal*xc*xc*(3*R-xc) - 0.5*interpenVolume:
				df:= (1/3)*pieVal*xc*(6*R-3*xc):
				xc:= xc-f/df:
				i := i+1:
				
			end do:

			if(i >= iterMax) then
				pen:= 0.001*R:
				xc := 0.5*pen:
			fi:
		 # J in NormalDir with a unit value for rho, material density as specified by Gonthier
			interpenJ[1] :=(1/30)*1*pieVal*(pen*pen*pen)*(3*(pen*pen)-15*R*pen+20*R*R):
		 # J in Tangential with a unit value for rho, material density as specified by Gonthier
			interpenJ[2] :=(1/60)*1*pieVal*(pen*pen)*(-9*(pen*pen*pen)+15*R*(pen*pen)+30*(pen*pen)*xc-20*pen*xc*xc+20*R*R*pen-80*pen*R*xc+60*xc*xc*R):			
		#############################################
		#
		#		Calculate normal force, rolling resistance, linear and rotational frictional forces
		#
		#############################################	
		
		#Really need vcn, the velocity of the centroid. I originally erroroneously
		#scaled vcn by xc/pen. It was a mistaken, so I've corrected it.
		vcn[1]	:= -velNormal[1]:#*(xc/(pen)): 
		vcn[2]	:= -velNormal[2]:#*(xc/(pen)):
		vcn[3]	:= -velNormal[3]:#*(xc/(pen)):
		
		#CONTACT FORCE	
		#
		
		 if(c*(vcn[1]+vcn[2]+vcn[3]) < -1) then
			contactForce[1]:= 0.0:
			contactForce[2]:= 0.0:
			contactForce[3]:= 0.0:
		else
			contactForce[1]:= k*interpenVolume*(1 + vcn[1]*c)*normal[1]:
			contactForce[2]:= k*interpenVolume*(1 + vcn[2]*c)*normal[2]:
			contactForce[3]:= k*interpenVolume*(1 + vcn[3]*c)*normal[3]:
		fi:
			
		contactForceMag:=(contactForce[1]*contactForce[1] + contactForce[2]*contactForce[2] + contactForce[3]*contactForce[3])^0.5:
			
			

		#ROLLING RESISTANCE
		
		if(wTangentMag > 0) then
			rollingResistance[1]	:=-(k)*c*interpenJ[2]*wTangent[1]:
			rollingResistance[2]	:=-(k)*c*interpenJ[2]*wTangent[2]:
			rollingResistance[3]	:=-(k)*c*interpenJ[2]*wTangent[3]:
		else
			rollingResistance[1]	:= 0.0:
			rollingResistance[2]	:= 0.0:
			rollingResistance[3]	:= 0.0:
		
		fi:
		
		#TANGENTIAL (linear) FRICTION: 
		#			   
		#	Calculate mu using a Coulomb friction model
			muVel := velTangentMag:
			if (0.0e0 < muVel) then

			    #In dynamic region muT
			    if (muVel > velD) then
					muT := muD:
				    
			    #Between dynamic and static		    
			    elif (muVel > velS and muVel < velD) then
				    delta 	:= (muVel - velS)/(velD-velS):
				    muT 	:= muS +(muD-muS)*delta*delta*(3.0-2.0*delta):
			    #In static region
			    else
				    delta	:= (muVel + velS)/(2.0*velS):
				    muT 	:= -muS + (2*muS)*delta*delta*(3.0-2.0*delta):
			    fi:

			fi:	
		
			#Viscous friction model, un comment lines below to get a dry friction model
			if(velTangentMag > 0) then
				linearFrictionForce[1]:= -muT*contactForceMag*velTangent[1]:#/(velTangentMag+1e-6):	
				linearFrictionForce[2]:= -muT*contactForceMag*velTangent[2]:#/(velTangentMag+1e-6):
				linearFrictionForce[3]:= -muT*contactForceMag*velTangent[3]:#/(velTangentMag+1e-6):
			else
			   	linearFrictionForce[1] := 0.0:
			   	linearFrictionForce[2] := 0.0:
			   	linearFrictionForce[3] := 0.0:
			fi:
		
			
		
		#MOMENT CREATED BY TANGENTIAL FRICTION:
			linearFrictionMoment[1] := (normal[3]*linearFrictionForce[2] - normal[2]*linearFrictionForce[3])*(R-pen):
			linearFrictionMoment[2] := (normal[1]*linearFrictionForce[3] - normal[3]*linearFrictionForce[1])*(R-pen):
  				linearFrictionMoment[3] := (normal[2]*linearFrictionForce[1] - normal[1]*linearFrictionForce[2])*(R-pen):
		
		#SPINNING FRICTIONAL TORQUE
			radGy := sqrt(interpenJ[1]/interpenVolume):
			muVel := radGy*wNormalMag:


			if (0.0e0 < muVel) then

			    #In dynamic region muT
			    if (muVel > velD) then
					muSpin := muD:

			    #Between dynamic and static		    
			    elif (muVel > velS and muVel < velD) then
				    delta 	:= (muVel - velS)/(velD-velS):
				    muSpin 	:= muS +(muD-muS)*delta*delta*(3.0-2.0*delta):
			    #In static region
			    else
				    delta	:= (muVel + velS)/(2.0*velS):
				    muSpin 	:= -muS + (2*muS)*delta*delta*(3.0-2.0*delta):
			    fi:

			fi:

		
			if(wNormalMag > 0.0) then
			
				spinFrictionTorque[1]:= -muSpin*contactForceMag*radGy*wNormal[1]:#/(wNormalMag + 1e-6):
				spinFrictionTorque[2]:= -muSpin*contactForceMag*radGy*wNormal[2]:#/(wNormalMag + 1e-6):
				spinFrictionTorque[3]:= -muSpin*contactForceMag*radGy*wNormal[3]:#/(wNormalMag + 1e-6):				
			
			else
				spinFrictionTorque[1] := 0.0:
				spinFrictionTorque[2] := 0.0:
				spinFrictionTorque[3] := 0.0:
			fi:
		
			
		
		
		# NORMAL FORCE ON SPHERE
		#Contact, linear friction
		aOut[1] := contactForce[1] + linearFrictionForce[1]:
		aOut[2] := contactForce[2] + linearFrictionForce[2]:
		aOut[3] := contactForce[3] + linearFrictionForce[3]:
		#Rolling resistance, spinning friction
		aOut[4] := rollingResistance[1] + spinFrictionTorque[1] + linearFrictionMoment[1]:
		aOut[5] := rollingResistance[2] + spinFrictionTorque[2] + linearFrictionMoment[2]:
		aOut[6] := rollingResistance[3] + spinFrictionTorque[3] + linearFrictionMoment[3]:
		# FORCE ON SURFACE
		aOut[7] := -aOut[1]:
		aOut[8] := -aOut[2]:
		aOut[9] := -aOut[3]:
		aOut[10] := -aOut[4]:
		aOut[11] := -aOut[5]:
		aOut[12] := -aOut[6]:
	fi:
	return 0.0:
end proc: