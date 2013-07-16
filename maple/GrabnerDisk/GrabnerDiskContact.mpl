################################################################################
#
#===============================================================================
#   Title          :   Grabner Disk Contact with Coulomb Friction
#   Contact Pair   :   Between a Disk and a Plane. 
#                      The local disk normal is the z axis.
#   Contact Model  :   Hunt-Crossley, that is f = k*x^p*(1+ c*dx/dt)
#   Contact Model  :   Hunt-Crossley, that is f = k*x^p*(1+ c*dx/dt)
#   Friction Model :   Coulomb friction, smooth spline interpolated
#   Feature        :   Disk radius changes with angle between disk
#                      normal and plane normal vectors
#
#   @author     Matthew Millard
#   @date       June 30, 2013        
#===============================================================================
#   Notation Conventions
#    Frame Number
#       0   :   Ground frame
#       1   :   Plane frame
#       2   :   Disk Frame
#
#   Kinematics
#   r   :   position vector
#   R   :   rotation matrix
#   v   :   velocity vector
#   w   :   angular velocity
#   
#   Points
#   O: origin of frame 1, attached to the plane
#   M: origin of frame 2, attached to the disk
#   P: point on the disk in frame 2 that is the closest distance to 
#      the plane
#   Q: point of intersection between the disk and plane that
#      also intersects the vector from M to P.
#
#===============================================================================
#   @param r20  position of disk in the global reference frame
#   @param r10  position of plane in the global reference frame
#
#   @param v20   velocity of disk in the global reference frame
#   @param v10   velocity of plane in the global reference frame
#
#   @param R20   rotation matrix from global reference frame to the disk frame
#   @param R10   rotation matrix from global reference frame to the plane frame
#
#   @param w20   angular velocity of disk in the global reference frame
#   @param w10   angular velocity of plane in the global reference frame
#
#   @param leParams: [sx,sy,sz, nx,ny,nz, t1x,t1y,t1z, t2x,t2y,t2z
#                      r,a,k,p,c,mus,mud,stVel,dyVel]
#       1idx 0idx
#       1  - 0.   sx: Point on surface, expressed in body-fixed frame (m)
#       2  - 1.   sy: Point on surface, expressed in body-fixed frame (m)
#       3  - 2.   sz: Point on surface, expressed in body-fixed frame (m)
#
#       4  - 3.   nx: Normal vector to surface, expressed in body-fixed frame
#       5  - 4.   ny: Normal vector to surface, expressed in body-fixed frame
#       6  - 5.   nz: Normal vector to surface, expressed in body-fixed frame
#
#       7  - 6.   t1x: tangential surface vector 1
#       8  - 7.   t1y: tangential surface vector 1
#       9  - 8.   t1z: tangential surface vector 1
#
#       10 - 9.   t2x: tangential surface vector 2
#       11 - 10.  t2y: tangential surface vector 2
#       12 - 11.  t2z: tangential surface vector 2
#
#       13 - 12.  r:    radius of the disk (m)
#       14 - 13   a:    Grabner disk rounding factor   
#                       (Variable C in Eqn 4.21 of Grabner's thesis)
#
#       15 - 14.  k:    spring constant (N/m)
#       16 - 15.  p:    force law exponent 
#       17 - 16.  c:    damping factor
#
#       18 - 17.  mus:  coefficient of static friction.
#       19 - 18.  mud:  coefficient of dynamic friction
#       20 - 19.  stVel stiction transition velocity
#       21 - 20.  dyVel dynamic transition velocity
#
#   @params aOut: 12 element array of a concatination of the forces and
#                 torques applied to body 1, and then to body 2
#===============================================================================
#  References
#
#   Grabner, Gerald (2003). Efficient and Reliable Multibody Simulation of
#   Regularized Impacts between Elementary Contact Pairs. PhD Thesis, 
#   Instit fur Mechanik und Getriebelehre, Technische Universitat Graz.
#
#
################################################################################
#rSubFunc["m3DContact","GrabnerDiskContact"] := proc(r20::array, r10::array,
#                                                    v20::array, v10::array, 
#                                                    R20::array, R10::array, 
#                                                    w20::array, w10::array, 
#                                               leParams::array, aOut::array)
GrabnerDiskContact := proc(r20::array, r10::array,
                            v20::array, v10::array, 
                            R20::array, R10::array, 
                            w20::array, w10::array, 
                       leParams::array, aOut::array)

    local 
          sx::float,   
          sy::float,   
          sz::float,             
            nx::float,   
            ny::float,   
            nz::float,   
          t1x::float,  
          t1y::float,  
          t1z::float,            
            t2x::float,  
            t2y::float,  
            t2z::float,  
          r::float,    
          a::float,              
            k::float,    
            p::float,    
            c::float,    
          mus::float,  
          mud::float,  
          stVel::float,
          dyVel::float,    
            alpha::float,
            radius::float,
          tmpf1::float,  
          tmpf2::float,
          epsRoot::float,
          gN::float,
          i::integer,
          x::float,
          velNormal::float,
          mu::float,
          forceNormal::float,
          velTangentMag::float,
          velEps::float,
          forceTangentMag::float,
          daFactor::float,
          delta::float,
                    
          rO::array(1..3),
          rM::array(1..3),
          rMO::array(1..3),
          rPM::array(1..3),
          rPMuv::array(1..3),
          rPMO::array(1..3),
          rQO::array(1..3),
          rQM::array(1..3),          
          nC::array(1..3),
          nP::array(1..3),
          tP1::array(1..3),
          tP2::array(1..3),
          vQ10::array(1..3),
          vQ20::array(1..3),
          vP10::array(1..3),
          vP20::array(1..3),
          vPQ::array(1..3),
          velTangent::array(1..3),
          forceTangentUV::array(1..3),
          tmp1::array(1..3),
          tmp2::array(1..3):
          
        
        for i from 1 to 12 do
            aOut[i] := 0:
        od:

        epsRoot := 1.490116119384766e-008:
        
        #Put parameter vector into local variables
          sx     := leParams[1]:
          sy     := leParams[2]:
          sz     := leParams[3]:         
            nx   := leParams[4]:
            ny   := leParams[5]:
            nz   := leParams[6]:         
          t1x    := leParams[7]:
          t1y    := leParams[8]:
          t1z    := leParams[9]:         
            t2x  := leParams[10]:
            t2y  := leParams[11]:
            t2z  := leParams[12]:         
          r      := leParams[13]:
          a      := leParams[14]:         
            k    := leParams[15]:
            p    := leParams[16]:
            c    := leParams[17]:         
          mus    := leParams[18]:
          mud    := leParams[19]:
          stVel  := leParams[20]:
          dyVel  := leParams[21]:
        #printf("sx, sy, sz  : %f,%f,%f \n",sx,sy,sz);
        #printf("nx, ny, nz  : %f,%f,%f \n",nx,ny,nz);
        #printf("t1x,t1y,t1z : %f,%f,%f \n",t1x,t1y,t1z);
        #printf("t2x,t2y,t2z : %f,%f,%f \n",t2x,t2y,t2z);
        #printf("  r,  a,    : %f,%f \n",r,a);
        #printf("  k,  p,  c : %f,%f,%f \n",k,p,c);
        #printf("mus,mud     : %f,%f \n",mus,mud);
        #printf("stVel,dyVel : %f,%f \n\n",stVel,dyVel);
  
  
        #surface normal
          nP[1]  := R10[1]*nx  + R10[2]*ny  + R10[3]*nz:
          nP[2]  := R10[4]*nx  + R10[5]*ny  + R10[6]*nz:
          nP[3]  := R10[7]*nx  + R10[8]*ny  + R10[9]*nz:
          #printf("nP: %f,%f,%f\n",nP[1],nP[2],nP[3]);
          
        #Surface tangent vectors        
          tP1[1] := R10[1]*t1x + R10[2]*t1y + R10[3]*t1z:
          tP1[2] := R10[4]*t1x + R10[5]*t1y + R10[6]*t1z:
          tP1[3] := R10[7]*t1x + R10[8]*t1y + R10[9]*t1z:
          #printf("tP1: %f,%f,%f\n",tP1[1],tP1[2],tP1[3]);
          
          
          tP2[1] := R10[1]*t2x + R10[2]*t2y + R10[3]*t2z:
          tP2[2] := R10[4]*t2x + R10[5]*t2y + R10[6]*t2z:
          tP2[3] := R10[7]*t2x + R10[8]*t2y + R10[9]*t2z:
          #printf("tP2: %f,%f,%f\n",tP2[1],tP2[2],tP2[3]);          
                    
        #Origin of surface frame 
          rO[1] := R10[1]*sx + R10[2]*sy + R10[3]*sz + r10[1]:
          rO[2] := R10[4]*sx + R10[5]*sy + R10[6]*sz + r10[2]:
          rO[3] := R10[7]*sx + R10[8]*sy + R10[9]*sz + r10[3]:   
          #printf("r0: %f,%f,%f\n\n",rO[1],rO[2],rO[3]);
          
        #Origin of disk frame
          rM[1] := r20[1]:
          rM[2] := r20[2]:
          rM[3] := r20[3]:
          #printf("rM: %f,%f,%f\n",rM[1],rM[2],rM[3]);
          
        #vector between the origin of the disk frame and the plane   
          rMO[1] := (rM[1] - rO[1]): 
          rMO[2] := (rM[2] - rO[2]): 
          rMO[3] := (rM[3] - rO[3]): 
          #printf("rMO: %f,%f,%f\n\n",rMO[1],rMO[2],rMO[3]);
        #vector between the origin of the disk, and the lowest point on its 
        #edge P
        # 1. disk normal (local disk z axis)
          nC[1] := R20[7]:
          nC[2] := R20[8]:
          nC[3] := R20[9]:
          #printf("nC: %f,%f,%f\n",nC[1],nC[2],nC[3]);
          
        # 2. From Grabner's thesis eqn 4.18        
          # 2a. First cross product c = a x b, tmp1 = nc x np
          tmp1[1] :=  (nC[2]*nP[3] - nC[3]*nP[2]): # a[2]*b[3] - a[3]*b[2] 
          tmp1[2] := -(nC[1]*nP[3] - nC[3]*nP[1]): #-a[1]*b[3] + a[3]*b[1]
          tmp1[3] :=  (nC[1]*nP[2] - nC[2]*nP[1]): # a[1]*b[2] - a[2]*b[1]
          #printf("nC x nP: %f,%f,%f\n",tmp1[1],tmp1[2],tmp1[3]);
        
          # 2b. Second cross product c = a x b, tmp2 = nc x tmp1
          tmp2[1] :=  (nC[2]*tmp1[3] - nC[3]*tmp1[2]): # a[2]*b[3] - a[3]*b[2] 
          tmp2[2] := -(nC[1]*tmp1[3] - nC[3]*tmp1[1]): #-a[1]*b[3] + a[3]*b[1]
          tmp2[3] :=  (nC[1]*tmp1[2] - nC[2]*tmp1[1]): # a[1]*b[2] - a[2]*b[1]
          #printf("nC x nC x nP: %f,%f,%f\n\n",tmp2[1],tmp2[2],tmp2[3]);
         
          #Pick the low end, not the high end of the disk
          tmpf1 := tmp2[1]*nP[1] + tmp2[2]*nP[2] + tmp2[3]*nP[3]:
          if (tmpf1 > 0) then
            tmp2[1] := -1.0*tmp2[1]:
            tmp2[2] := -1.0*tmp2[2]:
            tmp2[3] := -1.0*tmp2[3]:
          fi:
         
          # 2c. Magnitude
          tmpf1:= (tmp1[1]*tmp1[1] + tmp1[2]*tmp1[2] + tmp1[3]*tmp1[3])^0.5:
            if (tmpf1 < epsRoot) then
                tmpf1 := epsRoot:
            fi:
          #printf("tmpf1: %f\n",tmpf1);
                     
          #3. Disk radius
          #3a. Acute angle between disk normal and surface normal
          alpha := arccos( abs( nP[1]*nC[1] + nP[2]*nC[2] + nP[3]*nC[3] )):
          #printf("alpha : %f\n", alpha);  
            
          #3b. Disk radius
          radius := r*(1 - exp(-a*sin(alpha)) ):
          #printf("radius: %f\n", radius);
          
          #4. Vector rPM
          rPMuv[1] := tmp2[1]/tmpf1:  
          rPMuv[2] := tmp2[2]/tmpf1:  
          rPMuv[3] := tmp2[3]/tmpf1:  
          #printf("rPMuv: %f,%f,%f\n", rPMuv[1],rPMuv[2],rPMuv[3]);          
          
          rPM[1] := radius*rPMuv[1]:  
          rPM[2] := radius*rPMuv[2]:
          rPM[3] := radius*rPMuv[3]:
          #printf("rPM: %f,%f,%f\n", rPM[1],rPM[2],rPM[3]);
          
        #Compute the shortest distance between vector rPO and the plane
          rPMO[1] := rMO[1] + rPM[1]:
          rPMO[2] := rMO[2] + rPM[2]:
          rPMO[3] := rMO[3] + rPM[3]:
          #printf("rPMO: %f,%f,%f\n\n", rPMO[1],rPMO[2],rPMO[3]);
          
          gN      := rPMO[1]*nP[1] + rPMO[2]*nP[2] + rPMO[3]*nP[3]:           
          #printf("gN: %f\n\n", gN);
          
        if (gN < 0.0) then
        
        ########################################################################
        # Compute the normal and tangential velocities between point P
        # on the disk and its projection on the plane
        ########################################################################
        
          #Calculate the velocity of point P on the disk. Note point P
          #is a point fixed on the disk, it is not the center of 
          #pressure.
          vP20[1] :=v20[1]+w20[2]*rPM[3]-w20[3]*rPM[2]:# a[2]*b[3]-a[3]*b[2]
          vP20[2] :=v20[2]-w20[1]*rPM[3]+w20[3]*rPM[1]:#-a[1]*b[3]+a[3]*b[1]
          vP20[3] :=v20[3]+w20[1]*rPM[2]-w20[2]*rPM[1]:# a[1]*b[2]-a[2]*b[1]
          #printf("vP20: %f,%f,%f\n", vP20[1],vP20[2],vP20[3]):
          
          #Calculate the vector rQO which spans from the origin of the
          #plane frame to point Q, the intersection point of vector rPM 
          #with the plane
          tmpf1  := rPMO[1]*tP1[1] + rPMO[2]*tP1[2] + rPMO[3]*tP1[3]:
          tmpf2  := rPMO[1]*tP2[1] + rPMO[2]*tP2[2] + rPMO[3]*tP2[3]:
          rQO[1] :=   tmpf1*tP1[1] + tmpf2*tP2[1]:
          rQO[2] :=   tmpf1*tP1[2] + tmpf2*tP2[2]:
          rQO[3] :=   tmpf1*tP1[3] + tmpf2*tP2[3]:                  
          #printf("rQO: %f,%f,%f\n", rQO[1],rQO[2],rQO[3]):          
                    
          #Calculate the velocity of point Q on the plane. 
          #Again, note that Q is a point fixed on the 
          #plane. It is not the center of pressure.
          vQ10[1] :=v10[1]+w10[2]*rQO[3]-w10[3]*rQO[2]:# a[2]*b[3]-a[3]*b[2]
          vQ10[2] :=v10[2]-w10[1]*rQO[3]+w10[3]*rQO[1]:#-a[1]*b[3]+a[3]*b[1]
          vQ10[3] :=v10[3]+w10[1]*rQO[2]-w10[2]*rQO[1]:# a[1]*b[2]-a[2]*b[1]
          #printf("vQ10: %f,%f,%f\n", vQ10[1],vQ10[2],vQ10[3]):
          
          #Calculate the difference in velocity of these two points
          #vPQ::array(1..3),
          vPQ[1] := vP20[1]-vQ10[1]:
          vPQ[2] := vP20[2]-vQ10[2]:
          vPQ[3] := vP20[3]-vQ10[3]:          
          #printf("vPQ: %f,%f,%f\n", vPQ[1],vPQ[2],vPQ[3]):
          
          #Calculate normal and tangental velocity magnitudes 
          velNormal := vPQ[1]*nP[1] 
                     + vPQ[2]*nP[2] 
                     + vPQ[3]*nP[3]:
          #printf("velNormal: %f\n", velNormal): 
                        
          velTangent[1] := vPQ[1]*(tP1[1] + tP2[1]):
          velTangent[2] := vPQ[2]*(tP1[2] + tP2[2]):
          velTangent[3] := vPQ[3]*(tP1[3] + tP2[3]):
          #printf("velTangent: %f,%f,%f \n", velTangent[1],velTangent[2],velTangent[3]):              
                        
          velTangentMag :=  (velTangent[1]*velTangent[1] 
                           + velTangent[2]*velTangent[2]
                           + velTangent[3]*velTangent[3])^0.5:                           
          #printf("velTangentMag: %f\n\n", velTangentMag): 
 
        ########################################################################
        # Compute the coefficient of friction using the cubic spline function
        ########################################################################
        mu := 0:
        delta:= 0:
        
        if (0.0e0 < velTangentMag) then            
            #In dynamic region
            if (velTangentMag > dyVel) then
                mu := mud:
                #printf("Dynamic Region\n");
            #Between dynamic and static            
            elif (velTangentMag > stVel and velTangentMag < dyVel) then
                delta := (velTangentMag - stVel)/(dyVel-stVel):
                mu := mus +(mud-mus)*delta*delta*(3.0-2.0*delta):
                #printf("Static to Dynamic Region\n");
            #In static region
            else
                delta := (velTangentMag + stVel)/(2.0*stVel):
                mu := -mus + (2*mus)*delta*delta*(3.0-2.0*delta):                
                #printf("Static Region\n");
            fi:
        fi:    
        #printf("mu(%f), mus(%f), mud(%f), vels(%f), veld(%f)\n\n", mu,mus,mud,stVel,dyVel):
        
        ########################################################################
        # Compute the contact force magnitudee using the 
        # Hunt-Crossely contact model
        ########################################################################
        forceNormal := k*(abs(gN)^p)*(1-c*velNormal):
        
        if(forceNormal < 0) then
            forceNormal:=0:
        fi:
        #printf("Normal Force: %f\n", forceNormal):
        
        ########################################################################
        # Compute the friction force magnitude using a Coulomb model
        ########################################################################
        forceTangentMag := -mu*forceNormal:
        #printf("Tangential Force: %f\n", forceTangentMag);

        ########################################################################
        # Compute the friction force direction using Y.Gonthier's numerically
        # stable method defined in Eqn. 20 of
        #
        # Gonthier,Y., McPhee,J., Lange,C., and Piedboeuf,J.C. (2004). 
        # A Regularized Contact Model with Asymmetric Damping and Dwell-Time
        # Dependent Friction. Multibody System Dynamics (11): 209-233.
        ########################################################################
        velEps := 1000.0 * epsRoot:
        
        if (velTangentMag > velEps) then
            forceTangentUV[1] := velTangent[1] / velTangentMag:
            forceTangentUV[2] := velTangent[2] / velTangentMag:
            forceTangentUV[3] := velTangent[3] / velTangentMag:
        else
            tmpf1 := (velTangent[1] / velEps):
            forceTangentUV[1] := tmpf1 * ((3.0/2.0)*abs(tmpf1) 
                                          -(1.0/2.0)*abs(tmpf1*tmpf1*tmpf1)):
            tmpf1 := (velTangent[2] / velEps):
            forceTangentUV[2] := tmpf1 * ((3.0/2.0)*abs(tmpf1)
                                          -(1.0/2.0)*abs(tmpf1*tmpf1*tmpf1)):
            tmpf1 := (velTangent[3] / velEps):
            forceTangentUV[3] := tmpf1 * ((3.0/2.0)*abs(tmpf1)
                                          -(1.0/2.0)*abs(tmpf1*tmpf1*tmpf1)):        
        fi:
                
        ########################################################################
        # Apply the forces to point P
        ########################################################################
        
        #Calculate the vector rQM which goes from the origin of the
        #disk frame to the point where the vector rPM intersects
        #the contact plane
        rQM[1] := rQO[1] - rMO[1]:
        rQM[2] := rQO[2] - rMO[2]:
        rQM[3] := rQO[3] - rMO[3]: 
        
        #Forces applied to point P on the disk
        aOut[1] := forceNormal*nP[1] + forceTangentMag*forceTangentUV[1]:
        aOut[2] := forceNormal*nP[2] + forceTangentMag*forceTangentUV[2]:
        aOut[3] := forceNormal*nP[3] + forceTangentMag*forceTangentUV[3]:
        #printf("Forces on disk: %f, %f, %f\n", aOut[1],aOut[2],aOut[3]):
        
        #Moments applied to the center of the disk
        #M = rPM x F
        aOut[4] := rPM[2]*aOut[3] - rPM[3]*aOut[2]:# a[2]*b[3]-a[3]*b[2] 
        aOut[5] :=-rPM[1]*aOut[3] + rPM[3]*aOut[1]:#-a[1]*b[3]+a[3]*b[1]
        aOut[6] := rPM[1]*aOut[2] - rPM[2]*aOut[1]:# a[1]*b[2]-a[2]*b[1]
        #printf("Moments on disk: %f, %f, %f\n", aOut[4],aOut[5],aOut[6]):
        
        #Forces applied to the plane (equal and opposite of the forces
        #applied to point P of the disk).
        aOut[7] := -aOut[1]:
        aOut[8] := -aOut[2]:
        aOut[9] := -aOut[3]:
        #printf("Forces on Ground: %f, %f, %f\n", aOut[7],aOut[8],aOut[9]):        
        
        #Moments applied to the origin of the plane. Forces are applied
        #at point P, which is actually inside the plane. Technically this
        #means that friction forces will create moments about the origin
        #of the plane which are not normal to the plane ... though the 
        #other alternatives are also unphysical in some way. Welcome
        #to compliant contact models.
        aOut[10] := rPMO[2]*aOut[9] - rPMO[3]*aOut[8]:# a[2]*b[3]-a[3]*b[2] 
        aOut[11] :=-rPMO[1]*aOut[9] + rPMO[3]*aOut[7]:#-a[1]*b[3]+a[3]*b[1]
        aOut[12] := rPMO[1]*aOut[8] - rPMO[2]*aOut[7]:# a[1]*b[2]-a[2]*b[1]
        #printf("Moments on Ground: %f, %f, %f\n", aOut[10],aOut[11],aOut[12]):
        
    fi:
    return 0.0:
        
    
end proc:
        
    
   
