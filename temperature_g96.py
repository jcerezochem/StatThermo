#!/usr/bin/env python

import numpy
import re

#Constants
class const:
    # Conversion factors
    BOHRtoANGS= numpy.float64(5.2917720859e-1)  
    UMAtoKG   = numpy.float64(1.66053873e-27)   
    UMAtoAU   = numpy.float64(1.82288839e3)    
    AMUtoAU   = numpy.float64(1.82288839e3)  
    AMUtoKG   = numpy.float64(1.66053873e-27)  
    AUtoKG    = numpy.float64(9.10938291e-31)  
    BOHRtoM   = numpy.float64(5.291772083e-11)  
    BOHRtoNM  = numpy.float64(5.291772083e-2)            
    ANGStoM   = numpy.float64(1.e-10)   
    HARTtoJ   = numpy.float64(4.3597482e-18)    
    HtoKCALM  = numpy.float64(627.5095e0)       
    CALtoJ    = numpy.float64(4.184)         
    HtoeV     = numpy.float64(27.2114)  
    autown    = numpy.float64(2.1947463068e5)
    AUtoKelv  = numpy.float64(3.1577464e5)
    AUtoSEC   = numpy.float64(2.418884326505e-17)
    # Universal constants
    PI      = numpy.float64(3.14159265358979323846e0)
    NAv     = numpy.float64(6.02214129e23)      
    clight  = numpy.float64(2.99792458e8)       
    SL      = numpy.float64(2.99792458e8)       
    plank   = numpy.float64(6.62606957e-34)    
    plankbar= numpy.float64(1.054571726e-34)     
    kboltz  = numpy.float64(1.3806488e-23)      
    boltz   = numpy.float64(1.3806488e-23)    
    atmass  = numpy.float64(1.660538921e-27)     
    cvelau  = numpy.float64(137.0369)
    ZERO    = numpy.float64(1.0e-8)

def parse_g96(g96name,section):
    
    # If we ask for 'atnames' or 'resnames', read atomnames and resnames from POSITION
    if section == 'atnames' or section == 'resnames':
        section_='POSITION'
    else:
        section_=section
        
    # Read data from file
    with open(g96name) as f:
        
        for line in f:
            if section_ in line:
                break
            
        if section_ == "POSITION" or section_ == "VELOCITY":
            aname= []
            rname= []
            xyz  = []
            line = f.next()
            while line.strip() != "END":
                data = line.split()
                rname.append(data[1])
                aname.append(data[2])
                for i in data[4:7]:
                    xyz.append(float(i)*10.0)
                #xyz.append([float(i)*10.0 for i in data[4:7]])
                line = f.next()
            
    if section == "POSITION" or section == "VELOCITY":
        return xyz
    elif section == "atnames":
        return aname
    elif section == "resnames":
        return rname
        
def read_massfile(massfile):
    mass=[]
    with open(massfile) as f:
        
        for line in f:
            if len(line.strip()) == 0:
                continue
            mass.append(float(line.split()[0]))
            
    return mass

def compute_temp(v,mass,Ndf):
	"""
        This functions expects:
        v    in Anng/ps
        mass in AMU
        """

        mass = numpy.array(mass)
        v    = numpy.array(v)

        #Unit conversion (to AU)
        mass *= const.UMAtoAU
        v    /= const.BOHRtoANGS   # [L]
        v    *= 1e12               # [T^-1]: ps-1 --> s-1
        v    *= const.AUtoSEC      # [T^-1]: s-1  --> au

        Ekin=0.0
        for m,x,y,z in zip(mass,v[0::3],v[1::3],v[2::3]):
            Ekin += m*(x**2+y**2+z**2)
        Ekin /= 2.0

        return 2.0/float(Ndf)*Ekin*const.AUtoKelv


def compute_Ekin(v,mass):
        """
        This functions expects:
        v    in Anng/ps
        mass in AMU
        """

        mass = numpy.array(mass)
        v    = numpy.array(v)

        #Unit conversion (to AU)
        mass *= const.UMAtoAU
        v    /= const.BOHRtoANGS   # [L]
        v    *= 1e12               # [T^-1]: ps-1 --> s-1
        v    *= const.AUtoSEC      # [T^-1]: s-1  --> au

        Ekin=0.0
        for m,x,y,z in zip(mass,v[0::3],v[1::3],v[2::3]):
            Ekin += m*(x**2+y**2+z**2)
        Ekin /= 2.0

        return Ekin*const.HtoKCALM*const.CALtoJ



def Cart2Eckart(xyz,mass):
    
    """
    Function to compute the D matrix, used to transform from Cartesian to Eckart frame
    
      q = D t
    where t are the Eckart coordinates and q are mwc
    """

    # Compute COM
    rcom=numpy.zeros(3)
    rcom[0] = (xyz[0::3]*mass).sum()/mass.sum()
    rcom[1] = (xyz[1::3]*mass).sum()/mass.sum()
    rcom[2] = (xyz[2::3]*mass).sum()/mass.sum()
    # Substract com
    xyz[0::3] -= rcom[0]
    xyz[1::3] -= rcom[1]
    xyz[2::3] -= rcom[2]   
    
    # Get moment of inertia
    MI = numpy.zeros((3,3))
    MI[0,0] = (mass*(xyz[1::3]**2+xyz[2::3]**2)).sum()
    MI[1,1] = (mass*(xyz[0::3]**2+xyz[2::3]**2)).sum()
    MI[2,2] = (mass*(xyz[0::3]**2+xyz[1::3]**2)).sum()
    MI[0,1] = -(mass*xyz[0::3]*xyz[1::3]).sum()
    MI[1,0] = MI[0,1]
    MI[0,2] = -(mass*xyz[0::3]*xyz[2::3]).sum()
    MI[2,0] = MI[0,2]
    MI[1,2] = -(mass*xyz[1::3]*xyz[2::3]).sum()
    MI[2,1] = MI[1,2]
    # Get principal axes
    Idiag,Xrot = numpy.linalg.eigh(MI)
    print '-------------------------------'
    print ' DIAGONAL INTERTIA MOMENTS'
    print '-------------------------------'
    print 'Ix = {:12.5f}'.format(Idiag[0])
    print 'Iy = {:12.5f}'.format(Idiag[1])
    print 'Iz = {:12.5f}'.format(Idiag[2])
    print ''
    # Rotate to principal axes
    # ...
    
    D = numpy.zeros((3*Nat,3*Nat+6))
    # Translations
    D[0::3,0] = numpy.sqrt(mass)
    D[1::3,1] = numpy.sqrt(mass)
    D[2::3,2] = numpy.sqrt(mass)
    # Rotation (note: Xrot is used to account for the transpose, this is different wrt to Fortran code)
    for i in range(Nat):
        R = Xrot.transpose().dot(xyz[3*i:3*i+3])
        D[3*i  ,3] = (R[1]*Xrot[0,2] - R[2]*Xrot[0,1])*numpy.sqrt(mass[i])
        D[3*i+1,3] = (R[1]*Xrot[1,2] - R[2]*Xrot[1,1])*numpy.sqrt(mass[i]) 
        D[3*i+2,3] = (R[1]*Xrot[2,2] - R[2]*Xrot[2,1])*numpy.sqrt(mass[i])
        D[3*i  ,4] = (R[2]*Xrot[0,0] - R[0]*Xrot[0,2])*numpy.sqrt(mass[i])
        D[3*i+1,4] = (R[2]*Xrot[1,0] - R[0]*Xrot[1,2])*numpy.sqrt(mass[i]) 
        D[3*i+2,4] = (R[2]*Xrot[2,0] - R[0]*Xrot[2,2])*numpy.sqrt(mass[i])
        D[3*i  ,5] = (R[0]*Xrot[0,1] - R[1]*Xrot[0,0])*numpy.sqrt(mass[i])
        D[3*i+1,5] = (R[0]*Xrot[1,1] - R[1]*Xrot[1,0])*numpy.sqrt(mass[i]) 
        D[3*i+2,5] = (R[0]*Xrot[2,1] - R[1]*Xrot[2,0])*numpy.sqrt(mass[i])
    
    # Normalize
    idel=[]
    for i in range(6):
        norm = numpy.linalg.norm(D[:,i]) #numpy.sqrt(D[:,i].dot(D[:,i]))
        if norm > const.ZERO:
            D[:,i] /= norm
        else:
            idel.append(i)
    
    # Set number of rotations and trasnlations
    ntr = 6  
    if len(idel)==1:
        print 'Linear molecule detected'
        D = numpy.delete(D,idel[0],1)
        ntr -= 1
        is_linear=True
    elif len(idel)>1:
        raise Exception('Too many zero rotational vectors')
    else:
        is_linear=False
        
    # Now get internal coordinates by Grand-Schmidt orthogonalization
    D[:,ntr:] = numpy.identity(3*Nat)
    idel=[]
    for i in range(ntr,3*Nat+ntr):
        v = D[:,i]
        for j in range(i):
            projection = numpy.inner(v,D[:,j])
            D[:,i] -= projection*D[:,j]
    
        norm = numpy.linalg.norm(D[:,i])
        #print '%3i %23.15g'%(i+1, norm**2)
        if norm**2 > const.ZERO:
            D[:,i] /= norm
        else:
            idel.append(i)
            # Set to zero to avoid numerical errors
            #D[:,i] = 0.
            
    for i in idel[::-1]:
        D = numpy.delete(D,i,1)
            
            
    if D.shape[0] != D.shape[1]:
        raise Exception('Wrong number of zero roots during GS orthogonalization: '+str(len(idel)))
        
        
    return D, is_linear


if __name__ == '__main__':
    import sys

    g96file0 = sys.argv[1]
    massfile = sys.argv[2]

    
    # Get data for the original structure
    xyz      = numpy.array(parse_g96(g96file0,'POSITION'))
    v        = numpy.array(parse_g96(g96file0,'VELOCITY'))
    atnames  = parse_g96(g96file0,'atnames')
    Nat      = len(atnames)
    mass     = numpy.array(read_massfile(massfile))
    # And compute its D matrix
    D, is_linear = Cart2Eckart(xyz,mass)
    if is_linear:
        ntr=5
    else:
        ntr=6
 
    Tall  = compute_temp(v,mass,3*Nat)
    Ekall = compute_Ekin(v,mass)

    # Transform vel to mwc
    v_mwc = numpy.zeros(3*Nat)
    v_mwc[0::3] = v[0::3] * numpy.sqrt(mass)
    v_mwc[1::3] = v[1::3] * numpy.sqrt(mass)
    v_mwc[2::3] = v[2::3] * numpy.sqrt(mass)

    # Project to Tr+Rot
    vtr  = D[:,:ntr].transpose().dot(v_mwc)
    vint = D[:,ntr:].transpose().dot(v_mwc)
    print ''
    print 'Translation velocities'
    print '     x         y        z   '
    print '{:9.3f}{:9.3f}{:9.3f}'.format(vtr[0],vtr[1],vtr[2])
    print ''
    print 'Rotation velocities'
    if is_linear:
        print '     x         y   '
        print '{:9.3f}{:9.3f}'.format(vtr[3],vtr[4])
    else:
        print '     x         y        z   '
        print '{:9.3f}{:9.3f}{:9.3f}'.format(vtr[3],vtr[4],vtr[5])
    print ''


    # Compute Temperature of internal/external DOF 
    v_mwc = D[:,:ntr].dot(vtr)
    v[0::3] = v_mwc[0::3] / numpy.sqrt(mass)
    v[1::3] = v_mwc[1::3] / numpy.sqrt(mass)
    v[2::3] = v_mwc[2::3] / numpy.sqrt(mass)
    Ttr  = compute_temp(v,mass,ntr)
    Ektr = compute_Ekin(v,mass)
    v_mwc = D[:,ntr:].dot(vint)
    v[0::3] = v_mwc[0::3] / numpy.sqrt(mass)
    v[1::3] = v_mwc[1::3] / numpy.sqrt(mass)
    v[2::3] = v_mwc[2::3] / numpy.sqrt(mass)
    Tint  = compute_temp(v,mass,3*Nat-ntr)
    Ekint = compute_Ekin(v,mass)
  
    print '-------------------------------------------------------'
    print ' ENERGIES AND TEMPERATURES IN INTIAL CONFORMATION'
    print '-------------------------------------------------------'
    print 'Tr+Rot  : Ekin = {:12.3f}; T = {:9.3f}'.format(Ektr,Ttr)
    print 'Internal: Ekin = {:12.3f}; T = {:9.3f}'.format(Ekint,Tint)
    print 'Sum     : Ekin = {:12.3f}             '.format(Ekint+Ektr)
    print 'All DoF : Ekin = {:12.3f}; T = {:9.3f}'.format(Ekall,Tall)
    print ''
  

