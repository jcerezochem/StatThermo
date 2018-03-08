#!/usr/bin/env python

import numpy

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
    temp_au2K = numpy.float64(3.1577464e5)
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

        return 2.0/float(Ndf)*Ekin*const.temp_au2K



if __name__ == "__main__":
    import sys

    g96file  = sys.argv[1]
    massfile = sys.argv[2]

    v    = numpy.array(parse_g96(g96file,'VELOCITY'))
    mass = read_massfile(massfile)
    Nat = len(v)/3

    print '%8.2f'%compute_temp(v,mass,3*Nat-3)

