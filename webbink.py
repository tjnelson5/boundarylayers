# function: webbink.py
'''
This is a translation of Koji;'s Fortran code.

Given a white dwarf mass M (Msun), calculates its radius using Webbink's 
formula (see Pringle & Webbink 1975, MN 172, 493).  Returns radius in cm

Webbink's formula is not the latest, for example see Koester 1978, A&A 64, 289.
Moreover, mass-radius relationship depends on the composition and core 
temperature of the white dwarf. However, Webbink's formula is easy to use and is
probably good enough for many purposes.

Original program by Koji, 1988 Feb
function version by Koji, 1995 Aug
IDL version by Tommy, 2011 Feb
Python version by Tommy, 2012 Feb
'''

import sys,math

def radius(mass):
    if (mass <= 0.0):
        print "White dwarf mass must be positive and less than 1.44 Msol"
	return
    elif (mass >= 1.44):
        print "White dwarf mass must be positive and less than 1.44 Msol"
	return
    else:     
        x=(1.44/mass)-1.0
        temp1=0.3767-(0.00605*math.log10(x))
        radius=7.7e8*(math.pow(x,temp1))
        return radius
    
def logg(mass):
    if (mass <= 0.0):
        print "White dwarf mass must be positive and less than 1.44 Msol"
	return
    elif (mass >= 1.44):
        print "White dwarf mass must be positive and less than 1.44 Msol"
	return
    else:     
        x=(1.44/mass)-1.0
        temp1=0.3767-(0.00605*math.log10(x))
        radius=7.7e8*(math.pow(x,temp1))
	gravConst=6.67e-8
        solarMass=2.0e33
        g=gravConst*mass*solarMass/(math.pow(radius,2.0))
        logg=math.log10(g)
        return logg

