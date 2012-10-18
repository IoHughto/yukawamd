#!/usr/bin/python

import sys
import math 

z=zsq=m=0.0
try:
	rho=float(sys.argv[2])
except:
	s="Syntax Error: "+sys.argv[0]+" <ion.dat file> <density> (<Gamma/temperature>)"
        print s
        sys.exit()

try:
	f=open(sys.argv[1],"r")
except:
	s="Syntax Error: "+sys.argv[0]+" <ion.dat file> <density> (<Gamma/temperature>)"
        print s
        sys.exit()

lines=f.readlines()
for a in range(1,len(lines)):
	z+=float(lines[a].split()[1])
	m+=float(lines[a].split()[2])
	zsq+=float(lines[a].split()[1])**(5.0/3.0)

z/=len(lines)-1
m/=len(lines)-1
zsq/=len(lines)-1

ee=197.326938/137.03599
a=(4*math.pi*rho/3)**(-1.0/3.0)

gtbase=zsq*z**(1.0/3.0)*ee/a
omega=(4*math.pi*ee*z**2.0*rho/(m*931))**0.5
tp=1/omega
waa=omega*a**2.0
lamb=(2.0*(math.pi*137.035999)**(-0.5)*(3*math.pi**2.0*z*rho)**(1.0/3.0))
lamb=1/lamb
halfl=((len(lines)-1)/rho)**(1.0/3.0)/2.0

try:
	temp=float(sys.argv[3])
except IndexError:
	print 'gtbase = %f' % (gtbase)
else:
        print 'Gamma/T = %f (MeV)' % (gtbase/temp)
print 'omega = %f c/fm\n1/omega = %f fm/c\nomega*a^2 = %f fm c\nlambda = %f fm\nAverage Z = %f\nL/2 = %f fm\n' % (omega,tp,waa,lamb,z,halfl)
