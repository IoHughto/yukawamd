#!/usr/bin/python

import sys
import math 
import os

PART=int(raw_input('Number of particles: '))
spec=int(raw_input('Number of species: '))

total=0

CH=[]
MA=[]

for i in range(spec):
	print '\nSpecies {0}:'.format(i+1)
	charge=float(raw_input('Charge: '))
	mass=float(raw_input('Mass: '))
	if(spec>1):
		n=int(raw_input('Number: '))
	else:
		n=PART
	total+=n
	for j in range(n):
		CH.insert(len(CH),charge)
		MA.insert(len(MA),mass)
	
if total != PART:
	print 'Particle number mismatch.'
	sys.exit()

fname=raw_input('\nOutput file: ')
if(os.path.isfile(fname)):
	print '\nFile exists!'
	ans=raw_input('Overwrite? ')
	if(ans.lower() != 'y'):
		print '\nProgram exiting'
		sys.exit()
f=open(fname,'w')

f.write(str(PART)+'\n')

for i in range(PART):
	f.write('{0} {1} {2}\n'.format(i+1,CH[i],MA[i]))

f.close()




			
