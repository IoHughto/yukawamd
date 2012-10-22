#!/usr/bin/python

from random import*
from sys import*


if len(argv) != 5:
	exit("Syntax Error: <number of ions> <percentage of impurity> <charge of impurity> <mass>")

PART=int(argv[1])
frac=float(argv[2])
charge=float(argv[3])
mass=float(argv[4])

repnum=int(round(frac*PART))

restart=0
while(restart==0):
	reparr=[]

	for i in range(repnum):
		reparr.append(randint(1,PART))

	reparr.sort()

	counter=0

	for i in range(repnum-1):
		if(reparr[i]==reparr[i-1]):
			counter+=1

	if(counter==0):
		restart=1	

print(PART)

for i in range(1,PART+1):
	counter=0
	for j in range(repnum):
		if(i-1==reparr[j]):
			counter+=1
	if(counter==0):
		s=str(i)+' 34.00 80.00'
	else:
		s=str(i)+' '+str(charge)+' '+str(mass)
	print s
