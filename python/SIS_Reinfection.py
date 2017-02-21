#import matplotlib.pyplot as plt
import random, math
import numpy as np
from numpy import linalg
			
#def i(m,N):
#	return 1.0/m

#def ibar(m,N):
#	return ((float)(m-1))/((float)(m))
	
#def s(m,N):
#	return 1.0/(N-m)

#def sbar(m,N):
#	return ((float)((N-m-1)))/((float)((N-m)))

def Sj(gam,gamA,bet,betAc,betcA,j,hS,fS,pI,pS):
	gSj=np.asmatrix(np.zeros((N,1)))
	m=N-1
	thetam=bet*(m-1)*(N-m)+betcA*(N-m)+gam*(N-m)
	if j==M-1:
		prod=1
	else:
		prod=pI[m-1,j+1]
	gSj[m,0]=betcA*(N-m)*prod/thetam
	for m in reversed(range(1,N-1)):
		thetam=bet*(m-1)*(N-m)+betcA*(N-m)+gam*(N-m)
		if j==M-1:
			prod=1
		else:
			prod=pI[m-1,j+1]
		gSj[m,0]=betcA*(N-m)*prod/thetam+gam*(N-m)/thetam*gSj[m+1,0]/hS[m+1,0]
	m=1
	pS[m,j]=gSj[m,0]/hS[m,0]
	for m in range(2,N):
		pS[m,j]=(fS[m,0]*pS[m-1,j]+gSj[m,0])/hS[m,0]

def Ij(gam,gamA,bet,betAc,betcA,j,hI,fI,pS,pI):
	gIj=np.asmatrix(np.zeros((N,1)))
	m=N-1
	thetam=bet*m*(N-m-1)+betAc*m+gam*(N-m-1)+gamA
	gIj[m,0]=0
	for m in reversed(range(N-1)):
		thetam=bet*m*(N-m-1)+betAc*m+gam*(N-m-1)+gamA
		gIj[m,0]=gam*(N-m-1)/thetam*gIj[m+1,0]/hI[m+1,0]+gamA/thetam*pS[m+1,j]
	m=0
	pI[m,j]=gIj[m,0]/hI[m,0]
	for m in range(1,N):
		pI[m,j]=(fI[m,0]*pI[m-1,j]+gIj[m,0])/hI[m,0]	


def p_reinfection_SIS(gam,gamA,bet,betAc,betcA,N,M,i0,s0):
	m=N-1
	hS=np.asmatrix(np.zeros((N,1)))
	fS=np.asmatrix(np.zeros((N,1)))
	thetam=bet*(m-1)*(N-m)+betcA*(N-m)+gam*(N-m)
	hS[m,0]=1.0
	fS[m,0]=bet*(m-1)*(N-m)/thetam
	for m in reversed(range(1,N-1)):
		thetam=bet*(m-1)*(N-m)+betcA*(N-m)+gam*(N-m)
		hS[m,0]=1-gam*(N-m)/thetam*fS[m+1,0]/hS[m+1,0]
		fS[m,0]=bet*(m-1)*(N-m)/thetam

	hI=np.asmatrix(np.zeros((N,1)))
	fI=np.asmatrix(np.zeros((N,1)))
	m=N-1
	thetam=bet*m*(N-m-1)+betAc*m+gam*(N-m-1)+gamA
	hI[m,0]=1.0
	fI[m,0]=betAc*m/thetam
	for m in reversed(range(N-1)):
		thetam=bet*m*(N-m-1)+betAc*m+gam*(N-m-1)+gamA
		hI[m,0]=1-gam*(N-m-1)/thetam*fI[m+1,0]/hI[m+1,0]
		fI[m,0]=(bet*m*(N-m-1)+betAc*m)/thetam

	pS=np.asmatrix(np.zeros((N,M)))
	pI=np.asmatrix(np.zeros((N,M)))
	j=M-1
	Sj(gam,gamA,bet,betAc,betcA,j,hS,fS,np.asmatrix(np.ones((N,1))),pS)
	for j in reversed(range(M-1)):
		Ij(gam,gamA,bet,betAc,betcA,j+1,hI,fI,pS,pI)
		Sj(gam,gamA,bet,betAc,betcA,j,hS,fS,pI,pS)
		
	return (N-1.0)/N*pS[s0,0]+1.0/N*pI[s0,0]

#gam=1.0
#bet=0.02
#betAc=0.02
#betcA=0.02
#gamA=1.0


#N=100
#M=2

#s0=N-1
#i0=1

#probability=p_reinfection_SIS(gam,gamA,bet,betAc,betcA,N,M,i0,s0)
#print(probability)
