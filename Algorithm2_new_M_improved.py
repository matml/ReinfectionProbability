#import matplotlib.pyplot as plt
import random, math
import numpy as np
from numpy import linalg
from numpy.linalg import inv
#import scipy
#from scipy import linalg
import time

			
def i(m,n,N):
	return 1.0/m

def ibar(m,n,N):
	return ((float)(m-1))/((float)(m))
	
def e(m,n,N):
	return 1.0/(N-m-n)

def ebar(m,n,N):
	return ((float)((N-m-n-1)))/((float)((N-m-n)))

def r(m,n,N):
	return ((float)(1.0))/((float)(n))
	
def rbar(m,n,N):
	return ((float)(n-1))/((float)(n))
	
def d(m,n,N):
	return ((float)(1.0))/((float)(n))
	
def dbar(m,n,N):
	return ((float)(n-1.0))/((float)(n))
	
def w(m,n,N):
	return ((float)(1.0))/((float)(N-m-n))

def wbar(m,n,N):
	return ((float)(N-m-n-1.0))/((float)(N-m-n))

def BuildAkkSI(Akk,k):
	m=k-1
	n=1
	thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
	if k>2:
		Akk[0,1]=bet*(k-1)*ibar(k-1,1,N)/thetamn
	for i in range(1,k-2):
		m=k-(i+1)
		n=i+1
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		Akk[i,i-1]=alpI*n/thetamn
		Akk[i,i+1]=bet*m*n*ibar(m,n,N)/thetamn
	i=k-2
	m=k-(i+1)
	n=i+1
	thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
	if k>2:
		Akk[i,i-1]=alpI*n/thetamn

def BuildAkkMinus1SI(AkkMinus1,k):
	for i in range(1,k-1):
		m=k-(i+1)
		n=i+1
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		AkkMinus1[i,i-1]=gam*n/thetamn

def BuildAkkPlus1SI(AkkPlus1,k):
	for i in range(1,k):
		m=k-i
		n=i
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		AkkPlus1[i-1,i-1]=alpR*(N-k)/thetamn
		AkkPlus1[i-1,i]=sig*bet*(N-k)*n/thetamn

def BuildbkSI(bk,k):
	for j in range(1,k):
		m=k-j
		n=j
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		bk[j-1,0]=bet*m*n*i(m,n,N)/thetamn

def BuildAkkR(Akk,k):
	m=k-1
	n=1
	thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
	if k>2:
		Akk[0,1]=bet*(k-1)/thetamn
	for i in range(1,k-1):
		m=k-(i+1)
		n=i+1
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		Akk[i,i-1]=alpI*n/thetamn
		Akk[i,i+1]=bet*m*n/thetamn
	i=k-1
	m=k-(i+1)
	n=i+1
	thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
	if k>2:
		Akk[i,i-1]=alpI*n/thetamn

def BuildAkkMinus1R(AkkMinus1,k):
	for i in range(1,k):
		m=k-(i+1)
		n=i+1
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		AkkMinus1[i,i-1]=gam*n/thetamn

def BuildAkkPlus1R(AkkPlus1,k):
	for i in range(1,k+1):
		m=k-i
		n=i
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		AkkPlus1[i-1,i-1]=alpR*(N-k)*wbar(m,n,N)/thetamn
		AkkPlus1[i-1,i]=sig*bet*(N-k)*n*ebar(m,n,N)/thetamn

def BuildbkRj(j,bRj,k,pSj,pIjPlus1):
	for i in range(1,k+1):
		m=k-i
		n=i
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		if j==M-1:
			prod=1
		else:
			prod=pIjPlus1[k-1][n-1,0]
		bRj[k-1][i-1,0]=(alpR*(N-k)*w(m,n,N)*pSj[k][n-1,0]+sig*bet*(N-k)*n*e(m,n,N)*prod)/thetamn

def BuildAkkI(Akk,k):
	m=k-1
	n=1
	thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
	if k>2:
		Akk[0,1]=bet*(k-1)/thetamn
	for i in range(1,k-1):
		m=k-(i+1)
		n=i+1
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		Akk[i,i-1]=alpI*n*dbar(m,n,N)/thetamn
		Akk[i,i+1]=bet*m*n/thetamn
	i=k-1
	m=k-(i+1)
	n=i+1
	thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
	if k>2:
		Akk[i,i-1]=alpI*n*dbar(m,n,N)/thetamn

def BuildAkkMinus1I(AkkMinus1,k):
	for i in range(1,k):
		m=k-(i+1)
		n=i+1
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		AkkMinus1[i,i-1]=gam*n*rbar(m,n,N)/thetamn

def BuildAkkPlus1I(AkkPlus1,k):
	for i in range(1,k+1):
		m=k-i
		n=i
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		AkkPlus1[i-1,i-1]=alpR*(N-k)/thetamn
		AkkPlus1[i-1,i]=sig*bet*(N-k)*n/thetamn

def BuildbkIj(j,bIj,k,pSj,pRj):
	for i in range(2,k+1):
		m=k-i
		n=i
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		bIj[k-1][i-1,0]=(n*gam*r(m,n,N)*pRj[k-2][n-2,0]+n*alpI*d(m,n,N)*pSj[k-1][n-2,0])/thetamn

def BuildAkkS(Akk,k):
	m=k-1
	n=1
	thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
	if k>2:
		Akk[0,1]=bet*(k-1)*ibar(k-1,1,N)/thetamn
	for i in range(1,k-2):
		m=k-(i+1)
		n=i+1
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		Akk[i,i-1]=alpI*n/thetamn
		Akk[i,i+1]=bet*m*n*ibar(m,n,N)/thetamn
	i=k-2
	m=k-(i+1)
	n=i+1
	thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
	if k>2:
		Akk[i,i-1]=alpI*n/thetamn

def BuildAkkMinus1S(AkkMinus1,k):
	for i in range(1,k-1):
		m=k-(i+1)
		n=i+1
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		AkkMinus1[i,i-1]=gam*n/thetamn

def BuildAkkPlus1S(AkkPlus1,k):
	for i in range(1,k):
		m=k-i
		n=i
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		AkkPlus1[i-1,i-1]=alpR*(N-k)/thetamn
		AkkPlus1[i-1,i]=sig*bet*(N-k)*n/thetamn

def BuildbkSj(j,bSj,k,pIjPlus1):
	for p in range(1,k):
		m=k-p
		n=p
		thetamn=bet*m*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		if j==M-1:
			prod=1
		else:
			prod=pIjPlus1[k-1][n,0]
		bSj[k-1][n-1,0]=bet*m*n*i(m,n,N)*prod/thetamn

def AlgorithmSj(j,invHS,PS,bS,pI,pS,AkkS,AkkMinus1S,AkkPlus1S):
	k=2
	BuildbkSj(j,bS[j],k,pI[j])
	PS[j][k-1]=bS[j][k-1]
	for k in range(3,N+1):
		BuildbkSj(j,bS[j],k,pI[j])
		PS[j][k-1]=AkkMinus1S[k-2]*invHS[k-2]*PS[j][k-2]+bS[j][k-1]
	k=N
	pS[j][k-1]=invHS[k-1]*PS[j][k-1]
	for k in reversed(range(2,N)):
		pS[j][k-1]=invHS[k-1]*(AkkPlus1S[k-1]*pS[j][k]+PS[j][k-1])

		
def AlgorithmRj(j,invHR,PR,bR,pS,pI,pR,AkkR,AkkMinus1R,AkkPlus1R):
	k=1
	BuildbkRj(j,bR[j-1],k,pS[j],pI[j])
	PR[j-1][k-1]=bR[j-1][k-1]
	for k in range(2,N):
		BuildbkRj(j,bR[j-1],k,pS[j],pI[j])
		PR[j-1][k-1]=AkkMinus1R[k-2]*invHR[k-2]*PR[j-1][k-2]+bR[j-1][k-1]
	k=N-1
	pR[j-1][k-1]=invHR[k-1]*PR[j-1][k-1]
	for k in reversed(range(1,N-1)):
		pR[j-1][k-1]=invHR[k-1]*(AkkPlus1R[k-1]*pR[j-1][k]+PR[j-1][k-1])
		
def AlgorithmIj(j,invHI,PI,bI,pS,pR,pI,AkkI,AkkMinus1I,AkkPlus1I):
	k=1
	BuildbkIj(j,bI[j-1],k,pS[j],pR[j-1])
	PI[j-1][k-1]=bI[j-1][k-1]
	for k in range(2,N+1):
		BuildbkIj(j,bI[j-1],k,pS[j],pR[j-1])
		PI[j-1][k-1]=AkkMinus1I[k-2]*invHI[k-2]*PI[j-1][k-2]+bI[j-1][k-1]
	k=N
	pI[j-1][k-1]=invHI[k-1]*PI[j-1][k-1]
	for k in reversed(range(1,N)):
		pI[j-1][k-1]=invHI[k-1]*(AkkPlus1I[k-1]*pI[j-1][k]+PI[j-1][k-1])	

def probabilities(gam,bet,alpI,alpR,sig,N,M,i0,s0,r0):
	

	HS=[np.zeros((k-1,k-1)) for k in range(1,N+1)]
	invHS=[np.zeros((k-1,k-1)) for k in range(1,N+1)]
	AkkS=[np.zeros((k-1,k-1)) for k in range(1,N+1)]
	AkkMinus1S=[np.zeros((k-1,k-2)) for k in range(2,N+1)]
	AkkPlus1S=[np.zeros((k-2,k-1)) for k in range(2,N+1)]

	k=2
	BuildAkkSI(AkkS[k-1],k)
	HS[k-1]=np.eye(k-1)-AkkS[k-1]
	invHS[k-1]=inv(HS[k-1])

	for k in range(3,N+1):
		BuildAkkSI(AkkS[k-1],k)
		BuildAkkMinus1SI(AkkMinus1S[k-2],k)
		BuildAkkPlus1SI(AkkPlus1S[k-2],k-1)
		HS[k-1]=np.asmatrix(np.eye(k-1))-AkkS[k-1]-AkkMinus1S[k-2]*invHS[k-2]*AkkPlus1S[k-2]
		invHS[k-1]=inv(HS[k-1])

	HR=[np.zeros((k,k)) for k in range(1,N+1)]
	invHR=[np.zeros((k,k)) for k in range(1,N+1)]
	AkkR=[np.zeros((k,k)) for k in range(1,N)]
	AkkMinus1R=[np.zeros((k,k-1)) for k in range(2,N)]
	AkkPlus1R=[np.zeros((k-1,k)) for k in range(2,N)]

	k=1
	BuildAkkR(AkkR[k-1],k)
	HR[k-1]=np.asmatrix(np.eye(k))-AkkR[k-1]
	invHR[k-1]=inv(HR[k-1])

	for k in range(2,N):
		BuildAkkR(AkkR[k-1],k)
		BuildAkkMinus1R(AkkMinus1R[k-2],k)
		BuildAkkPlus1R(AkkPlus1R[k-2],k-1)
		HR[k-1]=np.asmatrix(np.eye(k))-AkkR[k-1]-AkkMinus1R[k-2]*invHR[k-2]*AkkPlus1R[k-2]
		invHR[k-1]=inv(HR[k-1])

	HI=[np.zeros((k,k)) for k in range(1,N+1)]
	invHI=[np.zeros((k,k)) for k in range(1,N+1)]
	AkkI=[np.zeros((k,k)) for k in range(1,N+1)]
	AkkMinus1I=[np.zeros((k,k-1)) for k in range(2,N+1)]
	AkkPlus1I=[np.zeros((k-1,k)) for k in range(2,N+1)]

	k=1
	BuildAkkI(AkkI[k-1],k)
	HI[k-1]=np.asmatrix(np.eye(k))-AkkI[k-1]
	invHI[k-1]=inv(HI[k-1])

	for k in range(2,N+1):
		BuildAkkI(AkkI[k-1],k)
		BuildAkkMinus1I(AkkMinus1I[k-2],k)
		BuildAkkPlus1I(AkkPlus1I[k-2],k-1)
		HI[k-1]=np.eye(k)-AkkI[k-1]-AkkMinus1I[k-2]*invHI[k-2]*AkkPlus1I[k-2]
		invHI[k-1]=inv(HI[k-1])

	pS=[[np.zeros((k-1,1)) for k in range(1,N+1)] for j in range(M)]
	pI=[[np.zeros((k,1)) for k in range(1,N+1)] for j in range(M)]
	pR=[[np.zeros((k,1)) for k in range(1,N+1)] for j in range(M)]
	PS=[[np.zeros((k-1,1)) for k in range(1,N+1)] for j in range(M)]
	PI=[[np.zeros((k,1)) for k in range(1,N+1)] for j in range(M)]
	PR=[[np.zeros((k,1)) for k in range(1,N+1)] for j in range(M)]
	bS=[[np.zeros((k-1,1)) for k in range(1,N+1)] for j in range(M)]
	bI=[[np.zeros((k,1)) for k in range(1,N+1)] for j in range(M)]
	bR=[[np.zeros((k,1)) for k in range(1,N+1)] for j in range(M)]

	j=M-1
	AlgorithmSj(j,invHS,PS,bS,pI,pS,AkkS,AkkMinus1S,AkkPlus1S)
	for j in reversed(range(M-1)):
		AlgorithmRj(j+1,invHR,PR,bR,pS,pI,pR,AkkR,AkkMinus1R,AkkPlus1R)
		AlgorithmIj(j+1,invHI,PI,bI,pS,pR,pI,AkkI,AkkMinus1I,AkkPlus1I)
		AlgorithmSj(j,invHS,PS,bS,pI,pS,AkkS,AkkMinus1S,AkkPlus1S)

	return pS#[0][s0+i0-1][i0-1,0]
