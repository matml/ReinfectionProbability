#import matplotlib.pyplot as plt
import random, math
import numpy as np
# from numpy import linalg
# from numpy.linalg import inv
import scipy
from scipy import linalg
from scipy.linalg import inv
import time

			
def BuildAkkSI(Akk,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	m=k-1
	n=1
	thetamnSI=bet*(m-1)*n+betcA*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
	if k>2:
		Akk[0,1]=bet*(k-2)/thetamnSI
	for i in range(1,k-2):
		m=k-(i+1)
		n=i+1
		thetamnSI=bet*(m-1)*n+betcA*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		Akk[i,i-1]=alpI*n/thetamnSI
		Akk[i,i+1]=bet*(m-1)*n/thetamnSI
	i=k-2
	m=k-(i+1)
	n=i+1
	thetamnSI=bet*(m-1)*n+betcA*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
	if k>2:
		Akk[i,i-1]=alpI*n/thetamnSI


def BuildAkkMinus1SI(AkkMinus1,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	for i in range(1,k-1):
		m=k-(i+1)
		n=i+1
		thetamnSI=bet*(m-1)*n+betcA*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		AkkMinus1[i,i-1]=gam*n/thetamnSI

def BuildAkkPlus1SI(AkkPlus1,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	for i in range(1,k):
		m=k-i
		n=i
		thetamnSI=bet*(m-1)*n+betcA*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		AkkPlus1[i-1,i-1]=alpR*(N-k)/thetamnSI
		AkkPlus1[i-1,i]=sig*bet*(N-k)*n/thetamnSI

def BuildbkSI(bk,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	for j in range(1,k):
		m=k-j
		n=j
		thetamnSI=bet*(m-1)*n+betcA*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		bk[j-1,0]=betcA*n/thetamnSI

def BuildAkkR(Akk,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	m=k-1
	n=1
	thetamnR=bet*m*n+gam*n+alpI*n+alpR*(N-m-n-1)+alpRA+sig*bet*n*(N-m-n-1)+sigA*betcA*n
	if k>2:
		Akk[0,1]=bet*(k-1)/thetamnR
	for i in range(1,k-1):
		m=k-(i+1)
		n=i+1
		thetamnR=bet*m*n+gam*n+alpI*n+alpR*(N-m-n-1)+alpRA+sig*bet*n*(N-m-n-1)+sigA*betcA*n
		Akk[i,i-1]=alpI*n/thetamnR
		Akk[i,i+1]=bet*m*n/thetamnR
	i=k-1
	m=k-(i+1)
	n=i+1
	thetamnR=bet*m*n+gam*n+alpI*n+alpR*(N-m-n-1)+alpRA+sig*bet*n*(N-m-n-1)+sigA*betcA*n
	if k>2:
		Akk[i,i-1]=alpI*n/thetamnR

def BuildAkkMinus1R(AkkMinus1,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	for i in range(1,k):
		m=k-(i+1)
		n=i+1
		thetamnR=bet*m*n+gam*n+alpI*n+alpR*(N-m-n-1)+alpRA+sig*bet*n*(N-m-n-1)+sigA*betcA*n
		AkkMinus1[i,i-1]=gam*n/thetamnR

def BuildAkkPlus1R(AkkPlus1,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	for i in range(1,k+1):
		m=k-i
		n=i
		thetamnR=bet*m*n+gam*n+alpI*n+alpR*(N-m-n-1)+alpRA+sig*bet*n*(N-m-n-1)+sigA*betcA*n
		AkkPlus1[i-1,i-1]=alpR*(N-k-1)/thetamnR
		AkkPlus1[i-1,i]=sig*bet*(N-k-1)*n/thetamnR

def BuildbkRj(j,bRj,k,pSj,pIjPlus1,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M):
	for i in range(1,k+1):
		m=k-i
		n=i
		thetamnR=bet*m*n+gam*n+alpI*n+alpR*(N-m-n-1)+alpRA+sig*bet*n*(N-m-n-1)+sigA*betcA*n
		if j==M-1:
			prod=1
		else:
			prod=pIjPlus1[k-1][n-1,0]
		bRj[k-1][i-1,0]=(alpRA*pSj[k][n-1,0]+sigA*betcA*n*prod)/thetamnR

def BuildAkkI(Akk,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	m=k-1
	n=1
	thetamnI=bet*m*(n-1)+betAc*m+gam*(n-1)+gamA+alpI*(n-1)+alpIA+alpR*(N-m-n)+sig*bet*(n-1)*(N-m-n)+sig*betAc*(N-m-n)
	if k>2:
		Akk[0,1]=(bet*m*(n-1)+betAc*m)/thetamnI
	for i in range(1,k-1):
		m=k-(i+1)
		n=i+1
		thetamnI=bet*m*(n-1)+betAc*m+gam*(n-1)+gamA+alpI*(n-1)+alpIA+alpR*(N-m-n)+sig*bet*(n-1)*(N-m-n)+sig*betAc*(N-m-n)
		Akk[i,i-1]=alpI*(n-1)/thetamnI
		Akk[i,i+1]=(bet*m*(n-1)+betAc*m)/thetamnI
	i=k-1
	m=k-(i+1)
	n=i+1
	thetamnI=bet*m*(n-1)+betAc*m+gam*(n-1)+gamA+alpI*(n-1)+alpIA+alpR*(N-m-n)+sig*bet*(n-1)*(N-m-n)+sig*betAc*(N-m-n)
	if k>2:
		Akk[i,i-1]=alpI*(n-1)/thetamnI

def BuildAkkMinus1I(AkkMinus1,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	for i in range(1,k):
		m=k-(i+1)
		n=i+1
		thetamnI=bet*m*(n-1)+betAc*m+gam*(n-1)+gamA+alpI*(n-1)+alpIA+alpR*(N-m-n)+sig*bet*(n-1)*(N-m-n)+sig*betAc*(N-m-n)
		AkkMinus1[i,i-1]=gam*(n-1)/thetamnI

def BuildAkkPlus1I(AkkPlus1,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	for i in range(1,k+1):
		m=k-i
		n=i
		thetamnI=bet*m*(n-1)+betAc*m+gam*(n-1)+gamA+alpI*(n-1)+alpIA+alpR*(N-m-n)+sig*bet*(n-1)*(N-m-n)+sig*betAc*(N-m-n)
		AkkPlus1[i-1,i-1]=alpR*(N-k)/thetamnI
		AkkPlus1[i-1,i]=(sig*bet*(N-k)*(n-1)+sig*betAc*(N-k))/thetamnI

def BuildbkIj(j,bIj,k,pSj,pRj,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	for i in range(2,k+1):
		m=k-i
		n=i
		thetamnI=bet*m*(n-1)+betAc*m+gam*(n-1)+gamA+alpI*(n-1)+alpIA+alpR*(N-m-n)+sig*bet*(n-1)*(N-m-n)+sig*betAc*(N-m-n)
		bIj[k-1][i-1,0]=(gamA*pRj[k-2][n-2,0]+alpIA*pSj[k-1][n-2,0])/thetamnI

def BuildAkkS(Akk,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	m=k-1
	n=1
	thetamnS=bet*(m-1)*n+betcA*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
	if k>2:
		Akk[0,1]=bet*(k-2)/thetamnS
	for i in range(1,k-2):
		m=k-(i+1)
		n=i+1
		thetamnS=bet*(m-1)*n+betcA*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		Akk[i,i-1]=alpI*n/thetamnS
		Akk[i,i+1]=bet*(m-1)*n/thetamnS
	i=k-2
	m=k-(i+1)
	n=i+1
	thetamnS=bet*(m-1)*n+betcA*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
	if k>2:
		Akk[i,i-1]=alpI*n/thetamnS

def BuildAkkMinus1S(AkkMinus1,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	for i in range(1,k-1):
		m=k-(i+1)
		n=i+1
		thetamnS=bet*(m-1)*n+betcA*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		AkkMinus1[i,i-1]=gam*n/thetamnS

def BuildAkkPlus1S(AkkPlus1,k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N):
	for i in range(1,k):
		m=k-i
		n=i
		thetamnS=bet*(m-1)*n+betcA*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		AkkPlus1[i-1,i-1]=alpR*(N-k)/thetamnS
		AkkPlus1[i-1,i]=sig*bet*(N-k)*n/thetamnS

def BuildbkSj(j,bSj,k,pIjPlus1,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M):
	for p in range(1,k):
		m=k-p
		n=p
		thetamnS=bet*(m-1)*n+betcA*n+gam*n+alpI*n+alpR*(N-m-n)+sig*bet*n*(N-m-n)
		if j==M-1:
			prod=1
		else:
			prod=pIjPlus1[k-1][n,0]
		bSj[k-1][n-1,0]=betcA*n*prod/thetamnS

def AlgorithmSj(j,invHS,PS,bS,pI,pS,AkkS,AkkMinus1S,AkkPlus1S,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M):
	k=2
	BuildbkSj(j,bS[j],k,pI[j],gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M)
	PS[j][k-1]=bS[j][k-1]
	for k in range(3,N+1):
		BuildbkSj(j,bS[j],k,pI[j],gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M)
		PS[j][k-1]=AkkMinus1S[k-2].dot(invHS[k-2].dot(PS[j][k-2]))+bS[j][k-1]
	k=N
	pS[j][k-1]=invHS[k-1].dot(PS[j][k-1])
	for k in reversed(range(2,N)):
		pS[j][k-1]=invHS[k-1].dot(AkkPlus1S[k-1].dot(pS[j][k])+PS[j][k-1])

		
def AlgorithmRj(j,invHR,PR,bR,pS,pI,pR,AkkR,AkkMinus1R,AkkPlus1R,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M):
	k=1
	BuildbkRj(j,bR[j-1],k,pS[j],pI[j],gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M)
	PR[j-1][k-1]=bR[j-1][k-1]
	for k in range(2,N):
		BuildbkRj(j,bR[j-1],k,pS[j],pI[j],gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M)
		PR[j-1][k-1]=AkkMinus1R[k-2].dot(invHR[k-2].dot(PR[j-1][k-2]))+bR[j-1][k-1]
	k=N-1
	pR[j-1][k-1]=invHR[k-1].dot(PR[j-1][k-1])
	for k in reversed(range(1,N-1)):
		pR[j-1][k-1]=invHR[k-1].dot(AkkPlus1R[k-1].dot(pR[j-1][k])+PR[j-1][k-1])
		
def AlgorithmIj(j,invHI,PI,bI,pS,pR,pI,AkkI,AkkMinus1I,AkkPlus1I,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M):
	k=1
	BuildbkIj(j,bI[j-1],k,pS[j],pR[j-1],gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
	PI[j-1][k-1]=bI[j-1][k-1]
	for k in range(2,N+1):
		BuildbkIj(j,bI[j-1],k,pS[j],pR[j-1],gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
		PI[j-1][k-1]=AkkMinus1I[k-2].dot(invHI[k-2].dot(PI[j-1][k-2]))+bI[j-1][k-1]
	k=N
	pI[j-1][k-1]=invHI[k-1].dot(PI[j-1][k-1])
	for k in reversed(range(1,N)):
		pI[j-1][k-1]=invHI[k-1].dot(AkkPlus1I[k-1].dot(pI[j-1][k])+PI[j-1][k-1])	

def p_reinfection_SIR(gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M,i0,s0,r0):
	

	HS=[np.zeros((k-1,k-1)) for k in range(1,N+1)]
	invHS=[np.zeros((k-1,k-1)) for k in range(1,N+1)]
	AkkS=[np.zeros((k-1,k-1)) for k in range(1,N+1)]
	AkkMinus1S=[np.zeros((k-1,k-2)) for k in range(2,N+1)]
	AkkPlus1S=[np.zeros((k-2,k-1)) for k in range(2,N+1)]

	k=2
	BuildAkkSI(AkkS[k-1],k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
	HS[k-1]=np.eye(k-1)-AkkS[k-1]
	invHS[k-1]=inv(HS[k-1])

	for k in range(3,N+1):
		BuildAkkSI(AkkS[k-1],k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
		BuildAkkMinus1SI(AkkMinus1S[k-2],k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
		BuildAkkPlus1SI(AkkPlus1S[k-2],k-1,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
		HS[k-1]=np.asmatrix(np.eye(k-1))-AkkS[k-1]-AkkMinus1S[k-2].dot(invHS[k-2].dot(AkkPlus1S[k-2]))
		invHS[k-1]=inv(HS[k-1])

	HR=[np.zeros((k,k)) for k in range(1,N+1)]
	invHR=[np.zeros((k,k)) for k in range(1,N+1)]
	AkkR=[np.zeros((k,k)) for k in range(1,N)]
	AkkMinus1R=[np.zeros((k,k-1)) for k in range(2,N)]
	AkkPlus1R=[np.zeros((k-1,k)) for k in range(2,N)]

	k=1
	BuildAkkR(AkkR[k-1],k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
	HR[k-1]=np.asmatrix(np.eye(k))-AkkR[k-1]
	invHR[k-1]=inv(HR[k-1])

	for k in range(2,N):
		BuildAkkR(AkkR[k-1],k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
		BuildAkkMinus1R(AkkMinus1R[k-2],k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
		BuildAkkPlus1R(AkkPlus1R[k-2],k-1,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
		HR[k-1]=np.asmatrix(np.eye(k))-AkkR[k-1]-AkkMinus1R[k-2].dot(invHR[k-2].dot(AkkPlus1R[k-2]))
		invHR[k-1]=inv(HR[k-1])

	HI=[np.zeros((k,k)) for k in range(1,N+1)]
	invHI=[np.zeros((k,k)) for k in range(1,N+1)]
	AkkI=[np.zeros((k,k)) for k in range(1,N+1)]
	AkkMinus1I=[np.zeros((k,k-1)) for k in range(2,N+1)]
	AkkPlus1I=[np.zeros((k-1,k)) for k in range(2,N+1)]

	k=1
	BuildAkkI(AkkI[k-1],k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
	HI[k-1]=np.asmatrix(np.eye(k))-AkkI[k-1]
	invHI[k-1]=inv(HI[k-1])

	for k in range(2,N+1):
		BuildAkkI(AkkI[k-1],k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
		BuildAkkMinus1I(AkkMinus1I[k-2],k,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
		BuildAkkPlus1I(AkkPlus1I[k-2],k-1,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N)
		HI[k-1]=np.eye(k)-AkkI[k-1]-AkkMinus1I[k-2].dot(invHI[k-2].dot(AkkPlus1I[k-2]))
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
	AlgorithmSj(j,invHS,PS,bS,pI,pS,AkkS,AkkMinus1S,AkkPlus1S,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M)
	for j in reversed(range(M-1)):
		AlgorithmRj(j+1,invHR,PR,bR,pS,pI,pR,AkkR,AkkMinus1R,AkkPlus1R,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M)
		AlgorithmIj(j+1,invHI,PI,bI,pS,pR,pI,AkkI,AkkMinus1I,AkkPlus1I,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M)
		AlgorithmSj(j,invHS,PS,bS,pI,pS,AkkS,AkkMinus1S,AkkPlus1S,gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M)

	return (N-1.0)/N*pS[0][s0+i0-1][i0-1,0]+1.0/N*pI[0][s0+i0-1][i0-1,0]
	

# gam=0.5
# gamA=gam
# bet=23.0/(2*284)
# betcA=bet
# betAc=bet
# alpI=0.0
# alpIA=0.0
# alpR=0.0
# alpRA=0.0
# sig=(1-0.98)
# sigA=sig

# N=284
# M=2

# i0=1
# s0=N-1
# r0=N-i0-s0
# start_time = time.time()
# probability=p_reinfection_SIR(gam,gamA,bet,betAc,betcA,alpI,alpIA,alpR,alpRA,sig,sigA,N,M,i0,s0,r0)
# elapsed_time = time.time() - start_time

# print(probability)
# print(elapsed_time)
