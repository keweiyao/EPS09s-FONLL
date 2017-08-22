#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys, os, h5py
from subprocess import Popen, PIPE, STDOUT, call
import fortranformat as ff 
from itertools import repeat
from multiprocessing import Pool, cpu_count
line = ff.FortranRecordReader('(F20.0)')


Adict = {'1': 'p',	'2':'d', '197': 'Au', '208': 'Pb'}

# 2-31 are errorsets
def fonll( sqrts=2760, A=[208, 208], b=[0.0, 0.0], 
			pT=10.0, y=0.0,
			M=1.3, scale=1.0, 
			EPS09_errorset = 1, 
			order = 'FONLL'):
	order_index = 1 if order == 'FONLL' else 2
	tag = Adict[str(A[0])] + '-' + Adict[str(A[1])]
	print ("system: ", tag)
	if not os.path.exists("log"):
		os.makedirs("log")
	call('rm -r ./log/{}-{}*'.format(tag, y), shell=True)
	beam1 = A[0] if A[0]==1 else -A[0]*1000-EPS09_errorset
	beam2 = A[1] if A[1]==1 else -A[1]*1000-EPS09_errorset
	b1, b2 = b
	
	print ("pt = {} [GeV], y = {}".format(pT, y) )
	inputs = """{}
{} {} {} 0 0 10550 
{} {} {} 0 0 10550
{}
-1
{} {}
{} {}
{}
""".format(	"{}-{}".format(tag, y), 
			beam1, b1, 0.5*sqrts, 
			beam2, b2, 0.5*sqrts, 
			M, 
			scale, scale, 
			pT, y, 
			order_index)
			
	p = Popen(["./fonlllha"], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
	p.communicate(input=inputs.encode('utf-8'))
	# read out fonll results
	dsigma_dpT2_dy = 0.0
	with open('./log/{}-{}.outlog'.format(tag, y), 'r') as f:
		for l in f:
			nl = l.split()
			if len(nl) < 1:
				continue
			if nl[0] == 'pt,y,fonll':	
				dsigma_dpT2_dy = line.read(nl[3])[0]*1e-9 # from pb/GeV2 to mb/Gev2
	return {"y": y, "pT": pT, "X": dsigma_dpT2_dy}

def calculate_at_y(y, sqrts, A, B, b1, b2, pT, M):
	return fonll(  sqrts=sqrts, A=[A, B], b=[b1, b2], 
				pT=pT, y=y,
				M=M, scale=1.0, EPS09_errorset=1, order='FONLL')

def calculate_at_b1_b2(M, sqrts, A, B, b1, b2, NpT):	
	pT = M*np.linspace(0.1**0.25, 100**0.25, NpT)**4
	mT = np.sqrt(pT**2 + M**2)
	Ny = 7
	ymax = np.arccosh(sqrts/2./mT)
	print("b1 = {}, b2 = {}".format(b1, b2))
	X = np.zeros([Ny, pT.shape[0]])
	for pTindex, (ipT, iymax) in enumerate(zip(pT, ymax)):
		yarray = np.linspace(0., iymax, Ny)
		Njobs = yarray.shape[0]
		Nproc = cpu_count()-1 or 1
		print("# pT, Njobs, Nproc = ", ipT, Njobs, Nproc)
		with Pool(processes=Nproc) as pool:
			all_ds = pool.starmap(calculate_at_y, 
							zip(yarray, repeat(sqrts),
								repeat(A), repeat(B),
								repeat(b1), repeat(b2),
								repeat(ipT), repeat(M))
								)
		
		for ds in all_ds:
			y = ds['y']
			yindex = int( (y - yarray[0])/(yarray[1]-yarray[0]) )
			X[yindex, pTindex] = ds['X']

	return ymax, pT, X

def main():
	M, sqrts = 1.3, 5020
	A, B = 208, 208
	tag = Adict[str(A)] + '-' + Adict[str(B)]
	NpT, Ny_half = 20, 7
	bmin, bmax, Nb = 0.0, 10.0, 11
	barray = np.linspace(bmin, bmax, Nb)

	# The averaged dsigam_dpT_dy per collison
	ymax, pT, X0 = calculate_at_b1_b2(M, sqrts, A, B, -1, -1, NpT)

	# The averaged dsigam_dpT_dy per collison
	ymax, pTfiner, Xfiner = calculate_at_b1_b2(M, sqrts, A, B, -1, -1, NpT*3)

	# pp-ref at the same energy and pT point
	ymax, pTfiner, Xpp = calculate_at_b1_b2(M, sqrts, 1, 1, -1, -1, NpT*3)

	# This b1,b2-dependent dsigma_dpTdy per collision dims = [Nb1, Nb2, Ny, NpT]
	# This grid is coarse
	Rarray = np.zeros([Nb, Nb, 2*Ny_half-1, NpT])
	for i1, b1 in enumerate(barray):
		for i2, b2 in enumerate(barray):
			ymax, pT, X = calculate_at_b1_b2(M, sqrts, A, B, b1, b2, NpT)
			Rarray[i1, i2, 6:] = X/X0
	for i1 in range(Nb):
		for i2 in range(Nb):
			cut = Rarray[i1, i2, 7:]
			Rarray[i1, i2, :6] = cut[::-1]
	
	# save to file
	path = './spectra/'
	filename = 'lhc-PbPb-b1-b2-y-pT.hdf5'
	if not os.path.exists(path):
		os.makedirs(path)
	f = h5py.File(path+filename, 'a')
	gp = f.create_group("{}-{}".format(sqrts,tag))
	gp.create_dataset('b1', data=barray)
	gp.create_dataset('b2', data=barray)
	gp.create_dataset('ymax', data=ymax)
	gp.create_dataset('pT', data=pT)
	gp.create_dataset('pTfiner', data=pTfiner)
	gp.attrs.create('M', M)
	gp.attrs.create('A', A)
	gp.attrs.create('B', B)
	gp.attrs.create('sqrts', sqrts)
	gp.create_dataset('X-IP/avg', data=Rarray)
	gp.create_dataset('X-avg', data=Xfiner)
	gp.create_dataset('X-pp', data=Xpp)
	f.close()


if __name__ == '__main__':
	main()

