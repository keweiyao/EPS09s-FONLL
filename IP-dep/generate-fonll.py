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
def fonll( sqrts, A, b, 
			pT, y,
			M, scale, 
			EPS09_errorset, 
			order, PDF):
	order_index = 1 if order == 'FONLL' else 2
	tag = Adict[str(A[0])] + '-' + Adict[str(A[1])]
	print ("system: ", tag)
	if not os.path.exists("log"):
		os.makedirs("log")
	call('rm -r ./log/{}-{}*'.format(tag, y), shell=True)
	beam1 = A[0] if A[0]==1 else -A[0]*1000-EPS09_errorset
	beam2 = A[1] if A[1]==1 else -A[1]*1000-EPS09_errorset
	b1, b2 = b
	
	print ("pt = {} [GeV], y = {}, PDF = {}".format(pT, y, PDF) )
	inputs = """{}
{} {} {} 0 0 {:d}
{} {} {} 0 0 {:d}
{}
-1
{} {}
{} {}
{}
""".format(	"{}-{}".format(tag, y), 
			beam1, b1, 0.5*sqrts, PDF,
			beam2, b2, 0.5*sqrts, PDF,
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

def calculate_at_y(y, sqrts, A, B, b1, b2, pT, M, tPDF):
	print(tPDF)
	return fonll(  sqrts=sqrts, A=[A, B], b=[b1, b2], 
					pT=pT, y=y,
					M=M, scale=1.0, EPS09_errorset=1, 
					order='FONLL', PDF=tPDF)

def calculate_at_b1_b2(M, sqrts, A, B, b1, b2, NpT, PDF):	
	pT = M*np.linspace(0.1, 100, NpT)
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
								repeat(ipT), repeat(M),
								repeat(PDF)	)
								)
		
		for ds in all_ds:
			y = ds['y']
			yindex = int( (y - yarray[0])/(yarray[1]-yarray[0]) )
			X[yindex, pTindex] = ds['X']

	return ymax, pT, X

def main():
	M, sqrts = 1.3, 2760
	A, B = 208, 208
	tag = Adict[str(A)] + '-' + Adict[str(B)]
	NpT, Ny_half = 100, 7
	f = h5py.File("spectra0.hdf5",'a')
	XAB = np.zeros([2*Ny_half-1, NpT])
	Xpp = np.zeros([2*Ny_half-1, NpT])
	# The averaged dsigam_dpT_dy per collison -- EPS09 + CTEQ6.6 / CTEQ6.6
	#if 'nCTEQ15np' in f:
	#	del f['nCTEQ15np']
	#gp0 = f.create_group("nCTEQ15np")
	#gp = gp0.create_group("Pb-Pb-5020")
	gp0 = f["nCTEQ15np"]
	del f["nCTEQ15np/Pb-Pb-2760"]
	gp = gp0.create_group("Pb-Pb-2760")

	ymax, pT, XAB_half = calculate_at_b1_b2(M, sqrts, 1, 1, -1, -1, NpT, 105100)
	XAB[Ny_half:] = XAB_half[:-1]
	XAB[:Ny_half] = XAB_half[::-1]
	gp.create_dataset('AB', data=XAB)

	ymax, pT, Xpp_half = calculate_at_b1_b2(M, sqrts, 1, 1, -1, -1, NpT, 104000)
	Xpp[Ny_half:] = Xpp_half[:-1]
	Xpp[:Ny_half] = Xpp_half[::-1]
	gp.create_dataset('pp', data=Xpp)
	

	f.close()

if __name__ == '__main__':
	main()

