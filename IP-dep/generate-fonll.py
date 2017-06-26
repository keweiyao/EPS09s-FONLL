#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from subprocess import Popen, PIPE, STDOUT, call
import fortranformat as ff
import h5py
line = ff.FortranRecordReader('(F20.0)')


Adict = {'1': 'p',	'2':'d', '197': 'Au', '208': 'Pb'}

# 2-31 are errorsets
def fonll( sqrts=2760, A=[208, 208], b=[0.0, 8.0], 
			pTrange=[0.1, 100, 20], yrange=[-2,2,5],
			M=1.3, scale=1.0, EPS09_errorset = 1, order = 'FONLL'):
	order_index = 1 if order == 'FONLL' else 2
	tag = Adict[str(A[0])] + '-' + Adict[str(A[1])]
	print ("system: ", tag)
	call('rm -r ./log/%s*'%tag, shell=True)
	beam1 = A[0] if A[0]==1 else -A[0]*1000-EPS09_errorset
	beam2 = A[1] if A[1]==1 else -A[1]*1000-EPS09_errorset
	b1, b2 = b
	pT = np.linspace(*pTrange)
	y = np.linspace(*yrange)

	for iy in y:
		for ipT in pT:
			print ("pt = {} [GeV], y = {}".format(ipT, iy) )
			inputs = """{}
{} {} {} 0 0 10550 
{} {} {} 0 0 10550
{}
-1
{} {}
{} {}
{}
""".format(	tag, 
			beam1, b1, 0.5*sqrts, 
			beam2, b2, 0.5*sqrts, 
			M, 
			scale, scale, 
			ipT, iy, 
			order_index)
			#print(inputs)
			p = Popen(["./fonlllha"], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
			p.communicate(input=inputs.encode('utf-8'))
	# read out fonll results
	shape = [yrange[-1], pTrange[-1]]
	pT, y, dsigma_dpT2_dy = [],[],[]
	with open('./log/%s.outlog'%tag, 'r') as f:
		for l in f:
			nl = l.split()
			if len(nl) < 1:
				continue
			if nl[0] == 'pt,y,fonll':	
				pT.append(line.read(nl[1])[0])
				y.append(line.read(nl[2])[0])
				dsigma_dpT2_dy.append(line.read(nl[3])[0]*1e-9) # from pb/GeV2 to mb/Gev2
	pT = np.array(pT).reshape(*shape)
	y = np.array(y).reshape(*shape)
	dsigma_dpT2_dy = np.array(dsigma_dpT2_dy).reshape(*shape)
	# save to file
	path = './spectra/'
	filename = 'temp.hdf5'
	if not os.path.exists(path):
		os.makedirs(path)
	f = h5py.File(path+filename, 'a')
	ds = f.create_group("{}-{}-{}-{}".format(sqrts,tag,b1,b2))
	ds.attrs.create("b1b2", np.array([b1, b2]))
	ds.create_dataset('pT', data=pT)
	ds.create_dataset('y', data=y)
	ds.create_dataset('dsigma_dpT2_dy', data=dsigma_dpT2_dy)
	f.close()


def main():
	db = 0.5
	for i in range(17):
		b1 = i*db
		for j in range(i,17):
			b2 = j*db
			fonll(  sqrts=200, A=[197, 197], b=[b1, b2], 
				pTrange=[0.5, 30, 20], yrange=[-0,0,1],
				M=1.3, scale=1.0, EPS09_errorset = 1, order = 'FONLL')


if __name__ == '__main__':
	main()

