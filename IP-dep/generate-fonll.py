#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from subprocess import Popen, PIPE, STDOUT, call
import fortranformat as ff
import h5py

line = ff.FortranRecordReader('(F20.0)')


M = 1.3
scale = 1.
sqrts = 2760
A1 = 208
A2 = 208
ptmin = 0.1
ptmax = 10.0
npt = 2
y = 0.0
b1 = 0.0
b2 = 9.0
EPS09_errorset = 1 # 2-31 are errorsets

order = 'FONLL'
if order == 'FONLL':
    order_index = 1
if order == 'FO':
    order_index = 2

tag = ''
if A1 == 1:
	tag += 'p'
if A1 == 197:
	tag += 'Au'
if A1 == 208:
	tag += 'Pb'
if A2 == 1:
	tag += 'p'
if A2 == 197:
	tag += 'Au'
if A2 == 208:
	tag += 'Pb'
print ("system: ", tag)


if A1 != 1:
	beam1 = -A1*1000-EPS09_errorset
else:
	beam1 = 1
if A2 != 1:
	beam2 = -A2*1000-EPS09_errorset
else:
	beam2 = 1


call('rm -r ./%s*'%tag, shell=True)

pt_list = np.linspace(ptmin, ptmax, npt)
for pt in pt_list:
	print ("pt, y = ", pt, y)
	inputs = "{} \n {} {} {} 0 0 10550 \n {} {} {} 0 0 10550 \n {} \n -1 \n {} {} \n {} {} \n {}\n".format(
		tag, 
		beam1, b1, 0.5*sqrts, 
		beam2, b2, 0.5*sqrts, 
		M, 
		scale, scale, 
		pt, y, 
		order_index)
	p = Popen(["./fonlllha"], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
	p.communicate(input=inputs.encode('utf-8'))

pt = []
y = []
dsigma1 = []
f = open('%s.outlog'%tag, 'r')
for l in f:
	nl = l.split()
	if len(nl) < 1:
		continue
	if nl[0] == 'pt,y,fonll':	
		pt.append(line.read(nl[1])[0])
		y.append(line.read(nl[2])[0])
		dsigma1.append(line.read(nl[3])[0]*1e-9) # from pb/GeV2 to mb/Gev2
f.close()
dsigma1 = np.array(dsigma1)
path = './spectra/'.format(order, sqrts)
filename = 'spectra-IP.hdf5'
if not os.path.exists(path):
	os.makedirs(path)
f = h5py.File(path+filename, 'a')
ds = f.create_group("{}-{}-4".format(sqrts,tag))
ds.attrs.create("b1b2", np.array([b1, b2]))
ds.create_dataset('pT', data=pt)
ds.create_dataset('dsigma_dpt2_dy', data=dsigma1)
f.close()

