import numpy as np
import h5py
import hdf5storage 
import pickle
import pandas as pd
import os
import re

path = "/mnt/ratmemory_hippocampus/Ito/data/210108_yuta2/"

f = open(path+"shuffle_parameters.txt")

lines = f.readlines()
f.close()

names = ["TC","fixed_music2r","fixed_music3","fixed_music4s","fixed_music4L","fixed_music_ex2","periodic_click"]
dnum = [100,6,6,4,4,8,44]
dnum2 = [1,10,10,10,3,10,5]

dindex = []
dstrt = []
data = {}

for n,x in enumerate(lines):
	if len(x) <= 1:
		continue
	temp = [1 if y in x else 0 for y in names]
	if sum(temp) > 0:
		index = np.argmax(temp)
		dstrt_ = []
		dstrt = []
		continue
	else:
		tmp = [m.rstrip() for m in re.split('[ ,]', x) if len(m) >0]
		if index == 0:
			dstrt.extend(tmp)
			tdst = np.zeros((1,len(dstrt)))
			tdst[0,:] = dstrt
			if not names[index] in data:
				data[names[index]] = tdst
			else:
				data[names[index]+"_2"] = tdst
		else:
			if len(dstrt_) == 0:
				dstrt_.extend(tmp[1:])
			else:
				dstrt_.extend(tmp)
			if len(dstrt_) >= dnum[index]:
				dstrt.extend(dstrt_)
				dstrt_ = []
			if len(dstrt) >= dnum[index]*dnum2[index]:
				tdst = np.zeros((1,len(dstrt)))
				tdst[0,:] = dstrt
				if not names[index] in data:
					data[names[index]] = tdst
				else:
					data[names[index]+"_2"] = tdst
for x,y in data.items():
	print(x,len(y[0]))

hdf5storage.write(data, '.', path+"shuf_param.mat", matlab_compatible=True)



