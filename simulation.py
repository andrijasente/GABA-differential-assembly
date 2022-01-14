# Author: Katerina Naydenova
# 
# Differential GABA-A receptor assembly diversifies structures and signaling
#
# Andrija Sente1*, Rooma Desai2, Katerina Naydenova1, Tomas Malinauskas3, 
# Youssef Jounaidi2, Jonas Miehling1, Xiaojuan Zhou2, Simonas Masiulis1,4, 
# Steven W. Hardwick5, Dimitri Y. Chirgadze5, Keith W. Miller2*, A. Radu Aricescu1*
#
# 1 - MRC Laboratory of Molecular Biology, Francis Crick Avenue, Cambridge, CB2 0QH, UK.
# 2 - Department of Anesthesia, Critical Care and Pain Medicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA, USA. 
# 3 - Division of Structural Biology, Wellcome Centre for Human Genetics, University of Oxford, Roosevelt Drive, Oxford, OX3 7BN, UK.
# 4 - Current address: Materials and Structural Analysis Division, Thermo Fisher Scientific, Achtseweg Noord, Eindhoven, 5651 GG, Netherlands.
# 5 - Department of Biochemistry, University of Cambridge, Tennis Court Road, Cambridge, CB2 1GA, UK.
#
# * - Correspondence to: asente@mrc-lmb.cam.ac.uk or k_miller@helix.mgh.harvard.edu or radu@mrc-lmb.cam.ac.uk

import time
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing.pool import Pool
from numba import jit, cuda

subunits = ['a4','b3','d']

all_possible_combinations = np.array(np.meshgrid(subunits,subunits,subunits,subunits,subunits)).T.reshape(-1,5)
double_conc = np.concatenate((all_possible_combinations, all_possible_combinations), axis=1)

subtype_labels = np.zeros(np.shape(all_possible_combinations)[0])
next_label = 1

for i in range(1,len(subtype_labels)):
	for j in range(0, i):
		if all_possible_combinations[i,:].tostring() in double_conc[j,:].tostring():
			subtype_labels[i] = subtype_labels[j]
			break
	else:
		subtype_labels[i] = next_label
		next_label += 1
		#print(all_possible_combinations[i,:], subtype_labels[i])

def SimulatePentamer(subunits, probs, weighted_interface_likelihoods_p, weighted_interface_likelihoods_pc):
	np_subunits = np.asarray(subunits, dtype=np.str_)
	pentamer = np.empty(5, dtype = object)
	pentamer[0] = np.random.choice(subunits,size=None,p=probs)
	for i in (1,2,3):
		principal = pentamer[i-1]
		which_principal = np.where(np_subunits == principal)[0]
		factor = weighted_interface_likelihoods_p[which_principal,:]
		mod_p = factor[0,:]
		pentamer[i] = np.random.choice(subunits,size=None,p=mod_p)
	principal = pentamer[3]
	complem = pentamer[0]
	which_principal = np.where(np_subunits == principal)[0]
	which_complem = np.where(np_subunits == complem)[0]
	factor_pc = weighted_interface_likelihoods_pc[which_complem,which_principal,:]
	mod_pc = factor_pc[0,:]
	pentamer[4] = np.random.choice(subunits,size=None,p=mod_pc)
	return (pentamer)	

def LoopOverAll(probs, subunits, interface_likelihoods, number_of_types, all_possible_combinations, subtype_labels):
	np.random.seed()	
	#important for randomness
	#start_time = time.time()
	counters = np.zeros(number_of_types)
	weighted_interface_likelihoods_p = (interface_likelihoods.T * probs[:, None]).T
	normaliz2 = 1/np.sum(weighted_interface_likelihoods_p, axis=1)
	weighted_interface_likelihoods_p = weighted_interface_likelihoods_p * normaliz2[:, None]
	weighted_interface_likelihoods_pc = np.stack((weighted_interface_likelihoods_p*interface_likelihoods[:,0],weighted_interface_likelihoods_p*interface_likelihoods[:,1],weighted_interface_likelihoods_p*interface_likelihoods[:,2]))
	normaliz3 = 1/np.sum(weighted_interface_likelihoods_pc, axis=2)
	weighted_interface_likelihoods_pc = weighted_interface_likelihoods_pc * normaliz3[:, :, None]
	for i in range(1,1000):
		pentamer = SimulatePentamer(subunits,probs,weighted_interface_likelihoods_p,weighted_interface_likelihoods_pc)
		ind_pentamer = np.where((all_possible_combinations == pentamer).all(axis=1))[0]
		subtype_number = int(subtype_labels[ind_pentamer])
		counters[subtype_number] += 1
		#print(pentamer, subtype_number, counters[subtype_number])
	#print(time.time()-start_time)
	return(*probs, *interface_likelihoods.flatten('C'), *counters)
	#plt.plot(counters,'ro')
	#plt.xlabel('subtype number')
	#plt.ylabel('occurence per 10000')
	#plt.show()

set_of_probs = []
max_transfection = 64.0
for i in np.geomspace(1/max_transfection,max_transfection,num=13):
	for j in np.geomspace(1/max_transfection,max_transfection,num=13):
		norm = float(1+i+j)
		set_of_probs.append(np.array([1/norm,i/norm,j/norm]))
#set_of_probs = [np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0]),np.array([1/3.0,1/3.0,1/3.0])]

start_time = time.time()
#interface_likelihoods_list = [np.array([[0.01,1,1],[1,1,1],[1,1,1]]),np.array([[0.01,1,1],[1,10,1],[1,1,1]]),np.array([[0.01,1,1],[1,100,1],[1,1,1]]),np.array([[0.01,1,1],[1,1000,1],[1,1,1]]),np.array([[0.01,1,1],[1,10000,1],[1,1,1]]),np.array([[0.01,1,1],[1,100000,1],[1,1,1]]),np.array([[0.01,10,1],[1,1,1],[1,1,1]]),np.array([[0.01,100,1],[1,1,1],[1,1,1]]),np.array([[0.01,1000,1],[1,1,1],[1,1,1]]),np.array([[0.01,10000,1],[1,1,1],[1,1,1]]),np.array([[0.01,100000,1],[1,1,1],[1,1,1]]),np.array([[0.01,1,1],[10,1,1],[1,1,1]]),np.array([[0.01,1,1],[100,1,1],[1,1,1]]),np.array([[0.01,1,1],[1000,1,1],[1,1,1]]),np.array([[0.01,1,1],[10000,1,1],[1,1,1]]),np.array([[0.01,1,1],[100000,1,1],[1,1,1]]),np.array([[0.01,10,1],[10,1,1],[1,1,1]]),np.array([[0.01,100,1],[100,1,1],[1,1,1]]),np.array([[0.01,1000,1],[1000,1,1],[1,1,1]]),np.array([[0.01,10000,1],[10000,1,1],[1,1,1]]),np.array([[0.01,100000,1],[100000,1,1],[1,1,1]])]
interface_likelihoods_list = []
#for aa in (0,0.1):
#	for ab in (1, 10, 100, 1000, 10000):
#		for ad in (1, 10, 100, 1000, 10000):
#			for ba in (1, 10, 100, 1000, 10000):
#				for bb in (1, 10, 100, 1000, 10000):
#					for bd in (1, 10, 100, 1000, 10000):
#						for da in (1, 10, 100, 1000, 10000):
#							for db in (1, 10, 100, 1000, 10000):
#								for dd in (1, 10, 100, 1000, 10000):
#									interface_likelihoods_list.append(np.array([[aa,ab,ad],[ba,bb,bd],[da,db,dd]]))

interface_likelihoods_list= [np.array([[0,1,1],[1,1,1],[1,1,1]])]

for interface_test in interface_likelihoods_list:
	#print(*interface_test)
	with Pool(processes=40) as pool:
		LoopingFunction = partial(LoopOverAll, subunits = subunits, interface_likelihoods=interface_test, number_of_types=next_label, all_possible_combinations=all_possible_combinations, subtype_labels=subtype_labels)
		results = pool.map(LoopingFunction, set_of_probs)
		for result in results:
			print(*result)
		pool.close()
#print(time.time()-start_time)
