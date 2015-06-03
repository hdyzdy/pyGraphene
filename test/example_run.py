#!/bin/env python

#EXAMPLE run for pyGraphene
#The below options will recreate the 1500/pc7/BaseCase runs from JPC paper
#Comments below indicate how to modify the script for different cases
#Environment variable $PYTHONPATH needs to include directory where pyGraphene classes are
#Environment variable $TINKERPATH needs to point to path where optimize executable and mm3.prm file are

import matplotlib
# matplotlib.use('Agg')
from pyGSimSet import *

simName = 'linear20_BC_1500'
#The starting substrate
#Change basetype and nrings as appropriate
pgn = pyGNetwork(baseType='linear',nRings=20)

#Excluded reactions: numbers match those from JPC paper
#exr = [14,15,16,17,18,19,20,22,29,30,31,32,33,34,35,36,37,38,39,41] #E1
#exr = [35,36,37,38] #E2
#exr = [29,30,31,32] #E3
#exr = [29,30,31,32,35,36,37,38] #E4
exr = [44]  #BC
rxns = pyGRLoad(excludedReactions=exr)
# print rxns
#Environment
pyge = pyGEnvironment(temp=1500)
#Eamples for other F cases
#pyge = pyGEnvironment(temp=1500,xOfC2H2=0.01) #F1 1500 K
#pyge = pyGEnvironment(temp=2000,xOfH2=0.01)   #F2 2000 K
#pyge = pyGEnvironment(temp=2500,xOfH=0.001)   #F3 2500 K
#pyge = pyGEnvironment(temp=1500,xOfCH3=0.01)  #F5 1500 K

#The simulation that will be run with different seeds
pgs = pyGSimulation(network=pgn,envs=[pyge],reactions=rxns)

#The set of simulations:
#seeds is a list of starting random seeds
#baseSim is the starting simulation defined above
#name is used for outputting filenames and is also defined above
pgss = pyGSimSet(seeds=range(15),baseSim=pgs,name=simName)

startTime = time.time()  #Keep track of how long simulations take
#Run simulations
#endTime is the cut off time for simulations
#first nOpt runs will be optimized with tinker the rest will not
#if saveSims==True then simulation variables will be saved
#if a run is optimized then an xyz file will be output for it
pgss.run(endTime=0.001,nOpt=15,saveSims=True) 

totalTime = time.time() - startTime
print 'Total time: ', totalTime/60., ' minutes'
#Write out pgss variable
pgss.save()

#Write out simulation mean quantities to csv for plotting in matlab, excel, kaleidagraph, ...
pgss.toCSV()

#Plot mean quantities using pyplot and save to disk
plt.plot(pgss.t_bins,pgss.nCMean,linewidth=4.0,color='black')
gr = pgss.calcGrowthRate(nAvgSteps=None)
fr5 = pgss.calcFR5()
plt.savefig(simName+'_growth')
plt.clf()
# plt.semilogy(pgss.t_bins,gr,pgss.t_bins,gr)
total_num=0
tmp_sum_gr=0
for (index, tmp_gr) in enumerate(gr):
	if not numpy.isnan(tmp_gr):
		tmp_sum_gr+=tmp_gr
		total_num+=1

mean_gr=tmp_sum_gr/float(total_num)
array_mean_gr=deepcopy(pgss.t_bins)
for i in range(len(array_mean_gr)):
	array_mean_gr[i]=mean_gr

# print array_mean_gr
plt.semilogy(pgss.t_bins,gr,'g')
# print len(pgss.t_bins),len(array_mean_gr)
plt.semilogy(pgss.t_bins,array_mean_gr,'r')
# print pgss.t_bins,gr
plt.savefig(simName+'_growth_rate')
plt.clf()
plt.plot(pgss.t_bins,fr5)
plt.savefig(simName+'_fr5')

plt.clf()
plt.savefig(simName+'_Rxn')
pgss.plotRxnBar()
plt.show
print 'Done.'
