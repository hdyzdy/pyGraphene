#!/bin/env python

from pyGSimulation import *

class pyGSimSet():
    def __init__(self,baseSim=pyGSimulation(),seeds=range(15),name='pgss_default'):
        Nbins = 1000
        self.name = name
        self.seeds = seeds
        self.baseSim = baseSim
        self.nRan = 0
        self.nCMean = numpy.zeros(Nbins,float)
        self.r5Mean = numpy.zeros(Nbins,float)
        self.r6Mean = numpy.zeros(Nbins,float)
        self.nCmfe5Mean = numpy.zeros(Nbins,float)
        self.site_acMean = numpy.zeros(Nbins,float)
        self.site_ac5Mean = numpy.zeros(Nbins,float)
        self.site_zzMean = numpy.zeros(Nbins,float)
        self.site_zz5Mean = numpy.zeros(Nbins,float)
        self.site_feMean = numpy.zeros(Nbins,float)
        self.site_fe5Mean = numpy.zeros(Nbins,float)
        self.site_bay5Mean = numpy.zeros(Nbins,float)
        self.site_bay5lMean = numpy.zeros(Nbins,float)
        self.site_bay6Mean = numpy.zeros(Nbins,float)
        self.site_bay6lMean = numpy.zeros(Nbins,float)
        self.site_otherMean = numpy.zeros(Nbins,float)
        self.rxnCounts = numpy.zeros(len(baseSim.reactions),float)

    def run(self,endTime=0.001,nOpt=numpy.Inf,saveSims=False,plotNC=False):
        self.t_bins = numpy.linspace(0,endTime,num=len(self.nCMean))
        for k in range(len(self.seeds)):
            print 'Performing run with seed: ', self.seeds[k]
            if k + 1 <= nOpt:
                doOpt = True
                outputGeos = True
            else:
                doOpt = False
                outputGeos = False
            runSim = deepcopy(self.baseSim)
            random.seed(self.seeds[k])
            runSim.name = self.name+'_sim'+'{0:03d}'.format(k)
            runSim.randState = random.getstate()
            runSim.reactApplyDict = self.baseSim.reactApplyDict #All sims share same reactApplyDict
            runSim.run(endMethod='time',endValue=endTime,doOpt=doOpt,outputGeos=outputGeos)
            self.baseSim.reactApplyDict = runSim.reactApplyDict
            if plotNC: plt.plot(runSim.time,runSim.nC,linewidth=0.2)
            if saveSims:
                runSim.reactApplyDict = None 
                runSim.save()
            self.updateMeans(runSim)
            self.nRan += 1

    def updateMeans(self,pgs):
        Nbins = len(self.t_bins)
        #Initialize
        nCBinned = numpy.zeros(Nbins,float)  #Carbon atoms
        r5Binned = numpy.zeros(Nbins,float)  #Five-member rings
        r6Binned = numpy.zeros(Nbins,float) #Six-member rings
        nCmfe5Binned = numpy.zeros(Nbins,float) #Adjusted carbons
        site_acBinned = numpy.zeros(Nbins,float)
        site_ac5Binned = numpy.zeros(Nbins,float)
        site_zzBinned = numpy.zeros(Nbins,float)
        site_zz5Binned = numpy.zeros(Nbins,float)
        site_feBinned = numpy.zeros(Nbins,float)
        site_fe5Binned = numpy.zeros(Nbins,float)
        site_bay5Binned = numpy.zeros(Nbins,float)
        site_bay5lBinned = numpy.zeros(Nbins,float)
        site_bay6Binned = numpy.zeros(Nbins,float)
        site_bay6lBinned = numpy.zeros(Nbins,float)
        site_otherBinned = numpy.zeros(Nbins,float)
        rxnCounts = pgs.calcRxnCounts()
        #Starting value
        nCBinned[0] += pgs.nC[0]
        r5Binned[0] += pgs.nR5[0]
        r6Binned[0] += pgs.nR6[0]
        nCmfe5Binned[0] += pgs.nC[0] - 2*pgs.siteCounts[0]['fe5']
        site_acBinned[0] += pgs.siteCounts[0]['ac']
        site_ac5Binned[0] += pgs.siteCounts[0]['ac5']
        site_zzBinned[0] += pgs.siteCounts[0]['zz']
        site_zz5Binned[0] += pgs.siteCounts[0]['zz5']
        site_feBinned[0] += pgs.siteCounts[0]['fe']
        site_fe5Binned[0] += pgs.siteCounts[0]['fe5']
        site_bay5Binned[0] += pgs.siteCounts[0]['bay5']
        site_bay5lBinned[0] += pgs.siteCounts[0]['bay5l']
        site_bay6Binned[0] += pgs.siteCounts[0]['bay6']
        site_bay6lBinned[0] += pgs.siteCounts[0]['bay6l']
        site_otherBinned[0] += pgs.siteCounts[0]['other']
        #Bin
        it = 0
        it_start = it
        for m in range(1,len(self.t_bins)):
            while pgs.time[it] <= self.t_bins[m]:
                nCBinned[m] += pgs.nC[it]
                r5Binned[m] += pgs.nR5[it]
                r6Binned[m] += pgs.nR6[it]
                nCmfe5Binned[m] += pgs.nC[it] - 2*pgs.siteCounts[it]['fe5']
                site_acBinned[m] += pgs.siteCounts[it]['ac']
                site_ac5Binned[m] += pgs.siteCounts[it]['ac5']
                site_zzBinned[m] += pgs.siteCounts[it]['zz']
                site_zz5Binned[m] += pgs.siteCounts[it]['zz5']
                site_feBinned[m] += pgs.siteCounts[it]['fe']
                site_fe5Binned[m] += pgs.siteCounts[it]['fe5']
                site_bay5Binned[m] += pgs.siteCounts[it]['bay5']
                site_bay5lBinned[m] += pgs.siteCounts[it]['bay5l']
                site_bay6Binned[m] += pgs.siteCounts[it]['bay6']
                site_bay6lBinned[m] += pgs.siteCounts[it]['bay6l']
                site_otherBinned[m] += pgs.siteCounts[it]['other']
                it += 1
                if it == len(pgs.time):
                    break
            numVals = it - it_start
            it_start = it
            if numVals == 0:
                nCBinned[m] = nCBinned[m-1]
                r5Binned[m] = r5Binned[m-1]
                r6Binned[m] = r6Binned[m-1]
                nCmfe5Binned[m] = nCmfe5Binned[m-1] 
                site_acBinned[m] = site_acBinned[m-1]
                site_ac5Binned[m] = site_ac5Binned[m-1]
                site_zzBinned[m] = site_zzBinned[m-1]
                site_zz5Binned[m] = site_zz5Binned[m-1]
                site_feBinned[m] = site_feBinned[m-1]
                site_fe5Binned[m] = site_fe5Binned[m-1]
                site_bay5Binned[m] = site_bay5Binned[m-1]
                site_bay5lBinned[m] = site_bay5lBinned[m-1]
                site_bay6Binned[m] = site_bay6Binned[m-1]
                site_bay6lBinned[m] = site_bay6lBinned[m-1]
                site_otherBinned[m] = site_otherBinned[m-1]
            else:
                nCBinned[m] /= numVals
                r5Binned[m] /= numVals
                r6Binned[m] /= numVals
                nCmfe5Binned[m] /= numVals
                site_acBinned[m] /= numVals
                site_ac5Binned[m] /= numVals
                site_zzBinned[m] /= numVals
                site_zz5Binned[m] /= numVals
                site_feBinned[m] /= numVals
                site_fe5Binned[m] /= numVals
                site_bay5Binned[m] /= numVals
                site_bay5lBinned[m] /= numVals
                site_bay6Binned[m] /= numVals
                site_bay6lBinned[m] /= numVals
                site_otherBinned[m] /= numVals
            if it == len(pgs.time):
                break
        if m < len(self.t_bins)-1:
            self.t_bins[m+1:] =  numpy.nan
            self.nCMean[m+1:] =  numpy.nan
            self.r5Mean[m+1:] =  numpy.nan
            self.r6Mean[m+1:] =  numpy.nan
            self.nCmfe5Mean[m+1:] =  numpy.nan
            self.site_acMean[m+1:] =  numpy.nan
            self.site_ac5Mean[m+1:] =  numpy.nan
            self.site_zzMean[m+1:] =  numpy.nan
            self.site_zz5Mean[m+1:] =  numpy.nan
            self.site_feMean[m+1:] =  numpy.nan
            self.site_fe5Mean[m+1:] =  numpy.nan
            self.site_bay5Mean[m+1:] =  numpy.nan
            self.site_bay5lMean[m+1:] =  numpy.nan
            self.site_bay6Mean[m+1:] =  numpy.nan
            self.site_bay6lMean[m+1:] =  numpy.nan
            self.site_otherMean[m+1:] =  numpy.nan
        self.nCMean = 1.0/(self.nRan+1)*((self.nRan)*self.nCMean + nCBinned)
        self.r5Mean = 1.0/(self.nRan+1)*((self.nRan)*self.r5Mean + r5Binned)
        self.r6Mean = 1.0/(self.nRan+1)*((self.nRan)*self.r6Mean + r6Binned)
        self.nCmfe5Mean = 1.0/(self.nRan+1)*((self.nRan)*self.nCmfe5Mean + nCmfe5Binned)
        self.site_acMean = 1.0/(self.nRan+1)*((self.nRan)*self.site_acMean + site_acBinned)
        self.site_ac5Mean = 1.0/(self.nRan+1)*((self.nRan)*self.site_ac5Mean + site_ac5Binned)
        self.site_zzMean = 1.0/(self.nRan+1)*((self.nRan)*self.site_zzMean + site_zzBinned)
        self.site_zz5Mean = 1.0/(self.nRan+1)*((self.nRan)*self.site_zz5Mean + site_zz5Binned)
        self.site_feMean = 1.0/(self.nRan+1)*((self.nRan)*self.site_feMean + site_feBinned)
        self.site_fe5Mean = 1.0/(self.nRan+1)*((self.nRan)*self.site_fe5Mean + site_fe5Binned)
        self.site_bay5Mean = 1.0/(self.nRan+1)*((self.nRan)*self.site_bay5Mean + site_bay5Binned)
        self.site_bay5lMean = 1.0/(self.nRan+1)*((self.nRan)*self.site_bay5lMean + site_bay5lBinned)
        self.site_bay6Mean = 1.0/(self.nRan+1)*((self.nRan)*self.site_bay6Mean + site_bay6Binned)
        self.site_bay6lMean = 1.0/(self.nRan+1)*((self.nRan)*self.site_bay6lMean + site_bay6lBinned)
        self.site_otherMean = 1.0/(self.nRan+1)*((self.nRan)*self.site_otherMean + site_otherBinned)
        self.rxnCounts = 1.0/(self.nRan+1)*((self.nRan)*self.rxnCounts + rxnCounts)

    def calcGrowthRate(self,nAvgSteps=None):
        if nAvgSteps is None:
            nSteps = len(self.t_bins)
            nAvgSteps = int(numpy.floor(nSteps/10))
        gr = numpy.zeros(len(self.nCMean),float)*numpy.nan  #Growth rate
        for k in range(nAvgSteps/2,len(self.nCMean)-nAvgSteps/2):
            #nCmfe5 = self.nCMean - self.site_fe5Mean
            gr[k] = (self.nCMean[k+nAvgSteps/2] - self.nCMean[k-nAvgSteps/2]) / \
                    (self.t_bins[k+nAvgSteps/2] - self.t_bins[k-nAvgSteps/2])
        return gr

    def calcGrowthRate_mfe5(self,nAvgSteps=None):
        if nAvgSteps is None:
            nSteps = len(self.t_bins)
            nAvgSteps = int(numpy.floor(nSteps/10))
        gr = numpy.zeros(len(self.nCmfe5Mean),float)*numpy.nan  #Growth rate
        for k in range(nAvgSteps/2,len(self.nCmfe5Mean)-nAvgSteps/2):
            #nCmfe5 = self.nCMean - self.site_fe5Mean
            gr[k] = (self.nCmfe5Mean[k+nAvgSteps/2]-self.nCmfe5Mean[k-nAvgSteps/2]) / \
                    (self.t_bins[k+nAvgSteps/2] - self.t_bins[k-nAvgSteps/2])
        return gr

    def calcFR5(self):
        fr5 = 1/(0.375)*(self.r5Mean - self.site_fe5Mean) / \
                (self.r5Mean - self.site_fe5Mean + self.r6Mean)   
        return fr5

    def calcNBF(self): #Edge bifurcations as defined in JPC paper
        return self.site_bay6lMean + self.site_otherMean

    def plotRxnBar(self):
        left = numpy.arange(0.5,len(self.rxnCounts),1)
        plt.bar(left,self.rxnCounts,width=1.0,log=True)
        plt.xticks(left+0.5,map(str,range(1,len(self.rxnCounts)+1)),rotation=90)

    def save(self):
        pFile = open(self.name+'.p','w')
        cPickle.dump(self,pFile)
        pFile.close()

    def toCSV(self):
        csvFile = open(self.name+'.csv','w')
        gr = self.calcGrowthRate()
        fr5 = self.calcFR5()
        nbf = self.calcNBF()
        csvFile.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18}\n'.format(
                "Time","nC","r5","r6","nCmfe5","ac","ac5","zz","zz5","fe","fe5","bay5","bay5l","bay6","bay6l","other",
                "growth rate","fR5","bifurcations"))
        for k in range(len(self.t_bins)):
            csvFile.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18}\n'.format(
                self.t_bins[k],
                self.nCMean[k],
                self.r5Mean[k],
                self.r6Mean[k],
                self.nCmfe5Mean[k],
                self.site_acMean[k],
                self.site_ac5Mean[k],
                self.site_zzMean[k],
                self.site_zz5Mean[k],
                self.site_feMean[k],
                self.site_fe5Mean[k],
                self.site_bay5Mean[k],
                self.site_bay5lMean[k],
                self.site_bay6Mean[k],
                self.site_bay6lMean[k],
                self.site_otherMean[k],
                gr[k],
                fr5[k],
                nbf[k]))
        csvFile.close()

if __name__ == '__main__':
    pgn = pyGNetwork(baseType='pc',nRings=7)
    pgs = pyGSimulation(network=pgn,envs=[pyGEnvironment(temp=2000)])
    pgss = pyGSimSet(seeds=range(15),baseSim=pgs)
    startTime = time.time()
    pgss.run(endTime=0.0005,nOpt=0)
    totalTime = time.time() - startTime
    print 'Total time: ', totalTime/60., ' minutes'
    pgss.save()
    pgss.toCSV()
    #Plotting
    #for k in range(len(pgss.sims)):
#        plt.plot(pgss.sims[k].time,pgss.sims[k].nC,linewidth=0.2)
    #plt.figure()
    plt.plot(pgss.t_bins,pgss.nCmfe5Mean,linewidth=4.0,color='black')
    gr = pgss.calcGrowthRate(nAvgSteps=None)
    grmfe5 = pgss.calcGrowthRate_mfe5(nAvgSteps=None)
    plt.figure()
    plt.semilogy(pgss.t_bins,gr,pgss.t_bins,grmfe5)
    plt.figure()
    pgss.plotRxnBar()
    plt.show()
    #plt.savefig('simset.png')

    #pgss.sims[-1].networks[-1].plot()
    #plt.show()
    print 'Done.'
