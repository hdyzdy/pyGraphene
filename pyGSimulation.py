#!/usr/bin/env python

import time
import sys, re
import numpy
import random
import cPickle
from copy import deepcopy 
import matplotlib.pyplot as plt

from pyGNode import pyGNode
from pyGEdge import pyGEdge, pyGCircList
from pyGEnvironment import pyGEnvironment
from pyGNetwork import pyGNetwork
from pyGReaction import pyGReaction, pyGRLoad

class pyGSimulation():
    def __init__(self,t_env=[0],envs=[pyGEnvironment()],network=pyGNetwork(),reactions=pyGRLoad(),name='default_pgs'):
        #random.seed(randSeed)
        self.name = name
        self.t_env = t_env
        self.environments = envs
        self.time = [0]
        self.networks = [network]
        self.nR5 = [network.nR5]
        self.nR6 = [network.nR6]
        self.nC = [network.nC]
        self.siteCounts = [network.countSites()]
        self.reactionsApplied = [(None,None)]
        self.reactions = reactions
        self.reactionRates = []
        self.randState = random.getstate()
        for m in range(len(reactions)):
            self.reactionRates.append(self.reactions[m].rate(self.environments[0]))
        self.reactApplyDict = {}

    def evolve(self,makeCopies=False):
        random.setstate(self.randState)
        #Determine where we are in the specified Time/Species conc history
        #and update environment
        if len(self.time) > 1:
            #print "branch1 self.time",self.time
            for m in range(len(self.t_env)):
                if self.time[-2] >= self.t_env[m]:
                    kEnvOld = m
                    # print "branch1 kEnvOld",kEnvOld
                    # print "branch1 t_env[m]",self.t_env[m],"m=",m
                if self.time[-1] >= self.t_env[m]:
                    # print "hi!"
                    kEnv = m
                    # print "kEnv=",kEnv,"m=",m
                    break
        else:
            # print "branch2 self.time",self.time
            kEnvOld = -1
            for m in range(len(self.t_env)):
                if self.time[-1] >= self.t_env[m]:
                    kEnv = m
                    # print "branch2 kEnv",kEnv
                    break
        #Update reaction rates (only if needed)
        if not kEnv == kEnvOld:
            # print "reactionNumbers=",len(self.reactions)
            for m in range(len(self.reactions)):
                # if self.reactions[m].number == 43 and self.nC[-1] < 50:
                #     continue
                self.reactionRates[m] = self.reactions[m].rate(self.environments[kEnv])
                # print "m=",m,"reactionRate=",self.reactionRates[m]

        for k in range(self.networks[-1].edge.length):
            idStr = self.networks[-1].edge.idStrings[k]
            if '-' in idStr:
                ismolecular = False
                break
            else:
                ismolecular = True

        #Determine tNext for each site
        tstep = numpy.Inf
        for k in range(self.networks[-1].edge.length):
            idStr = self.networks[-1].edge.idStrings[k]
            #Inactive nodes have idStrings == ''
            if idStr:
                #Identify applicable reactions for an unknown site
                if not (idStr in self.reactApplyDict):
                    applicableReactions = []
                    for m in range(len(self.reactions)):
                        # me = self.reactions[m] #reaction index
                        #if self.reactions[m].number == 43 and self.nC[-1] < 50:
                        #     continue
                        # if self.reactions[m].isApplicable(idStr):
                        #        applicableReactions.append(m)

                        if self.reactions[m].number != 43:
                            if self.reactions[m].isApplicable(idStr):
                               applicableReactions.append(m)
                        #        print applicableReactions
                        else:
                            # print "wrong in"
                            if self.reactions[m].isApplicable(idStr) and ismolecular == True: #and self.nC[-1]>90:
                               print "ismolecular!"
                               applicableReactions.append(m)

                    self.reactApplyDict[idStr] = applicableReactions
                    # print "k=",k,"idStr=",idStr,"applicableReactions=",applicableReactions
                netRate = 0
                applicableReactions = self.reactApplyDict[idStr]
                #Calculate timesteps for each site
                for m in range(len(applicableReactions)):
                    netRate += self.reactionRates[applicableReactions[m]]
                if netRate > 0:
                    tmptstep = -numpy.log(random.random())/netRate
                    if tmptstep < tstep:
                        tstep = tmptstep
                        kEdge = k
                        # print "kEdge=",kEdge
                        netRateApply = netRate
                        idStrApply = idStr
        if tstep == numpy.Inf:
            print 'No reactive nodes.  Exiting.'
            return (-1, None)
        applicableReactions = self.reactApplyDict[idStrApply]
        reactionChooser = random.random()
        prob = [self.reactionRates[applicableReactions[0]]/netRateApply]
        if reactionChooser < prob[0]:
            reactionChosen = 0
        else:
            for m in range(1,len(applicableReactions)):
                prob.append(self.reactionRates[applicableReactions[m]]/netRateApply + prob[m-1])
                if reactionChooser < prob[m]:
                    reactionChosen = m
                    break
        if makeCopies: self.networks.append(deepcopy(self.networks[-1])) 
        reactionNumber = self.reactions[applicableReactions[reactionChosen]].number
        reactionName = self.reactions[applicableReactions[reactionChosen]].name
        print "reaction.number=",reactionNumber,"reaction.name=",reactionName
        try:
            kEdge = self.reactions[applicableReactions[reactionChosen]].apply(self.networks[-1],kEdge)
        except Exception, e:
            i_0 = self.networks[-1].edge.nodes[kEdge]
            iN0 = self.networks[-1].nodes[i_0].cns[0]
            iN1 = self.networks[-1].nodes[i_0].cns[1]
            print 'Reaction ' + reactionName + ' failed to complete at edge node ', kEdge, 
            print 'corresponding to network node ', i_0,
            print 'at step ', len(self.time)-1
            print 'node type: ', self.networks[-1].nodes[i_0].type
            print 'node ncn: ', self.networks[-1].nodes[i_0].ncn
            print 'node cns: ', [iN0,iN1]
            print 'cns0 type: ', self.networks[-1].nodes[iN0].type
            print 'cns0 ncn: ', self.networks[-1].nodes[iN0].ncn
            print 'cns1 type: ', self.networks[-1].nodes[iN1].type
            print 'cns1 ncn: ', self.networks[-1].nodes[iN1].ncn
            print 'last 10 reactions: ', self.reactionsApplied[-10:]
            print "Reaction error: ",e.__class__," ",e
            if makeCopies: self.networks.pop()
            return (-1,None)
            #raise
        self.time.append(self.time[-1]+tstep)
        self.reactionsApplied.append((reactionNumber,reactionName))
        self.randState = random.getstate()
        self.nR5.append(self.networks[-1].nR5)
        self.nR6.append(self.networks[-1].nR6)
        # self.nC.append(self.networks[-1].nC)
        # tmporarily modification for the CAHM addition
        # self.nC.append(self.networks[-1].nC+self.networks[-1].nCCAHM)
        self.nC.append(self.networks[-1].nC)
        self.siteCounts.append(self.networks[-1].countSites())
        if re.search(r'R[5,6]B',reactionName): return (2,kEdge)
        elif self.reactions[applicableReactions[reactionChosen]].isStructural: return (1,kEdge)
        else: return (0,kEdge)

    def toXYZ(self,nets='all',filename="pyg.xyz"):
        if (nets=='all'): nets = range(len(self.networks))
        outF=open(filename,'w')
        for k in nets:
#            if (self.reactionsApplied[k]=='None' or eval(self.reactionsApplied[k]+'().isStructural')):
            self.networks[k].toXYZ(filename=filename,title=str(self.time[k]))

    def save(self,filename=None):
        if filename is None:
            filename = self.name + '.p'
        pFile = open(filename,'w')
        cPickle.dump(self,pFile)
        pFile.close()

    def run(self,endMethod='time',endValue=5e-3,outputGeos=False,doOpt=True,debug=False):
        endRun = False
        if (outputGeos): self.networks[-1].toXYZ(title=str(self.time[-1]),filename=self.name+'.xyz')
        if endMethod == 'time': endTest = 'self.time[-1]'
        elif endMethod == 'steps': endTest = 'len(self.time) - 1'
        while (eval(endTest) < endValue and not endRun):
            (evoStatus,kEdge) = self.evolve()
            if debug:
                if not self.networks[-1].nC == len(self.networks[-1].nodes):
                    print 'Inconsistent length and nC after reaction ', self.reactionsApplied[-1]
                    print 'evolve step ', len(self.time)-1
                    endRun = True
                if not self.networks[-1].edge.length == len(self.networks[-1].edge.nodes):
                    print 'Inconsistent edge.length and len(edge.nodes) after reaction ', self.reactionsApplied[-1]
                    print 'evolve step ', len(self.time)-1
                    endRun = True
                for node in self.networks[-1].nodes:
                    if any(numpy.isnan(node.pos)):
                        print 'Node at NaN pos after reaction ', self.reactionsApplied[-1]
                        print 'evolve step ', len(self.time)-1
                        endRun = True
                if endMethod == 'steps' and (len(self.time)-1 > endValue - 10):
                    self.networks[-1].plot()
            if evoStatus < 0:
                if (outputGeos): self.networks[-1].toXYZ(title=str(self.time[-1]),append=True)
                break
            elif evoStatus == 2 or ((self.nR6[-1]-self.nR6[-2]) > 0):
                if doOpt:
                    if (self.networks[-1].isPlanar and re.search(r'^.*R5S.*t',self.reactionsApplied[-1][1])):
                        firstCurve = True
                    else: firstCurve = False 
                    geoStatus = self.networks[-1].checkGeo(optNodes=[self.networks[-1].edge.nodes[kEdge]],
                                                           promoteCurvature=firstCurve,loose=False,noFail=False)
                    if geoStatus < 0:
                        print 'Molecule geometry failed check.'
                        endRun = True
                #Cheating to reduce number of outputs
                if (outputGeos): self.networks[-1].toXYZ(title=str(self.time[-1]),append=True,
                                                                       filename=(self.name+'.xyz'))

    def calcRxnCounts(self,rxns=None,plot=False):
        if rxns is None:
            rxns = []
            for k in range(len(self.reactions)):
                rxns.append(self.reactions[k].number)
        rxns.sort()
        rxnCounts = numpy.zeros(len(rxns))
        for k in range(len(self.reactionsApplied)):
            rxnNumber = self.reactionsApplied[k][0]
            try: #Error is rxnNumber not in rxns (i.e. its a reaction we don't want to count)
                rxnCountIndex = rxns.index(rxnNumber)
                rxnCounts[rxnCountIndex] += 1
            except:
                pass
        if plot:
            left = numpy.arange(0.5,len(rxns),1)
            try:
                plt.bar(left,rxnCounts,width=1.0,log=True,figure=plt.figure())
            except ValueError:
                print('It appears none of the selected reactions occured, try again.')
                plt.close(plt.gcf())
                return rxnCounts
            plt.xticks(left+0.5,map(str,rxns),rotation=90)
        return rxnCounts

if __name__ == '__main__':
    random.seed(0)
    pgn = pyGNetwork(baseType='pc',nRings=7)
    pgs = pyGSimulation(network=pgn,envs=[pyGEnvironment(temp=2500)])
    startTime = time.time()
    pgs.run(endMethod='time',endValue=0.0005,outputGeos=False,doOpt=True,debug=False)
    #pgs.run(endMethod='steps',endValue=119,outputGeos=False,doOpt=True,debug=True)
#    pgs.networks[-1].toTinker(eRMS=0.5,optNodes=pgs.networks[-1].edge.nodes)
    totalTime = time.time() - startTime
    print "Simulation Time: ", pgs.time[-1]
    print "Total steps: ", len(pgs.time)-1
    print 'Total time: ', totalTime/60., ' minutes'
    print 'Last reaction: ', pgs.reactionsApplied[-1]
    #pgs.toXYZ()
    #pgs.networks[-2].plot()
    #pgs.networks[-2].plot()
    pgs.networks[-1].plot()
    sdl = pgs.siteCounts
    print 'nC = ', pgs.nC[-1]
    print pgs.calcRxnCounts(plot=True)
    #pgs.networks[-1].plot()
    if plt.get_fignums(): #If there are figures to plot then plot
        plt.show()
    #print pgs.reactionsApplied
    #print pgs.networks[-1].edge.nodes
    #print pgs.networks[-1].edge.idStrings
    #print pgs.networks[-1].nodes[27].cns
    print 'Done.'

