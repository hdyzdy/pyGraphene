
import re
import os, sys
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import tempfile

from pyGNode import pyGNode
from pyGEdge import pyGEdge, pyGCircList
from pyGEnvironment import pyGEnvironment
from pyGReaction import pyGReaction

class pyGNetwork:
    # tmporarily modification for the CAHM addition
    # nCCAHM = 0

    def __init__(self,cyclic=0,cyclicBoundary=[-float('inf'),float('inf')],nRings=4,baseType='linear'):
        self.cyclic=cyclic
        self.cyclicBoundary=cyclicBoundary
        self.lastOpt = 0
        # tmporarily modification for the CAHM addition
        self.nCCAHM = 0
        if baseType == 'pc':
            self.makePCSubstrate(nRings)
            self.isPlanar = 1
        elif baseType == 'linear':
            self.makeLinearSubstrate(nRings)
            self.isPlanar = 1
            self.classifyEdge()
        else:
            raise 'Invalid pyGNetwork base type.'

    def makePCSubstrate(self,nRings):
        self.nR6 = nRings
        self.nR5 = 0
        if ((nRings == 4) or (nRings == 7)):
            self.makeLinearSubstrate(2,solid=0)
            kEdge = pyGReaction().adsorb6corner(self,3,+1) 
            self.nodes[self.edge.nodes[kEdge+4]].ncn = 'H'
            self.classifyEdge()
            kEdge = pyGReaction().adsorb6boat(self,kEdge-2,+1,'S1','S1')
            self.classifyEdge()
            if (nRings == 7):
                kEdge = pyGReaction().adsorb6corner(self,kEdge+1,+1)
                self.classifyEdge()
                kEdge = pyGReaction().adsorb6boat(self,kEdge+4,+1,'S1','S1')
                self.classifyEdge()
                pyGReaction().adsorb6boat(self,kEdge-3,-1,'S1','S1')
                self.classifyEdge()
        else:
            print "Can't make peri-condensed substrate with ",nRings," rings. Making linear instead."
            self.makeLinearSubstrate(nRings,solid=0)

    def makeLinearSubstrate(self,nRings,solid=1):
        self.nR6 = nRings
        self.nR5 = 0
        posUp = numpy.array([0,0.5*1.42,0])
        mirrorMat = numpy.array([1,-1,1])
        posDown = mirrorMat * posUp
        zig = 1.42*numpy.array([numpy.cos(numpy.pi/6),numpy.sin(numpy.pi/6),0])
        zag = mirrorMat*zig
        if self.cyclic:
            blNodes = [pyGNode(pos=posUp,cns=[1,2,nRings*10-2],ncn='',type='S1')]
            blNodes.append(pyGNode(pos=posDown,cns=[0,3,rings*10-1],ncn='',type='S1'))
        else:
            blNodes = [pyGNode(pos=posUp,cns=[1,2],ncn='H',type='A1')]
            blNodes.append(pyGNode(pos=posDown,cns=[0,3],ncn='H',type='A1'))
        while len(blNodes)+1 < 4 * nRings:
            if numpy.mod(len(blNodes),4) != 0:
                posUp = posUp + zig
                posDown = mirrorMat * posUp
                cnsUp = [len(blNodes)-2,len(blNodes)+2]
                ncnUp = 'H'
                typeUp = 'A1'
                cnsDown = [len(blNodes)-1,len(blNodes)+3]
                if solid:
                    ncnDown = 'solid'
                    typeDown = 'S1'
                else:
                    ncnDown = 'H'
                    typeDown = 'A1'
            else:
                posUp = posUp + zag
                posDown = mirrorMat * posUp
                cnsUp = [len(blNodes)-2,len(blNodes)+1,len(blNodes)+2]
                ncnUp = ''
                typeUp = 'S1'
                cnsDown = [len(blNodes)-1,len(blNodes),len(blNodes)+3]
                ncnDown = ''
                typeDown = 'S1'
            blNodes.append(pyGNode(pos=posUp,cns=cnsUp,ncn=ncnUp,type=typeUp))
            blNodes.append(pyGNode(pos=posDown,cns=cnsDown,ncn=ncnDown,type=typeDown))
        if self.cyclic:
            blNodes[-2].cns = [len(blNodes)-4,0]
            blNodes[-2].type = 'S1'
            blNodes[-1].cns = [len(blNodes)-3,1]
            blNodes[-1].type = 'S1'
            edgeNodes = range(0,len(blNodes)-1,2)
        else:
            posUp = posUp + zag
            posDown = mirrorMat * posUp
            cnsUp = [len(blNodes)-2,len(blNodes)+1]
            ncnUp = 'H'
            typeUp = 'A1'
            cnsDown = [len(blNodes)-1,len(blNodes)]
            ncnDown = 'H'
            typeDown = 'A1'
            blNodes.append(pyGNode(pos=posUp,cns=cnsUp,ncn=ncnUp,type=typeUp))
            blNodes.append(pyGNode(pos=posDown,cns=cnsDown,ncn=ncnDown,type=typeDown))
            edgeNodes = range(0,len(blNodes)-1,2)+range(len(blNodes)-1,-1,-2)
        self.nodes = blNodes
        self.edge = pyGEdge(nodes=edgeNodes,idStrings=len(edgeNodes)*[''],directions=numpy.zeros(len(edgeNodes)))
        self.nC = len(self.nodes)
        self.classifyEdge()

    def plot(self,labels=True):
        plt.figure()
        plotDict = {'A1H':'ko','A2H':'bo','A1-':'ro','A2-':'go',
                    'S1':'ks','S2':'bs','A1C2H':'bD','A2C2H':'kD',
                    'A2H2':'mo','S1solid':'rs','A1C2H3':'yo'}
        xmin = -1
        xmax = 5
        ymin = -1.5
        ymax = 8
        k = 0
        for node in self.nodes:
            try:
                nodePlotStyle = plotDict[node.type+node.ncn]
            except: 
                nodePlotStyle = 'yD'
            plt.plot([node.pos[0]],[node.pos[1]],nodePlotStyle)
            if node.pos[0] < xmin:
                xmin = node.pos[0] - 1
            if node.pos[0] > xmax:
                xmax = node.pos[0] + 1.5
            if node.pos[1] < ymin:
                ymin = node.pos[1] - 1
            if node.pos[1] > ymax:
                ymax = node.pos[1] + 1
            if labels:
                plt.text(node.pos[0],node.pos[1]+0.1,str(k))
                k += 1
        for k in range(self.nC):
            pos1 = self.nodes[k].pos
            for k2 in range(len(self.nodes[k].cns)):
                pos2 = self.nodes[self.nodes[k].cns[k2]].pos
                plt.plot([pos1[0],pos2[0]],[pos1[1],pos2[1]],'k')
        for k in range(self.edge.length):
            nodeA = self.nodes[self.edge.nodes[k]]
            nodeB = self.nodes[self.edge.nodes[k+1]]
            Xa = nodeA.pos[0]; Ya = nodeA.pos[1]
            Xb = nodeB.pos[0]; Yb = nodeB.pos[1]
            plt.plot([Xa,Xb],[Ya,Yb],'--r')

        plt.axis([xmin,xmax,ymin,ymax])

    def plot3d(self,labels=True):
        print "This function still doesn't work..."
        fig = plt.figure()
        ax = Axes3D(fig)
        plotDict = {'A1H':'ko','A2H':'bo','A1-':'ro','A2-':'go',
                    'S1':'ks','S2':'bs','A1C2H':'bD','A2C2H':'kD',
                    'A2H2':'mo','S1solid':'rs'}
        xmin = -1
        xmax = 5
        ymin = -1.5
        ymax = 8
        k = 0
        scatterX = []
        scatterY = []
        scatterZ = []
        for node in self.nodes:
            try:
                nodePlotStyle = plotDict[node.type+node.ncn]
            except: 
                nodePlotStyle = 'yD'
            #ax.scatter([node.pos[0]],[node.pos[1]],[node.pos[2]],nodePlotStyle)
            #ax.scatter([node.pos[0]],[node.pos[1]])
            #print node.pos
            scatterX.append(node.pos[0])
            scatterY.append(node.pos[1])
            scatterZ.append(node.pos[2])
            if node.pos[0] < xmin:
                xmin = node.pos[0] - 1
            if node.pos[0] > xmax:
                xmax = node.pos[0] + 1.5
            if node.pos[1] < ymin:
                ymin = node.pos[1] - 1
            if node.pos[1] > ymax:
                ymax = node.pos[1] + 1
            if labels:
                #ax.text(node.pos[0],node.pos[1]+0.1,node.pos[2],str(k))
                k += 1
        for k in range(self.nC):
            pos1 = self.nodes[k].pos
            for k2 in range(len(self.nodes[k].cns)):
                pos2 = self.nodes[self.nodes[k].cns[k2]].pos
                #ax.plot([pos1[0],pos2[0]],[pos1[1],pos2[1]],[pos1[2],pos2[2]],'k')
                #ax.plot((pos1[0],pos2[0]),(pos1[1],pos2[1]))
        for k in range(self.edge.length):
            nodeA = self.nodes[self.edge.nodes[k]]
            nodeB = self.nodes[self.edge.nodes[k+1]]
            Xa = nodeA.pos[0]; Xb = nodeB.pos[1]
            Ya = nodeA.pos[1]; Yb = nodeB.pos[1]
            Za = nodeA.pos[2]; Zb = nodeB.pos[2]
            #ax.plot([Xa,Xb],[Ya,Yb],[Za,Zb],'--r')
            ax.plot([Xa,Xb],[Ya,Yb],[Za,Zb])
        #ax.axis([xmin,xmax,ymin,ymax])
        scatterX = numpy.array(scatterX)
        scatterY = numpy.array(scatterY)
        scatterZ = numpy.array(scatterZ)
        if any(isnan(scatterX)) or any(isnan(scatterY)) or any(isnan(scatterZ)):
            print "There IS a NaN in the scatter data."
        #ax.scatter(scatterX,scatterY,scatterZ)

    def toXYZ(self,filename='pyg.xyz',title='A Title',colors=False,tinker=False,append=False):
        if append: outF = open(filename,'a')
        else: outF = open(filename,'w')
        outF.write(str(self.nC)+'\n') #Number of atoms in spec
        if not tinker: outF.write(title+'\n')  #Title line (not used in tinker file)
        print '\n\n\n\ni am here!!!!\n'
        self.plot()
        plt.show()
        print '\n\n\n\n\n\nhey\n\n\n\n'
        for k in range(self.nC):
            if tinker:
                continue
                if len(self.nodes[k].cns) == 2:
                    outF.write('{0:d} C {1[0]:.3g} {1[1]:.3g} {1[2]:.3g} 2 {2[0]:d} {2[1]:d}\n'\
                               .format(k+1,self.nodes[k].pos,[x+1 for x in self.nodes[k].cns]))
                elif len(self.nodes[k].cns) == 1:
                    outF.write('{0:d} C {1[0]:.3g} {1[1]:.3g} {1[2]:.3g} 2 {2[0]:d}\n'\
                               .format(k+1,self.nodes[k].pos,[x+1 for x in self.nodes[k].cns]))
                else:                    
                    outF.write('{0:d} C {1[0]:.3g} {1[1]:.3g} {1[2]:.3g} 2 {2[0]:d} {2[1]:d} {2[2]:d}\n'\
                               .format(k+1,self.nodes[k].pos,[x+1 for x in self.nodes[k].cns]))
            else:
                outF.write('C {1[0]:.3g} {1[1]:.3g} {1[2]:.3g}\n'\
                           .format(k+1,self.nodes[k].pos))

    def toTinker(self,tempDir=None,eRMS=0.5,optNodes=None,sphere=6.0):
        tinkerPath = os.getenv('TINKERPATH')
        if tinkerPath is None:
            print "TINKERPATH environment variable not set."
            print "Quitting."
            return -1
        if tempDir is None:
            tempDir = tempfile.mkdtemp()
            rmDir = True
        else:
            rmDir = False
        xyzFileName = os.path.join(tempDir,'tinker_opt.xyz')
        self.toXYZ(tinker=True,filename=xyzFileName)
        keyFileName = os.path.join(tempDir,'tinker_opt.key')
        keyF = open(keyFileName,'w')
        keyF.write('MAXITER 100\n')
        if optNodes is None: optNodes = range(self.nC)
        for iNode in optNodes: 
            keyF.write('SPHERE '+str(iNode+1)+' '+str(sphere)+'\n')
        if self.cyclic:
            cbMag = self.cyclicBoundary(2) - self.cyclicBoundary(1)
            keyF.write('A-axis '+str(cbMag)+'\n')
        keyF.close()
        mm3FileName = os.path.join(tempDir,'tinker_opt_mm3.prm')
        if sys.platform == 'win32':
            print 'here1\n'
            os.system('copy "'+os.path.join(tinkerPath,'mm3.prm')+'" '+mm3FileName+' > NUL')
            logFileName = os.path.join(tempDir,'tinker_opt.log')
            print '"'+os.path.join(tinkerPath,'optimize')+'" '+xyzFileName+' '+\
                mm3FileName[:-4]+' '+str(eRMS)+' > '+logFileName
            tinkerStat = os.system('"'+os.path.join(tinkerPath,'optimize')+'" '+xyzFileName+' '+
                mm3FileName[:-4]+' '+str(eRMS)+' > '+logFileName)
        else:
            print 'hre2\n'
            os.system('cp '+os.path.join(tinkerPath,'mm3.prm')+' '+mm3FileName)
            logFileName = os.path.join(tempDir,'tinker_opt.log')
            tinkerStat = os.system(os.path.join(tinkerPath,'optimize')+' '+xyzFileName+' '+
                           mm3FileName[:-4]+' '+str(eRMS)+' &> '+logFileName)
        if tinkerStat != 0:
            print "TINKER returned non-zero status."
            print "Quitting."
            return -1
        outFileName = os.path.join(tempDir,'tinker_opt.xyz_2')
        optFile = open(outFileName,'r')
        line = optFile.readline()  #Discard number of atoms
        for k in range(self.nC):
            line = optFile.readline()
            line = re.sub('\s*\s', '|', line) #Replace whitespace with separator
            line = line.split('|')            #Split line at separator
            self.nodes[k].pos = numpy.array([float(line[3]), float(line[4]), float(line[5])]) #Update node positions
        optFile.close()
        if sys.platform == 'win32':
                os.system('del '+os.path.join(tempDir,'tinker_opt*'))
        else:
                os.system('rm '+os.path.join(tempDir,'tinker_opt*'))
        if rmDir:
            os.system('rmdir '+tempDir)
        return 0

    def checkGeo(self,nTrys=0,optNodes=[],promoteCurvature=False,loose=False,noFail=False):
        if promoteCurvature:
            interiorS2 = []
            for k in range(self.nC):
                if ((self.nodes[k].type == 'S2') and (not k in self.edge.nodes)):
                    self.nodes[k].pos[-1] += 0.25
                    interiorS2.append(k)
            print "Promoting curvature."
            print "Interior nodes: ", str(interiorS2)
            optNodes.extend(interiorS2)
            self.isPlanar = 0
        if (nTrys == 0):
            #print "Optimizing because of probs, just around node."
            self.toTinker(eRMS=0.5,optNodes=optNodes,sphere=8.0)
            self.lastOpt = 0
        elif (nTrys < 3):
            print "Optimizing because of probs, just around problem nodes."
            self.toTinker(eRMS=0.1,optNodes=optNodes)
            self.lastOpt = 0
        elif (nTrys == 3):
            print "Optimizing because of probs, whole edge."
            self.toTinker(eRMS=0.5,optNodes=self.edge.nodes,sphere=4.0)
            self.lastOpt = 0
        probs = 0
        probNodes = []
        for k in range(self.edge.length):
            i_0 = self.edge.nodes[k]
            if self.nodes[i_0].type[1] == '1':
                for m in range(len(self.nodes[i_0].cns)):
                    i_N = self.nodes[i_0].cns[m]
                    #Only check distances for R6 neighbors
                    if  self.nodes[i_N].type[1] == '1':
                        bondVect = self.nodes[i_0].pos - self.nodes[i_N].pos
                        bondDist = numpy.linalg.norm(bondVect)
                        if (not loose) and ((bondDist > 1.65) or (bondDist < 1.15)):
                            probs += 1
                            probNodes.append(i_0)
                            break
                        elif bondDist > 2.0:
                            probs += 1
                            probNodes.append(i_0)
                            break
        if (probs == 0): 
            return 0
        elif nTrys == 0:
            return self.checkGeo(nTrys=1,optNodes = probNodes,noFail=noFail,loose=loose)
        elif nTrys == 1:
            #print "Doing firstCurve again."
            return self.checkGeo(nTrys=2,optNodes= probNodes,promoteCurvature=True,noFail=noFail,loose=loose)
        elif nTrys == 2:
            return self.checkGeo(nTrys=3,noFail=noFail,loose=loose)
        else:
            if(noFail):
                print 'Geo failed check but in noFail mode so continuing.'
                return 0
            else:
                print "Problem nodes: ", str(probNodes)
                if len(probNodes) == 2:
                    print "Distance between nodes: ", numpy.linalg.norm(self.nodes[probNodes[0]].pos-self.nodes[probNodes[1]].pos)
                return -1

    def classifyEdge(self,kEdge=None):
        self.edge.length = len(self.edge.nodes)
        if (kEdge is None) or (self.edge.str == pyGCircList('')):
            self.edge.str = pyGCircList([])
            self.edge.actives = pyGCircList([])
            self.edge.idStrings = pyGCircList(['']*self.edge.length)
            self.edge.directions = pyGCircList([0]*self.edge.length)
            for k in range(self.edge.length):
                ke = self.edge.nodes[k]  #Edge node index
                self.edge.str.append(self.nodes[ke].type + self.nodes[ke].ncn)
                if self.nodes[ke].type[0] == 'A':
                    self.edge.actives.append(k) 
            reClassNodes = range(len(self.edge.actives))
        else:
            #print 'hi!'
            #print kEdge
            i_0 = self.edge.nodes[kEdge]
            self.edge.str[kEdge] = self.nodes[i_0].type + self.nodes[i_0].ncn
            actIndex = self.edge.actives.index(kEdge)
            reClassNodes = range(actIndex-3,actIndex+4)
            #print self.edge.actives
            #print i_0   
            #print actIndex
            # print reClassNodes
            # self.plot()

            # plt.show()
        if len(self.edge.actives) < 6:
            raise ValueError('Edge to small for classification.')
        for k in reClassNodes:
            activesK = self.edge.actives[k]
            downStr = self.edge.str.circgetslice(self.edge.actives[k-3],activesK)
            # print 'k\t',k,self.edge.actives[k-3]
            # print downStr
            downStr.reverse()
            downStr = ''.join(downStr)
            upStr = ''.join(self.edge.str.circgetslice(activesK+1,self.edge.actives[k+3]+1))
            downStr = re.sub(r'S1solid.*$','S1solid',downStr)
            upStr = re.sub(r'S1solid.*$','S1solid',upStr)
            mStr = self.edge.str[activesK]
            if upStr > downStr:
                self.edge.idStrings[activesK] = mStr + '_' + upStr + '_' + downStr
                self.edge.directions[activesK] = -1
            else:
                self.edge.idStrings[activesK] = mStr + '_' + downStr + '_' + upStr
                self.edge.directions[activesK] = +1

    def removeNode(self,kNode):
        self.nodes.pop(kNode)
        if kNode in self.edge.nodes: self.edge.nodes.remove(kNode)
        for k in range(len(self.nodes)):  
            for m in range(len(self.nodes[k].cns)):
                if self.nodes[k].cns[m] > kNode:
                    self.nodes[k].cns[m] -= 1
        for k in range(len(self.edge.nodes)):
            if self.edge.nodes[k] > kNode:
                self.edge.nodes[k] -= 1

    def countSites(self):
        siteDict = {'ac':0,'ac5':0,'zz':0,'zz5':0,'fe':0,'fe5':0,
                        'bay5':0,'bay5l':0,'bay6':0,'bay6l':0,'other':0}
        for m in range(self.edge.length):
            idStr = self.edge.idStrings[m]
            if idStr:  #Only actives have non-empty ids
                siteID = idStr.split('_')
                upStr = siteID[2 - (self.edge.directions[m]<0)]
                if re.match(r'S1S1A1',upStr): #If true armchair all will be ones
                    siteDict['ac'] += 1
                elif re.match(r'S[1-2]S[1-2]A[1-2]',upStr):  #Otherwise some may be twos
                    siteDict['ac5'] += 1
                elif re.match(r'S1A1',upStr):  #Zigzags have one solid in between.
                    siteDict['zz'] += 1
                elif re.match(r'S[1-2]A[1-2]',upStr): #One solid one active... non-reactive
                    siteDict['zz5'] += 1
                elif re.match(r'A1',upStr):  #A1's can't be next two A2's
                    siteDict['fe'] += 1
                elif re.match(r'A2',upStr):  #A1's can't be next two A2's
                    siteDict['fe5'] += 1
                elif re.match(r'S1S1S1A1',upStr):
                    siteDict['bay5'] += 1
                elif re.match(r'S[1-2]S1S[1-2]A[1-2]',upStr):
                    siteDict['bay5l'] += 1
                elif re.match(r'S1S1S1S1A1',upStr):
                    siteDict['bay6'] += 1
                elif re.match(r'S[1-2]S[1-2]S[1-2]S[1-2]A[1-2]',upStr):
                    siteDict['bay6l'] += 1
                else:
                    siteDict['other'] += 1
        return siteDict

if __name__ == '__main__':
    pyGN = pyGNetwork(baseType='pc',nRings=7) 
    pyGN.toTinker()
    print pyGN.nodes[5]
    #print pyGN.edge.idStrings
    #pyGN.classifyEdge()
#    pyGN.toXYZ()
    print pyGN.edge.nodes
    print pyGN.edge.idStrings
    pyGN.plot()
    plt.savefig('pgne.png')
    print 'Done.'

