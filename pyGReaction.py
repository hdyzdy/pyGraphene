
import re
import numpy
import random

from pyGNode import pyGNode

class pyGReaction:
    def __init__(self):
        pass
    
    def CAHM(self,pyGNe,kEdge,type0):
        if (kEdge < 2): nTurned = 2
        elif (kEdge > (pyGNe.edge.length-2) ): nTurned = -2
        else: nTurned = 0
        pyGNe.edge.turn(nTurned)
        kEdge += nTurned
        i_0 = pyGNe.edge.nodes[kEdge]
        ip1 = pyGNe.edge.nodes[kEdge-1]
        ip2 = pyGNe.edge.nodes[kEdge+1]
        (pos1new, pos2new) = pyGNe.nodes[i_0].setCAHMpos( 
                                pyGNe.nodes[ip1],
                                pyGNe.nodes[ip2])
        pyGNe.nodes.append(pyGNode(pos=pos1new,
                                   cns=[i_0,pyGNe.nC+1],
                                   ncn='H',type='C3'))
        pyGNe.nodes.append(pyGNode(pos=pos2new,
                                   cns=[pyGNe.nC],
                                   ncn='H2',type='C4')) 
        pyGNe.nodes[i_0].type = type0
        pyGNe.nodes[i_0].ncn = 'C2H3' 
        pyGNe.nodes[i_0].cns.append(pyGNe.nC)
        
        pyGNe.nC += 2
        #pyGNe.edge.length += 2
        return kEdge                       
    def adsorb6boat(self,pyGNe,kEdge,side,type0,type3):
        if (kEdge < 4 and side < 0): nTurned = 4
        elif (kEdge > (pyGNe.edge.length-5) and side > 0): nTurned = -4
        else: nTurned = 0
        pyGNe.edge.turn(nTurned)
        kEdge += nTurned
        i_0 = pyGNe.edge.nodes[kEdge]
        ip1 = pyGNe.edge.nodes[kEdge+1*side]
        ip2 = pyGNe.edge.nodes[kEdge+2*side]
        ip3 = pyGNe.edge.nodes[kEdge+3*side]
        #if side > 0:
        #    (i_0,ip1,ip2,ip3) = pyGNe.edge.nodes.circgetslice(kEdge,kEdge+4)
        #else:
        #    (ip3,ip2,ip1,i_0) = pyGNe.edge.nodes.circgetslice(kEdge-3,kEdge+1)
        (pos1new, pos2new) = pyGNe.nodes[i_0].setR6pos( 
                                pyGNe.nodes[ip1],
                                pyGNe.nodes[ip2], 
                                pyGNe.nodes[ip3])
        pyGNe.nodes.append(pyGNode(pos=pos1new,
                                   cns=[i_0,pyGNe.nC+1],
                                   ncn='H',type='A1'))
        pyGNe.nodes.append(pyGNode(pos=pos2new,
                                   cns=[ip3,pyGNe.nC],
                                   ncn='H',type='A1'))
        pyGNe.nodes[i_0].type = type0
        pyGNe.nodes[i_0].ncn = '' 
        pyGNe.nodes[i_0].cns.append(pyGNe.nC)
        
        pyGNe.nodes[ip3].type = type3    
        pyGNe.nodes[ip3].ncn = ''
        pyGNe.nodes[ip3].cns.append(pyGNe.nC+1)
        pyGNe.edge.nodes[kEdge+1*side] = pyGNe.nC
        pyGNe.edge.nodes[kEdge+2*side] = pyGNe.nC+1
        #pyGNe.edge.turn(-nTurned)
        pyGNe.nR6 += 1
        pyGNe.nC += 2
        #pyGNe.edge.length += 2
        return kEdge

    def adsorb5zigzag(self,pyGNe,kEdge,side,ncn1,ncn2):
        if (kEdge < 2 and side < 0): nTurned = 2
        elif (kEdge > (pyGNe.edge.length-3) and side > 0): nTurned = -2
        else: nTurned = 0
        pyGNe.edge.turn(nTurned)
        kEdge += nTurned
        i_0 = pyGNe.edge.nodes[kEdge]
        ip1 = pyGNe.edge.nodes[kEdge+1*side]
        ip2 = pyGNe.edge.nodes[kEdge+2*side]
        (pos1new, pos2new) = pyGNe.nodes[i_0].setR5pos(pyGNe.nodes[ip1],pyGNe.nodes[ip2])
        pyGNe.nodes.append(pyGNode(pos=pos1new,cns=[i_0,pyGNe.nC+1],ncn=ncn1,type='A2'))
        pyGNe.nodes.append(pyGNode(pos=pos2new,cns=[ip2,pyGNe.nC],ncn=ncn2,type='A2'))
        pyGNe.nodes[i_0].type = 'S2'
        pyGNe.nodes[ip2].type = 'S2'
        pyGNe.nodes[i_0].ncn = ''
        pyGNe.nodes[ip2].ncn = ''
        pyGNe.nodes[i_0].cns.append(pyGNe.nC)
        pyGNe.nodes[ip2].cns.append(pyGNe.nC+1)
        if side < 0:
            pyGNe.edge.nodes.pop(kEdge+side)
            pyGNe.edge.nodes.insert(kEdge-1,pyGNe.nC)
            pyGNe.edge.nodes.insert(kEdge-1,pyGNe.nC+1)
        else:
            pyGNe.edge.nodes.pop(kEdge+side)
            pyGNe.edge.nodes.insert(kEdge+1,pyGNe.nC+1)
            pyGNe.edge.nodes.insert(kEdge+1,pyGNe.nC)
        #pyGNe.edge.turn(-nTurned)
        pyGNe.nR5 += 1
        pyGNe.nC += 2
        pyGNe.edge.length += 1
        return kEdge

    def adsorb6corner(self,pyGNe,kEdge,side):
        #Add a 6 member ring at a corner.  Assumes reacting site is the opposite side as radical carbon
        if (kEdge < 2): nTurned = 2
        elif (kEdge > (pyGNe.edge.length-3)): nTurned = -2
        else: nTurned = 0
        pyGNe.edge.turn(nTurned)
        kEdge += nTurned
        im1 = pyGNe.edge.nodes[kEdge-1*side]
        i_0 = pyGNe.edge.nodes[kEdge]
        ip1 = pyGNe.edge.nodes[kEdge+1*side]
        ip2 = pyGNe.edge.nodes[kEdge+2*side]
        pyGNe.nodes[i_0].type = 'S1'
        pyGNe.nodes[ip1].type = 'S1'
        pyGNe.nodes[i_0].ncn = ''
        pyGNe.nodes[ip1].ncn = ''
        pyGNe.nodes[i_0].cns.append(pyGNe.nC)
        pyGNe.nodes[ip1].cns.append(pyGNe.nC+1)
        (pos1new,pos4new) = pyGNe.nodes[i_0].setR6pos(
                                pyGNe.nodes[im1],
                                pyGNe.nodes[ip2],
                                pyGNe.nodes[ip1])
        pyGNe.nodes.append(pyGNode(pos=pos1new,    #Adding lN 
                                   cns=[i_0,pyGNe.nC+2],
                                   ncn='H',type='A1'))
        pyGNe.nodes.append(pyGNode(pos=pos4new,    #Adding lN + 1
                                   cns=[ip1,pyGNe.nC+3],
                                   ncn='-',type='A1'))
        (pos2new,pos3new) = pyGNe.nodes[-2].setR6pos(
                                pyGNe.nodes[i_0],
                                pyGNe.nodes[ip1],
                                pyGNe.nodes[-1])
        pyGNe.nodes.append(pyGNode(pos=pos2new,   #Adding lN + 2
                                   cns=[pyGNe.nC,pyGNe.nC+3],
                                   ncn='H',type='A1'))
        pyGNe.nodes.append(pyGNode(pos=pos3new,   #Adding lN + 3
                                   cns=[pyGNe.nC+1,pyGNe.nC+2],
                                   ncn='H',type='A1'))
        if side < 0:
            pyGNe.edge.nodes.insert(kEdge  ,pyGNe.nC+1)
            pyGNe.edge.nodes.insert(kEdge+1,pyGNe.nC+3)
            pyGNe.edge.nodes.insert(kEdge+2,pyGNe.nC+2)
            pyGNe.edge.nodes.insert(kEdge+3,pyGNe.nC)
        else:
            pyGNe.edge.nodes.insert(kEdge+1,pyGNe.nC+1)
            pyGNe.edge.nodes.insert(kEdge+1,pyGNe.nC+3)
            pyGNe.edge.nodes.insert(kEdge+1,pyGNe.nC+2)
            pyGNe.edge.nodes.insert(kEdge+1,pyGNe.nC)
        #pyGNe.edge.turn(-nTurned)
        pyGNe.nR6 += 1
        pyGNe.nC += 4
        pyGNe.edge.length += 4
        return kEdge

    def desorb6corner(self,pyGNe,kEdge,side,ncn_0,ncnp5):
        if (kEdge < 5): nTurned = 5
        elif (kEdge > (pyGNe.edge.length-6)): nTurned = -5
        else: nTurned = 0
        pyGNe.edge.turn(nTurned)
        kEdge += nTurned
        #im1 = pyGNe.edge.nodes[kEdge-side]
        i_0 = pyGNe.edge.nodes[kEdge]
        ip1 = pyGNe.edge.nodes[kEdge+side]
        ip2 = pyGNe.edge.nodes[kEdge+side*2]
        ip3 = pyGNe.edge.nodes[kEdge+side*3]
        ip4 = pyGNe.edge.nodes[kEdge+side*4]
        ip5 = pyGNe.edge.nodes[kEdge+side*5]
        pyGNe.nodes[i_0].type = 'A1'
        pyGNe.nodes[ip5].type = 'A1'
        pyGNe.nodes[i_0].ncn = ncn_0
        pyGNe.nodes[ip5].ncn = ncnp5
        pyGNe.nodes[i_0].cns.remove(ip1)
        pyGNe.nodes[ip5].cns.remove(ip4)
        ii = sorted([ip1,ip2,ip3,ip4],reverse=True)
        for i in ii:
            pyGNe.removeNode(i)
        #if side < 0:
        #    iPop = kEdge - 3
            #pyGNe.edge.nodes = pyGNe.edge.nodes.circgetslice(kEdge+1,kEdge-3)
        #else:
        #    iPop = kEdge
            #pyGNe.edge.nodes = pyGNe.edge.nodes.circgetslice(kEdge+4,kEdge)
        #for m in range(4):
        #    pyGNe.edge.nodes.pop(iPop)
        #pyGNe.edge.turn(-nTurned)
        pyGNe.nR6 -= 1
        pyGNe.nC -= 4
        pyGNe.edge.length -= 4
        return kEdge

    def desorb6boat(self,pyGNe,kEdge,side,ncnm1,ncnp2):
        if (kEdge < 2): nTurned = 2
        elif (kEdge > (pyGNe.edge.length-3)): nTurned = -2
        else: nTurned = 0
        pyGNe.edge.turn(nTurned)
        kEdge += nTurned
        im1 = pyGNe.edge.nodes[kEdge-side]
        i_0 = pyGNe.edge.nodes[kEdge]
        ip1 = pyGNe.edge.nodes[kEdge+side]
        ip2 = pyGNe.edge.nodes[kEdge+side*2]
        for m in pyGNe.nodes[im1].cns:
            if (not m in pyGNe.edge.nodes):
                pyGNe.edge.nodes[kEdge] = m
                i_0new = m
                break
        for m in pyGNe.nodes[ip2].cns:
            if (not m in pyGNe.edge.nodes):
                pyGNe.edge.nodes[kEdge+side] = m
                ip1new = m
                break
        pyGNe.nodes[im1].ncn = ncnm1
        pyGNe.nodes[ip2].ncn = ncnp2
        pyGNe.nodes[im1].cns.remove(i_0)
        pyGNe.nodes[ip2].cns.remove(ip1)
        pyGNe.nodes[im1].type = re.sub('S','A',pyGNe.nodes[im1].type)
        pyGNe.nodes[ip2].type = re.sub('S','A',pyGNe.nodes[ip2].type)
        if pyGNe.nodes[im1].type[1] == '2':  #This is for five-and-six fall apart
            pyGNe.nodes[i_0new].type = 'S2'
        if pyGNe.nodes[ip2].type[1] == '2':  #This is for five-and-six fall apart
            pyGNe.nodes[ip1new].type = 'S2'
        ii = sorted([i_0,ip1],reverse=True)
        for i in ii:
            pyGNe.removeNode(i)
        #pyGNe.edge.turn(-nTurned)
        pyGNe.nR6 -= 1
        pyGNe.nC -= 2
        #pyGNe.edge.length -= 2
        return kEdge

    def desorb5zigzag(self,pyGNe,kEdge,side,ncnm1,ncnp2):
        if (kEdge < 2): nTurned = 2
        elif (kEdge > (pyGNe.edge.length-3)): nTurned = -2
        else: nTurned = 0
        pyGNe.edge.turn(nTurned)
        kEdge += nTurned
        im1 = pyGNe.edge.nodes[kEdge-side]
        i_0 = pyGNe.edge.nodes[kEdge]
        ip1 = pyGNe.edge.nodes[kEdge+side]
        ip2 = pyGNe.edge.nodes[kEdge+side*2]
        pyGNe.nodes[im1].type= 'A1'
        pyGNe.nodes[ip2].type= 'A1'
        pyGNe.nodes[im1].ncn = ncnm1
        pyGNe.nodes[ip2].ncn = ncnp2
        pyGNe.nodes[im1].cns.remove(i_0)
        pyGNe.nodes[ip2].cns.remove(ip1)
        if pyGNe.nodes[im1].cns[0] in pyGNe.edge.nodes:
            iUncovered = pyGNe.nodes[im1].cns[1]
        else:
            iUncovered = pyGNe.nodes[im1].cns[0]
        pyGNe.nodes[iUncovered].type = 'S1'
        if side < 0:
            pyGNe.edge.nodes.insert(kEdge-1,iUncovered)
        else:
            pyGNe.edge.nodes.insert(kEdge,iUncovered)
        ii = sorted([i_0,ip1],reverse=True)
        for i in ii:
            pyGNe.removeNode(i)
        #pyGNe.edge.turn(-nTurned)
        pyGNe.nR5 -= 1
        pyGNe.nC -= 2
        pyGNe.edge.length -= 1
        return kEdge

    def flipR5R6(self,pyGNe,kEdge,side):
        if (kEdge < 3): nTurned = 3
        elif (kEdge > (pyGNe.edge.length-4)): nTurned = -3
        else: nTurned = 0
        pyGNe.edge.turn(nTurned)
        kEdge += nTurned
        im3 = pyGNe.edge.nodes[kEdge-side*3]
        im2 = pyGNe.edge.nodes[kEdge-side*2]
        im1 = pyGNe.edge.nodes[kEdge-side]
        i_0 = pyGNe.edge.nodes[kEdge]
        ip1 = pyGNe.edge.nodes[kEdge+side]
        ip2 = pyGNe.edge.nodes[kEdge+side*2]
        ip3 = pyGNe.edge.nodes[kEdge+side*3]
        for m in pyGNe.nodes[ip1].cns:
           if not m in pyGNe.edge.nodes:
                i_intNeighb = m
                break
        v1 = pyGNe.nodes[ip1].pos - pyGNe.nodes[i_intNeighb].pos
        v2 = pyGNe.nodes[i_0].pos - pyGNe.nodes[i_intNeighb].pos
        v1norm = numpy.linalg.norm(v1)
        v2norm = numpy.linalg.norm(v2)
        v1 = v1 * v2norm/v1norm
        v2 = v2 * v1norm/v2norm
        pyGNe.nodes[i_0].pos = pyGNe.nodes[i_intNeighb].pos + v2
        pyGNe.nodes[ip1].pos = pyGNe.nodes[i_intNeighb].pos + v1
        pyGNe.nodes[i_0].cns.append(i_intNeighb)
        pyGNe.nodes[ip1].cns.remove(i_intNeighb)
        pyGNe.nodes[i_0].type = 'S2'
        pyGNe.nodes[ip1].type = 'A1'
        pyGNe.nodes[i_0].ncn = ''
        pyGNe.nodes[ip1].ncn = 'H'
        pyGNe.nodes[i_intNeighb].cns.remove(ip1)
        pyGNe.nodes[i_intNeighb].cns.append(i_0)
        if re.match(r'A',pyGNe.nodes[im1].type):
            pyGNe.nodes[im1].type = 'A2'
            if re.match(r'A',pyGNe.nodes[im2].type):
                pyGNe.nodes[im2].type = 'A2'
                pyGNe.nodes[im3].type = 'S2'
            else:
                pyGNe.nodes[im2].type = 'S2'
        else:
            pyGNe.nodes[im1].type = 'S2'
        if re.match(r'A',pyGNe.nodes[ip2].type):
            pyGNe.nodes[ip2].type = 'A1'
            pyGNe.nodes[ip3].type = 'S1'
        else:
            pyGNe.nodes[ip2].type = 'S1'
        #pyGNe.edge.turn(-nTurned)
        return kEdge

class pyGR_CstCsm(pyGReaction):
    #Dissertation Reaction 1
    isStructural = False
    number = 1
    name = 'CstCsm'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        i_0 = pyGNe.edge.nodes[kEdge]
        pyGNe.nodes[i_0].ncn = '-'
        pyGNe.classifyEdge(kEdge=kEdge)
        return kEdge
    def rate(self, pyGEn):    #Keifer et al. per site based on benzene
        T = pyGEn.temperature #J. Phys. Chem. 89 2013-2019 (1985)
        cOfH = pyGEn.concOf('H')
        return 4.2e13 * numpy.exp(-8052./T) * cOfH
    def isApplicable(self, siteType):
        if re.match(r'A1H',siteType): return True
        return False

class pyGR_CsmtCs(pyGReaction):
    #Dissertation Reaction 2
    isStructural = False
    number = 2
    name = 'CsmtCs'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        kNode = pyGNe.edge.nodes[kEdge]
        pyGNe.nodes[kNode].ncn = 'H'
        pyGNe.classifyEdge(kEdge=kEdge)
        return kEdge
    def rate(self,pyGEn):     #Reverse of abstraction plus H addition
        T = pyGEn.temperature #Keq for abstraction from B3LYP/6-311G(d,p) for benzene
        cOfH2 = pyGEn.concOf('H2')
        cOfH  = pyGEn.concOf('H')
        Keq_abstr = 7.58976 * numpy.exp(-2097.3/T) #Least squares fit of ln(K) = ln(A) + B/T
        coeff_H2 = pyGR_CstCsm().rate(pyGEn) / Keq_abstr / cOfH
        coeff_H = 2.0e13 #Frenklach et al. 30th Combustion Symposium
        return coeff_H*cOfH + coeff_H2*cOfH2
    def isApplicable(self, siteType):
        flippable = (pyGR_CsR5R6mtCsR6mR5E().isApplicable(siteType) or
                     pyGR_CsR5R6mtCsR6mR5I().isApplicable(siteType) or
                     pyGR_CsR5R6mtCsR6mR5EI().isApplicable(siteType) or
                     pyGR_CsR5R6mtCsR6mR5IE().isApplicable(siteType) )
        if (re.match(r'A1-',siteType) and (not flippable)): return True
        return False

class pyGR_CsmtCsR5H(pyGReaction):
    #Dissertation Reaction 3
    isStructural = True
    number = 3
    name = 'CsmtCsR5H'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if (re.match(r'S1A1H',siteID[1]) and re.match(r'S1A1H',siteID[2])): #Both sides are boat
            # print siteID
            if random.random() > 0.5: #Go up
                side = +1
            else:                        #Go down
                side = -1
        elif re.match(r'S1A1H',siteID[1]): #First side is boat 
            # print siteID
            side = siteDir * -1        
        elif re.match(r'S1A1H',siteID[2]): #Second side is boat
            # print siteID
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for adsorbing 5-member ring.')
        kEdge = self.adsorb5zigzag(pyGNe,kEdge,side,'H','H2')
        pyGNe.classifyEdge()
        # print pyGNe.edge.idStrings[kEdge]
        return kEdge
    def rate(self,pyGEn):
        T = pyGEn.temperature
        cOfC2H2 = pyGEn.concOf('C2H2')
        #Assumes steady state for radical intermediates with C2H2 addition from 
        #26th and 5ring cyclization from 27th.
        kap = 1.1e7 * T**1.71*numpy.exp(-16300.0/(8.314*T))  # Rates from Frenklach 26th for 
        kmap = 1.3e14 * numpy.exp(-174800.0/(8.314*T))       # reactions of Cs with C2H2
        kbp = 4.8e12*numpy.exp(-140300.0/(8.314*T))
        kcp = 6.8e11*numpy.exp(-11084.0/T)
        return cOfC2H2 * kap*kcp/(kmap + kbp + kcp)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1-',siteID[0]) and 
            (re.match(r'S1A1H',siteID[1]) or re.match(r'S1A1H',siteID[2]))):
            return True
        return False 

class pyGR_CsmtCsC2H(pyGReaction):
    #Dissertation Reaction 4
    isStructural = False
    number = 4
    name = 'CsmtCsC2H'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        kNode = pyGNe.edge.nodes[kEdge]
        pyGNe.nodes[kNode].ncn = 'C2H'
        pyGNe.classifyEdge(kEdge=kEdge)
        return kEdge
    def rate(self,pyGEn):     
        T = pyGEn.temperature 
        cOfC2H2 = pyGEn.concOf('C2H2')
        #Rates assuming steady state equation from Frenklach 26th Symposium
        kap = 1.1e7 * T**1.71 * numpy.exp(-1960./T)
        kmap = 1.3e14 * numpy.exp(-21024.8/T)
        kbp = 4.8e12 * numpy.exp(-16875./T)
        #kmbp = 1.5e10*T**0.85*numpy.exp(-601.4/T);
        kcp = 6.8e11*numpy.exp(-11084.0/T)
        coeff = kap*kbp/(kmap+kbp+kcp) 
        return coeff*cOfC2H2 
    def isApplicable(self, siteType):
        if re.match(r'A1-',siteType): return True
        return False

class pyGR_CsC2HtCsm(pyGReaction):
    #Dissertation Reaction 5
    isStructural = False
    number = 5
    name = 'CsC2HtCsm'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        kNode = pyGNe.edge.nodes[kEdge]
        pyGNe.nodes[kNode].ncn = '-'
        pyGNe.classifyEdge(kEdge=kEdge)
        return kEdge
    def rate(self,pyGEn):     
        T = pyGEn.temperature 
        cOfH = pyGEn.concOf('H')
        #Rates assuming steady state equation from Frenklach 26th Symposium
        #kap = 1.1e7 * T**1.71 * numpy.exp(-1960./T)
        kmap = 1.3e14 * numpy.exp(-21024.8/T)
        kbp = 4.8e12 * numpy.exp(-16875./T)
        kmbp = 1.5e10*T**0.85*numpy.exp(-601.4/T)
        kcp = 6.8e11*numpy.exp(-11084.0/T)
        coeff = kmap*kmbp/(kmap+kbp+kcp) 
        return coeff*cOfH
    def isApplicable(self, siteType):
        if re.match(r'A1C2H',siteType): return True
        return False

class pyGR_CsC2HtCsR5H(pyGReaction):
    #Dissertation Reaction 6
    isStructural = True
    number = 6
    name = 'CsC2HtCsR5H'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if (re.match(r'S1A1H',siteID[1]) and re.match(r'S1A1H',siteID[2])): #Both sides are boat
            if random.random() > 0.5: #Go up
                side = +1
            else:                        #Go down
                side = -1
        elif re.match(r'S1A1H',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S1A1H',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for adsorbing 5-member ring.')
        kEdge = self.adsorb5zigzag(pyGNe,kEdge,side,'H','H2')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        T = pyGEn.temperature
        cOfH = pyGEn.concOf('H')
        #Assumes steady state for radical intermediates with C2H2 addition from 
        #26th and 5ring cyclization from 27th.
        #kap = 1.1e7 * T**1.71*numpy.exp(-16300.0/(8.314*T))  # Rates from Frenklach 26th for 
        kmap = 1.3e14 * numpy.exp(-174800.0/(8.314*T))       # reactions of Cs with C2H2
        kmbp = 1.5e10 * T**0.85 * numpy.exp(-5000./(8.314*T)) 
        kbp = 4.8e12*numpy.exp(-140300.0/(8.314*T))
        kcp = 6.8e11*numpy.exp(-11084.0/T)
        return cOfH * kmbp*kcp/(kmap + kbp + kcp)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1C2H',siteID[0]) and 
            (re.match(r'S1A1H',siteID[1]) or re.match(r'S1A1H',siteID[2]))):
            return True
        return False 

class pyGR_CsR5HtCsm(pyGReaction):
    #Dissertation Reaction 7
    isStructural = True
    number = 7
    name = 'CsR5HtCsm'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]):
            side = siteDir * 1
        elif re.match(r'S2',siteID[2]):
            side = siteDir * -1
        else:
            raise ValueError('Invalid type for desorbing 5-member zigzag ring.')
        kEdge = self.desorb5zigzag(pyGNe,kEdge,side,'-','H')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        T = pyGEn.temperature
        rateDict = {1500:2.39E+03, 2000:2.47E+06, 2300:2.40E+07, 2500:7.15E+07}
        if T in rateDict:
            rate = rateDict[T]
        else:
            #Interpolate (See JPC KMC Paper).
            rate = 3.1e11*T**0.87*exp(-37403.0/T)
        if T < 1500 or T > 2500:
            print 'Using Arhennius fit outside fitted temperature range.'
        return rate
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A2H(?!2)',siteID[0]) and
           (re.match(r'A2H2',siteID[1]) or re.match(r'A2H2',siteID[2]))):
            return True
        return False

class pyGR_CsR5HtCsC2H(pyGReaction):
    #Dissertation Reaction 8
    isStructural = True
    number = 8
    name = 'CsR5HtCsC2H'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]):
            side = siteDir * 1
        elif re.match(r'S2',siteID[2]):
            side = siteDir * -1
        else:
            raise ValueError('Invalid type for desorbing 5-member zigzag ring.')
        kEdge = self.desorb5zigzag(pyGNe,kEdge,side,'C2H','H')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        T = pyGEn.temperature
        rateDict = {1500:1.35E+04, 2000:9.26E+06, 2300:8.52E+07, 2500:2.56E+08}
        if T in rateDict:
            rate = rateDict[T]
        else:
            #Interpolate (See JPC KMC Paper).
            rate = 6.7e11*T**0.84*numpy.exp(-35625./T)
        if T < 1500 or T > 2500:
            print 'Using Arhennius fit outside fitted temperature range.'
        return rate
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A2H(?!2)',siteID[0]) and
           (re.match(r'A2H2',siteID[1]) or re.match(r'A2H2',siteID[2]))):
            return True
        return False

class pyGR_CsR5HtR5HCs(pyGReaction):
    #Dissertation Reaction 9
    isStructural = True
    number = 9
    name = 'CsR5HtR5HCs'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]):
            side = siteDir * -1
        elif re.match(r'S2',siteID[2]):
            side = siteDir * 1
        else:
            raise ValueError('Invalid type for migrating 5-member ring.')
        kEdge = self.desorb5zigzag(pyGNe,kEdge,-side,'-','H')
        pyGNe.edge.length = len(pyGNe.edge.nodes) #Have to do this to keep consistent
        if side < 0:
            kEdge -= 1
        kEdge = self.adsorb5zigzag(pyGNe,kEdge,side,'H','H2')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        T = pyGEn.temperature
        rateDict = {1500:8.01E+04, 2000:5.35E+06, 2300:2.00E+07, 2500:3.87E+07}
        if T in rateDict:
            rate = rateDict[T]
        else:
            #Interpolate (See JPC KMC Paper).
            rate = 1.3e11*T**0.16*numpy.exp(-23099./T)
        if T < 1500 or T > 2500:
            print 'Using Arhennius fit outside fitted temperature range.'
        return rate
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A2H(?!2)',siteID[0]) and
           ((re.match(r'A2H2',siteID[1]) and re.match(r'S2S1A1H',siteID[2])) or
            (re.match(r'A2H2',siteID[2]) and re.match(r'S2S1A1H',siteID[1])))):
            return True
        return False

class pyGR_CsR5tCsR5m(pyGReaction):
    #Dissertation Reaction 10
    isStructural = False
    number = 10
    name = 'CsR5tCsR5m'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        kNode = pyGNe.edge.nodes[kEdge]
        pyGNe.nodes[kNode].ncn = '-'
        pyGNe.classifyEdge(kEdge=kEdge)
        return kEdge
    def rate(self,pyGEn):    
        T = pyGEn.temperature 
        cOfH = pyGEn.concOf('H')
        #Rate taken from Knyazev et al. JPC 100 11346-11354 (1996)
        coeff = 8.42e-17*T**1.93*numpy.exp(-6518./T)*6.023e23
        return coeff*cOfH
    def isApplicable(self, siteType):
        siteID = siteType.split('_')
        if (re.match(r'A2H(?!2)',siteID[0]) and
           ((re.match(r'A2H(?!2)',siteID[1]) or re.match(r'A2H(?!2)',siteID[2])) or
            (re.match(r'S2',siteID[1]) and re.match(r'S2',siteID[2]) ))):
            return True
        return False

class pyGR_CsR5mtCsR5(pyGReaction):
    #Dissertation Reaction 11
    isStructural = False
    number = 11
    name = 'CsR5mtCsR5'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        kNode = pyGNe.edge.nodes[kEdge]
        pyGNe.nodes[kNode].ncn = 'H'
        pyGNe.classifyEdge(kEdge=kEdge)
        return kEdge
    def rate(self,pyGEn):    
        T = pyGEn.temperature 
        cOfH = pyGEn.concOf('H')
        cOfH2 = pyGEn.concOf('H2')
        #Reverse of abstraction taken from Knyazev et al. JPC 100 11346-11354 (1996)
        #H addition to radical taken from GRI-Mech for ethylene
        coeff_H = 6.08e12*T**0.27*numpy.exp(-141./T)
        coeff_H2 = 1.57e-20*T**2.56*numpy.exp(-2529./T)*6.023e23
        return coeff_H*cOfH + coeff_H2*cOfH2
    def isApplicable(self, siteType):
        if re.match(r'A2-',siteType): return True
        return False

class pyGR_CsR5tCsR5H(pyGReaction):
    #Dissertation Reaction 12
    isStructural = False
    number = 12
    name = 'CsR5tCsR5H'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        kNode = pyGNe.edge.nodes[kEdge]
        pyGNe.nodes[kNode].ncn = 'H2'
        pyGNe.classifyEdge(kEdge=kEdge)
        return kEdge
    def rate(self,pyGEn):    
        T = pyGEn.temperature 
        cOfH = pyGEn.concOf('H')
        #From GRI-Mech for ethylene
        coeff = 5.40E+11*T**0.45*numpy.exp(-916./T)
        return coeff*cOfH
    def isApplicable(self, siteType):
        siteID = siteType.split('_')
        if (re.match(r'A2H(?!2)',siteID[0]) and
           ((re.match(r'A2H(?!2)',siteID[1]) or re.match(r'A2H(?!2)',siteID[2])) or
            (re.match(r'S2',siteID[1]) and re.match(r'S2',siteID[2]) ))):
            return True
        return False

class pyGR_CsR5HtCsR5(pyGReaction):
    #Dissertation Reaction 13
    isStructural = False
    number = 13
    name = 'CsR5HtCsR5'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        kNode = pyGNe.edge.nodes[kEdge]
        pyGNe.nodes[kNode].ncn = 'H'
        pyGNe.classifyEdge(kEdge=kEdge)
        return kEdge
    def rate(self,pyGEn):    
        T = pyGEn.temperature 
        cOfH = pyGEn.concOf('H')
        #H desorption taken as reverse of addition from GRI-Mech for ethylene
        #Using GRI-Mech Keq
        coeff_addition = 5.40E+11*T**0.45*numpy.exp(-916./T)
        rate_addition = coeff_addition*cOfH
                      #1500 K,   2000 K  , 2500 K 
        # KeqArray = [2.41E+05, 1.24E+04, 2.15E+03];  %Taking addition as forward reaction
        # 
        # switch T
        #     case 1500
        #         Keq = KeqArray(1);
        #     case 2000
        #         Keq = KeqArray(2);
        #     case 2500
        #         Keq = KeqArray(3);
        #     otherwise
        #         error('No data for reaction at this temp');
        # end
        #Exponential fit to GRI-Mech 3.0 data for T = 1500 - 2500;
        Keq = 1.791*numpy.exp(17708./T)
        #Also include abstraction by analogy to C2H5 + H -> C2H4 + H2 from GRI-Mech
        rate_abstraction = 2.0e12*cOfH
        return rate_addition/Keq + rate_abstraction
    def isApplicable(self, siteType):
        if re.match(r'A2H2',siteType): return True
        return False

class pyGR_CsR5HpR5tCsR6R5H(pyGReaction):
    #Dissertation Reaction 14
    isStructural = True
    number = 14
    name = 'CsR5HpR5tCsR6R5H'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]):
            side = siteDir * -1
        elif re.match(r'S2',siteID[2]):
            side = siteDir * 1
        else:
            raise ValueError('Invalid type for migrating 5-member ring.')
        addHsiteIndex = pyGNe.edge.nodes[kEdge + side*5]
        pyGNe.nodes[addHsiteIndex].ncn = 'H2'
        kEdge = self.desorb5zigzag(pyGNe,kEdge,-side,'-','H')
        if side < 0:
            kEdge -= 1
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S1','S2')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        T = pyGEn.temperature
        #Rate derived from high pressure limit rate constants from RRKM of JPC
        #flip paper, by direct integration treating 17, 18, and 21 as products.
        rateDict = {1500:2.41E+04, 1750:4.48E+05, 2000:3.78E+06, 2250:1.81E+07, 2500:5.73E+07}
        if T in rateDict:
            rate = rateDict[T]
        else:
            #Interpolate (See JPC KMC Paper).
            rate = 8.9e5*T**2.28*numpy.exp(-30944./T);
        if T < 1500 or T > 2500:
            print 'Using Arhennius fit outside fitted temperature range.'
        return rate
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A2H(?!2)',siteID[0]) and
           ((re.match(r'A2H2',siteID[1]) and re.match(r'S2S1S2A2HA2H(?!2)',siteID[2])) or
            (re.match(r'A2H2',siteID[2]) and re.match(r'S2S1S2A2HA2H(?!2)',siteID[1])))):
            return True
        return False

class pyGR_CsR6R5HtCsR5HpR5(pyGReaction):
    #Dissertation Reaction 15
    isStructural = True
    number = 15
    name = 'CsR6R5HtCsR5HpR5'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]):
            side = siteDir * 1
        elif re.match(r'S2',siteID[2]):
            side = siteDir * -1
        else:
            raise ValueError('Invalid type for migrating 5-member ring.')
        removeHsiteIndex = pyGNe.edge.nodes[kEdge - side*2]
        pyGNe.nodes[removeHsiteIndex].ncn = 'H'
        kEdge = self.desorb6boat(pyGNe,kEdge,side,'H','-')
        kEdge += side*2
        kEdge = self.adsorb5zigzag(pyGNe,kEdge,side,'H','H2')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        T = pyGEn.temperature
        #Rate derived from high pressure limit rate constants from RRKM of JPC
        #flip paper, by direct integration treating 17, 18, and 21 as products.
        rateDict = {1500:6.24E+00,1750:4.62E+02,2000:1.10E+04,2250:1.23E+05,2500:7.59E+05}
        if T in rateDict:
            rate = rateDict[T]
        else:
            #Interpolate (See JPC KMC Paper).
            rate = 2.1e9*T**1.14*numpy.exp(-41952./T);
        if T < 1500 or T > 2500:
            print 'Using Arhennius fit outside fitted temperature range.'
        return rate
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1H',siteID[0]) and
           ((re.match(r'S2A2H2S2',siteID[1]) and re.match(r'A1HS1S1A1H',siteID[2])) or
            (re.match(r'S2A2H2S2',siteID[2]) and re.match(r'A1HS1S1A1H',siteID[1])))):
            return True
        return False

class pyGR_CsR6R5HtCsC2HpR5(pyGReaction):
    #Dissertation Reaction 16
    isStructural = True
    number = 16
    name = 'CsR6R5HtCsC2HpR5'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]):
            side = siteDir * 1
        elif re.match(r'S2',siteID[2]):
            side = siteDir * -1
        else:
            raise ValueError('Invalid type for migrating 5-member ring.')
        removeHsiteIndex = pyGNe.edge.nodes[kEdge - side*2]
        pyGNe.nodes[removeHsiteIndex].ncn = 'H'
        kEdge = self.desorb6boat(pyGNe,kEdge,side,'H','C2H')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        T = pyGEn.temperature
        #Rate derived from high pressure limit rate constants from RRKM of JPC
        #flip paper, by direct integration treating 17, 18, and 21 as products.
        rateDict = {1500:4.20E-01,1750:8.37E+01,2000:4.14E+03,2250:8.17E+04,2500:8.25E+05}
        if T in rateDict:
            rate = rateDict[T]
        else:
            #Interpolate (See JPC KMC Paper).
            rate = 3.8e10*T**1.30*numpy.exp(-51929./T)
        if T < 1500 or T > 2500:
            print 'Using Arhennius fit outside fitted temperature range.'
        return rate
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1H',siteID[0]) and
           ((re.match(r'S2A2H2S2',siteID[1]) and re.match(r'A1HS1(?!solid)',siteID[2])) or
            (re.match(r'S2A2H2S2',siteID[2]) and re.match(r'A1HS1(?!solid)',siteID[1])))):
            return True
        return False

class pyGR_CsR6R5HtCsmpR5(pyGReaction):
    #Dissertation Reaction 17
    isStructural = True
    number = 17
    name = 'CsR6R5HtCsmpR5'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]):
            side = siteDir * 1
        elif re.match(r'S2',siteID[2]):
            side = siteDir * -1
        else:
            raise ValueError('Invalid type for migrating 5-member ring.')
        removeHsiteIndex = pyGNe.edge.nodes[kEdge - side*2]
        pyGNe.nodes[removeHsiteIndex].ncn = 'H'
        kEdge = self.desorb6boat(pyGNe,kEdge,side,'H','-')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        T = pyGEn.temperature
        #Rate derived from high pressure limit rate constants from RRKM of JPC
        #flip paper, by direct integration treating 17, 18, and 21 as products.
        rateDict = {1500:7.54E-02,1750:2.50E+01,2000:1.89E+03,2250:5.21E+04,2500:6.75E+05}
        if T in rateDict:
            rate = rateDict[T]
        else:
            #Interpolate (See JPC KMC Paper).
            rate = 4.0e10*T**1.53*numpy.exp(-57225./T)
        if T < 1500 or T > 2500:
            print 'Using Arhennius fit outside fitted temperature range.'
        return rate
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1H',siteID[0]) and
           ((re.match(r'S2A2H2S2',siteID[1]) and re.match(r'A1HS1(?!solid)',siteID[2])) or
            (re.match(r'S2A2H2S2',siteID[2]) and re.match(r'A1HS1(?!solid)',siteID[1])))):
            return True
        return False

class pyGR_CsR6R5HpR5tCsR5HR6pR5(pyGReaction):
    #Dissertation Reaction 18
    isStructural = True
    number = 18
    name = 'CsR6R5HpR5tCsR5HR6pR5'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]):
            side = siteDir * 1
        elif re.match(r'S2',siteID[2]):
            side = siteDir * -1
        else:
            raise ValueError('Invalid type for migrating 5-member ring.')
        removeHsiteIndex = pyGNe.edge.nodes[kEdge - side*2]
        pyGNe.nodes[removeHsiteIndex].ncn = 'H'
        kEdge = self.desorb6boat(pyGNe,kEdge,side,'H','-')
        kEdge += side*2
        addHsiteIndex = pyGNe.edge.nodes[kEdge+side*4]
        pyGNe.nodes[addHsiteIndex].ncn = 'H2'
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S1','S2')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsR6R5HtCsR5HpR5().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1H',siteID[0]) and
           ((re.match(r'S2A2H2S2',siteID[1]) and re.match(r'A1HS1S1S2A2HA2H(?!2)',siteID[2])) or
            (re.match(r'S2A2H2S2',siteID[2]) and re.match(r'A1HS1S1S2A2HA2H(?!2)',siteID[1])))):
            return True
        return False

class pyGR_CsR6R5HpR5R6tCsR5SR6pR5(pyGReaction):
    #Dissertation Reaction 19
    isStructural = True
    number = 19
    name = 'CsR6R5HpR5R6tCsR5SR6pR5'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]):
            side = siteDir * 1
        elif re.match(r'S2',siteID[2]):
            side = siteDir * -1
        else:
            raise ValueError('Invalid type for migrating 5-member ring.')
        removeHsiteIndex = pyGNe.edge.nodes[kEdge - side*2]
        pyGNe.nodes[removeHsiteIndex].ncn = 'H'
        kEdge = self.desorb6boat(pyGNe,kEdge,side,'H','-')
        kEdge += side*2
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S1','S2')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsR6R5HtCsR5HpR5().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1H',siteID[0]) and
           ((re.match(r'S2A2H2S2',siteID[1]) and re.match(r'A1HS1S1S2A2HS2',siteID[2])) or
            (re.match(r'S2A2H2S2',siteID[2]) and re.match(r'A1HS1S1S2A2HS2',siteID[1])))):
            return True
        return False

class pyGR_CsR6R5HpR6tCsR6R6pR5(pyGReaction):
    #Dissertation Reaction 20
    isStructural = True
    number = 20
    name = 'CsR6R5HpR6tCsR6R6pR5'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]):
            side = siteDir * 1
        elif re.match(r'S2',siteID[2]):
            side = siteDir * -1
        else:
            raise ValueError('Invalid type for migrating 5-member ring.')
        removeHsiteIndex = pyGNe.edge.nodes[kEdge - side*2]
        pyGNe.nodes[removeHsiteIndex].ncn = 'H'
        kEdge = self.desorb6boat(pyGNe,kEdge,side,'H','-')
        kEdge += side*2
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S1','S1')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsR6R5HtCsR5HpR5().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1H',siteID[0]) and
           ((re.match(r'S2A2H2S2',siteID[1]) and re.match(r'A1HS1S1S1A1H',siteID[2])) or
            (re.match(r'S2A2H2S2',siteID[2]) and re.match(r'A1HS1S1S1A1H',siteID[1])))):
            return True
        return False

class pyGR_CsR5HpC2HtCsR6m(pyGReaction):
    #Dissertation Reaction 21
    isStructural = True
    number = 21
    name = 'CsR5HpC2HtCsR6m'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'A2',siteID[1]):
            side = siteDir * 1
        elif re.match(r'A2',siteID[2]):
            side = siteDir * -1
        else:
            raise ValueError('Invalid type for migrating 5-member ring.')
        kEdge = self.desorb5zigzag(pyGNe,kEdge,-side,'H','H')
        if side < 0:
            kEdge -= 1
        kEdge = self.adsorb6corner(pyGNe,kEdge,side)
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsR5HtR5HCs().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A2H(?!2)',siteID[0]) and
           ((re.match(r'S2A1C2H',siteID[1]) and re.match(r'A2H2S2',siteID[2])) or
            (re.match(r'S2A1C2H',siteID[2]) and re.match(r'A2H2S2',siteID[1])))):
            return True
        return False

class pyGR_CsR6R5HpC2HtCsR6mpR5(pyGReaction):
    #Dissertation Reaction 22
    isStructural = True
    number = 22
    name = 'CsR6R5HpC2HtCsR6mpR5'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]):
            side = siteDir * 1
        elif re.match(r'S2',siteID[2]):
            side = siteDir * -1
        else:
            raise ValueError('Invalid type for migrating 5-member ring.')
        removeHsiteIndex = pyGNe.edge.nodes[kEdge - side*2]
        pyGNe.nodes[removeHsiteIndex].ncn = 'H'
        kEdge = self.desorb6boat(pyGNe,kEdge,side,'H','-')
        kEdge += side*2
        kEdge = self.adsorb6corner(pyGNe,kEdge,side)
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsR6R5HtCsR5HpR5().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1H',siteID[0]) and
           ((re.match(r'S2A2H2S2',siteID[1]) and re.match(r'A1HS1A1C2H',siteID[2])) or
            (re.match(r'S2A2H2S2',siteID[2]) and re.match(r'A1HS1A1C2H',siteID[1])))):
            return True
        return False

class pyGR_CsR6mtCsR5HpC2H(pyGReaction):
    #Dissertation Reaction 23
    isStructural = True
    number = 23
    name = 'CsR6mtCsR5HpC2H'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'A1HA1-',siteID[1]):
            side = siteDir * -1
        elif re.match(r'A1HA1-',siteID[2]):
            side = siteDir * 1
        else:
            raise ValueError('Invalid type for removing 6-member corner ring.')
        kEdge -= side*2
        kEdge = self.desorb6corner(pyGNe,kEdge,side,'-','C2H')
        if side < 0:
            kEdge -= 4
        kEdge = self.adsorb5zigzag(pyGNe,kEdge,-side,'H','H2')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        #Taken from 31st R17 (reaction (k)) multiplied by a branching ratio
        #from migration with simultaneous desorption from integrated equations.
        T = pyGEn.temperature 
        coeff = 1.3e11*T**1.08*numpy.exp(-35428./T)
        k7 = pyGR_CsR5HtCsm().rate(pyGEn)
        k8 = pyGR_CsR5HtCsC2H().rate(pyGEn)
        k9 = pyGR_CsR5HtR5HCs().rate(pyGEn)
        br_ratio = k9/(k7+k8+k9)
        if T < 1500 or T > 2500:
            print 'Using Arhennius fit outside fitted temperature range.'
        return coeff*br_ratio
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1H',siteID[0]) and
           ((re.match(r'A1HA1-S1',siteID[1]) and re.match(r'A1HS1S1A1H',siteID[2])) or
            (re.match(r'A1HA1-S1',siteID[2]) and re.match(r'A1HS1S1A1H',siteID[1])))):
            return True
        return False

class pyGR_CsR6mtCsC2HpC2H(pyGReaction):
    #Dissertation Reaction 24
    isStructural = True
    number = 24
    name = 'CsR6mtCsC2HpC2H'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'A1HA1-',siteID[1]):
            side = siteDir * -1
        elif re.match(r'A1HA1-',siteID[2]):
            side = siteDir * 1
        else:
            raise ValueError('Invalid type for removing 6-member corner ring.')
        kEdge -= side*2
        kEdge = self.desorb6corner(pyGNe,kEdge,side,'C2H','C2H')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        #Taken from 31st R17 (reaction (k)) multiplied by a branching ratio
        #from migration with simultaneous desorption from integrated equations.
        T = pyGEn.temperature 
        coeff = 1.3e11*T**1.08*numpy.exp(-35428./T)
        k7 = pyGR_CsR5HtCsm().rate(pyGEn)
        k8 = pyGR_CsR5HtCsC2H().rate(pyGEn)
        k9 = pyGR_CsR5HtR5HCs().rate(pyGEn)
        br_ratio = k8/(k7+k8+k9)
        if T < 1500 or T > 2500:
            print 'Using Arhennius fit outside fitted temperature range.'
        return coeff*br_ratio
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1H',siteID[0]) and
           ((re.match(r'A1HA1-S1',siteID[1]) and re.match(r'A1HS1',siteID[2])) or
            (re.match(r'A1HA1-S1',siteID[2]) and re.match(r'A1HS1',siteID[1])))):
            return True
        return False

class pyGR_CsR6mtCsmpC2H(pyGReaction):
    #Dissertation Reaction 25
    isStructural = True
    number = 25
    name = 'CsR6mtCsmpC2H'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'A1HA1-',siteID[1]):
            side = siteDir * -1
        elif re.match(r'A1HA1-',siteID[2]):
            side = siteDir * 1
        else:
            raise ValueError('Invalid type for removing 6-member corner ring.')
        kEdge -= side*2
        kEdge = self.desorb6corner(pyGNe,kEdge,side,'-','C2H')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        #Taken from 31st R17 (reaction (k)) multiplied by a branching ratio
        #from migration with simultaneous desorption from integrated equations.
        T = pyGEn.temperature 
        coeff = 1.3e11*T**1.08*numpy.exp(-35428./T)
        k7 = pyGR_CsR5HtCsm().rate(pyGEn)
        k8 = pyGR_CsR5HtCsC2H().rate(pyGEn)
        k9 = pyGR_CsR5HtR5HCs().rate(pyGEn)
        br_ratio = k7/(k7+k8+k9)
        if T < 1500 or T > 2500:
            print 'Using Arhennius fit outside fitted temperature range.'
        return coeff*br_ratio
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1H',siteID[0]) and
           ((re.match(r'A1HA1-S1',siteID[1]) and re.match(r'A1HS1',siteID[2])) or
            (re.match(r'A1HA1-S1',siteID[2]) and re.match(r'A1HS1',siteID[1])))):
            return True
        return False

class pyGR_CsmpC2HtCsR6m(pyGReaction):
    #Dissertation Reaction 26
    isStructural = True
    number = 26
    name = 'CsmpC2HtCsR6m'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if (re.match(r'A1C2H',siteID[1]) and re.match(r'A1C2H',siteID[2])): #Both sides can adsorb 6 corner
            if random.random() > 0.5: #Go up
                side = +1
            else:                        #Go down
                side = -1
        elif re.match(r'A1C2H',siteID[1]):
            side = siteDir * -1
        elif re.match(r'A1C2H',siteID[2]):
            side = siteDir * 1
        else:
            raise ValueError('Invalid type for adsorbing 6-member corner ring.')
        kEdge = self.adsorb6corner(pyGNe,kEdge,side)
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        #Rates and steady state equation from Frenklach 26th Symposium 
        #(Reaction 11)
        T = pyGEn.temperature 
        cOfC2H2 = pyGEn.concOf('C2H2')
        kap = 1.1e7 * T**1.71*numpy.exp(-1960./T)
        kmap = 1.3e14*numpy.exp(-21024.8/T)
        kbp = 4.8e12*numpy.exp(-16875./T)
        kd = 2.5e12*T**-0.13*numpy.exp(-7902./T)
        coeff = kap*kd/(kmap+kbp+kd)
        return coeff*cOfC2H2
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1-',siteID[0]) and
           (re.match(r'A1C2H',siteID[1]) or
            re.match(r'A1C2H',siteID[2]))):
            return True
        return False

class pyGR_CsmtCsR6(pyGReaction):
    #Dissertation Reaction 27
    isStructural = True
    number = 27
    name = 'CsmtCsR6'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if (re.match(r'S1S1A1H',siteID[1]) and re.match(r'S1S1A1H',siteID[2])): #Both sides are boat
            if random.random() > 0.5: #Go up
                side = +1
            else:                        #Go down
                side = -1
        elif re.match(r'S1S1A1H',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S1S1A1H',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for adsorbing 6-member boat ring.')
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S1','S1')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        T = pyGEn.temperature
        cOfC2H2 = pyGEn.concOf('C2H2')
        #Coeff from 31st Symposium
        coeff = 8.0e7*T**1.56*numpy.exp(-1912./T);
        return coeff*cOfC2H2;
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1-',siteID[0]) and 
            (re.match(r'S1S1A1H',siteID[1]) or re.match(r'S1S1A1H',siteID[2]))):
            return True
        return False 

class pyGR_CsR5HpR6tCsR6R6(pyGReaction):
    #Dissertation Reaction 28
    isStructural = True
    number = 28
    name = 'CsR5HpR6tCsR6R6'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2S1S1A1H',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S2S1S1A1H',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for adsorbing 6-member boat ring.')
        kEdge = self.desorb5zigzag(pyGNe,kEdge,-side,'-','H')
        if side < 0:
           kEdge -= 1
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S1','S1')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsR5HtR5HCs().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A2H(?!2)',siteID[0]) and 
           ((re.match(r'S2S1S1A1H',siteID[1]) and re.match(r'A2H2S2',siteID[2])) or
            (re.match(r'S2S1S1A1H',siteID[2]) and re.match(r'A2H2S2',siteID[1])))):
            return True
        return False 

class pyGR_CsR5R6mtCsR6mR5E(pyGReaction):
    #Dissertation Reaction 29
    isStructural = True
    number = 29
    name = 'CsR5R6mtCsR6mR5E'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S2',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for external ring flip.')
        prob = 0.5
        if random.random() < prob:
            kEdge = pyGR_CsmtCs().apply(pyGNe,kEdge)
        else: 
            kEdge = self.flipR5R6(pyGNe,kEdge,side)
            pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsmtCs().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1-',siteID[0]) and 
           ((re.match(r'S2A2HS2',siteID[1]) and re.match(r'A1HS1(?!solid)',siteID[2])) or
            (re.match(r'S2A2HS2',siteID[2]) and re.match(r'A1HS1(?!solid)',siteID[1])))):
            return True
        return False 

class pyGR_CsR5R6mtCsR6mR5I(pyGReaction):
    #Dissertation Reaction 30
    isStructural = True
    number = 30
    name = 'CsR5R6mtCsR6mR5I'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S2',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for internal ring flip.')
        prob = 0.5
        if random.random() < prob:
            kEdge = pyGR_CsmtCs().apply(pyGNe,kEdge)
        else: 
            kEdge = self.flipR5R6(pyGNe,kEdge,side)
            pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsmtCs().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1-',siteID[0]) and 
           ((re.match(r'S2S2',siteID[1]) and re.match(r'S1A1',siteID[2])) or
            (re.match(r'S2S2',siteID[2]) and re.match(r'S1A1',siteID[1])))):
            return True
        return False 

class pyGR_CsR5R6mtCsR6mR5EI(pyGReaction):
    #Dissertation Reaction 31
    isStructural = True
    number = 31
    name = 'CsR5R6mtCsR6mR5EI'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S2',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for external to internal ring flip.')
        prob = 2./3. 
        if random.random() < prob:
            kEdge = pyGR_CsmtCs().apply(pyGNe,kEdge)
        else: 
            kEdge = self.flipR5R6(pyGNe,kEdge,side)
            pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsmtCs().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1-',siteID[0]) and 
           ((re.match(r'S2A2HS2',siteID[1]) and re.match(r'S1A1',siteID[2])) or
            (re.match(r'S2A2HS2',siteID[2]) and re.match(r'S1A1',siteID[1])))):
            return True
        return False 

class pyGR_CsR5R6mtCsR6mR5IE(pyGReaction):
    #Dissertation Reaction 32
    isStructural = True
    number = 32
    name = 'CsR5R6mtCsR6mR5IE'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S2',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for internal to external ring flip.')
        prob = 1./3. 
        if random.random() < prob:
            kEdge = pyGR_CsmtCs().apply(pyGNe,kEdge)
        else: 
            kEdge = self.flipR5R6(pyGNe,kEdge,side)
            pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsmtCs().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1-',siteID[0]) and 
           ((re.match(r'S2S2',siteID[1]) and re.match(r'A1HS1(?!solid)',siteID[2])) or
            (re.match(r'S2S2',siteID[2]) and re.match(r'A1HS1(?!solid)',siteID[1])))):
            return True
        return False 

class pyGR_CsmpR5tCsR6R5H(pyGReaction):
    #Dissertation Reaction 33
    isStructural = True
    number = 33
    name = 'CsmpR5tCsR6R5H'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if (re.match(r'S1S2A2H(?!2)',siteID[1]) and re.match(r'S1S2A2H(?!2)',siteID[2])): #Both sides are boat
            if random.random() > 0.5: #Go up
                side = +1
            else:                        #Go down
                side = -1
        elif re.match(r'S1S2A2H(?!2)',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S1S2A2H(?!2)',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for adsorbing 6-member boat ring.')
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S1','S2')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        T = pyGEn.temperature
        cOfC2H2 = pyGEn.concOf('C2H2')
        #Rate of ring adsorption to Csm + R5; taken as reverse
        #of ring separation to Csm.
        back_rate = pyGR_CsR6R5HtCsmpR5().rate(pyGEn)
        # From B3LYP/6-311G(d,p) data using Multiwell 2008.2 and scaling
        # frequencies by 0.9668 (units are Mol/CC)
        KeqDict = {1500:2.14E-12,1750:3.46E-10,2000:1.48E-08,2250:2.61E-07,2500:2.50E-06}
        if T in KeqDict:
            Keq = KeqDict[T]
        else:
            Keq = 3.24e3*numpy.exp(-52539./T)
            if T < 1500 or T > 2500:
                print 'Using Arhennius fit outside fitted temperature range.'
        return back_rate/Keq*cOfC2H2
        #return back_rate/Keq #Testing effect of this
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1-',siteID[0]) and 
            (re.match(r'S1S2A2H(?!2)',siteID[1]) or re.match(r'S1S2A2H(?!2)',siteID[2]))):
            return True
        return False 

class pyGR_CsR5mtCsR6R5H(pyGReaction):
    #Dissertation Reaction 34
    isStructural = True
    number = 34
    name = 'CsR5mtCsR6R5H'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if (re.match(r'S2S1A1H',siteID[1]) and re.match(r'S2S1A1H',siteID[2])): #Both sides are boat
            if random.random() > 0.5: #Go up
                side = +1
            else:                        #Go down
                side = -1
        elif re.match(r'S2S1A1H',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S2S1A1H',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for adsorbing 6-member boat ring.')
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S2','S1')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsmpR5tCsR6R5H().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A2-',siteID[0]) and 
            (re.match(r'S2S1A1H',siteID[1]) or re.match(r'S2S1A1H',siteID[2]))):
            return True
        return False 

class pyGR_CsmpR5StCsR6(pyGReaction):
    #Dissertation Reaction 35
    isStructural = True
    number = 35
    name = 'CsmpR5StCsR6'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if (re.match(r'S2S2A1H',siteID[1]) and re.match(r'S2S2A1H',siteID[2])): #Both sides are boat
            if random.random() > 0.5: #Go up
                side = +1
            else:                        #Go down
                side = -1
        elif re.match(r'S2S2A1H',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S2S2A1H',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for adsorbing 6-member boat ring.')
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S1','S1')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsmtCsR5H().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1-',siteID[0]) and 
            (re.match(r'S2S2A1H',siteID[1]) or re.match(r'S2S2A1H',siteID[2]))):
            return True
        return False 

class pyGR_CsC2HpR5StCsR6(pyGReaction):
    #Dissertation Reaction 36
    isStructural = True
    number = 36
    name = 'CsC2HpR5StCsR6'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if (re.match(r'S2S2A1H',siteID[1]) and re.match(r'S2S2A1H',siteID[2])): #Both sides are boat
            if random.random() > 0.5: #Go up
                side = +1
            else:                        #Go down
                side = -1
        elif re.match(r'S2S2A1H',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S2S2A1H',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for adsorbing 6-member boat ring.')
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S1','S1')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsC2HtCsR5H().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1C2H',siteID[0]) and 
            (re.match(r'S2S2A1H',siteID[1]) or re.match(r'S2S2A1H',siteID[2]))):
            return True
        return False 

class pyGR_CsR5HpR5StCsR6(pyGReaction):
    #Dissertation Reaction 37
    isStructural = True
    number = 37
    name = 'CsR5HpR5StCsR6'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2S2S2A1H',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S2S2S2A1H',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for adsorbing 6-member boat ring.')
        kEdge = self.desorb5zigzag(pyGNe,kEdge,-side,'-','H')
        if side < 0:
           kEdge -= 1
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S1','S1')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsR5HtR5HCs().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A2H(?!2)',siteID[0]) and 
           ((re.match(r'S2S2S2A1H',siteID[1]) and re.match(r'A2H2S2',siteID[2])) or
            (re.match(r'S2S2S2A1H',siteID[2]) and re.match(r'A2H2S2',siteID[1])))):
            return True
        return False 

class pyGR_CsR6R5HpR5StCsR6pR5(pyGReaction):
    #Dissertation Reaction 38
    isStructural = True
    number = 38
    name = 'CsR6R5HpR5StCsR6pR5'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S2',siteID[1]):
            side = siteDir * 1
        elif re.match(r'S2',siteID[2]):
            side = siteDir * -1
        else:
            raise ValueError('Invalid type for migrating 5-member ring.')
        removeHsiteIndex = pyGNe.edge.nodes[kEdge - side*2]
        pyGNe.nodes[removeHsiteIndex].ncn = 'H'
        kEdge = self.desorb6boat(pyGNe,kEdge,side,'H','-')
        kEdge += side*2
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S1','S1')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsR6R5HtCsR5HpR5().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1H',siteID[0]) and
           ((re.match(r'S2A2H2S2',siteID[1]) and re.match(r'A1HS1S2S2A1H',siteID[2])) or
            (re.match(r'S2A2H2S2',siteID[2]) and re.match(r'A1HS1S2S2A1H',siteID[1])))):
            return True
        return False

class pyGR_CsR5mpR5tCsR6(pyGReaction):
    #Dissertation Reaction 39
    isStructural = True
    number = 39
    name = 'CsR5mpR5tCsR6'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if (re.match(r'S2S2A2HA2H(?!2)',siteID[1]) and re.match(r'S2S2A2HA2H(?!2)',siteID[2])): #Both sides are boat
            if random.random() > 0.5: #Go up
                side = +1
            else:                        #Go down
                side = -1
        elif re.match(r'S2S2A2HA2H(?!2)',siteID[1]): #First side is boat 
            side = siteDir * -1        
        elif re.match(r'S2S2A2HA2H(?!2)',siteID[2]): #Second side is boat
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for adsorbing 6-member boat ring.')
        kEdge = self.adsorb6boat(pyGNe,kEdge,side,'S2','S2')
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        return pyGR_CsmpR5tCsR6R5H().rate(pyGEn)
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A2-',siteID[0]) and 
            (re.match(r'S2S2A2HA2H(?!2)',siteID[1]) or re.match(r'S2S2A2HA2H(?!2)',siteID[2]))):
            return True
        return False 

class pyGR_CsR6BtCsR6S(pyGReaction):
    #Dissertation Reaction 40
    isStructural = True
    number = 40
    name = 'CsR6BtCsR6S'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S1S1S1S1A1H',siteID[1]): #First side is bay
            side = siteDir * -1        
        elif re.match(r'S1S1S1S1A1H',siteID[2]): #Second side is bay
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for 6-member ring bay closure.')
        if (kEdge < 5 and side < 0): nTurned = 5
        elif (kEdge > (pyGNe.edge.length-6) and side > 0): nTurned = -5
        else: nTurned = 0
        pyGNe.edge.turn(nTurned)
        kEdge += nTurned
        i_0 = pyGNe.edge.nodes[kEdge]
        ip5 = pyGNe.edge.nodes[kEdge+side*5]
        pyGNe.nodes[i_0].cns.append(ip5)
        pyGNe.nodes[ip5].cns.append(i_0)
        pyGNe.nodes[i_0].ncn = ''
        pyGNe.nodes[ip5].ncn = ''
        pyGNe.nodes[i_0].type = 'S1'
        pyGNe.nodes[ip5].type = 'S1'
        if side < 0:
            i_pop = kEdge - 4
        else:
            i_pop = kEdge + 1
        for k in range(4):
            pyGNe.edge.nodes.pop(i_pop)
        #pyGNe.edge.turn(-nTurned)
        pyGNe.nR6 += 1
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        #Rate for R6 Bay Closure
        #Taken from M.Kraft Preprint 52 
        T = pyGEn.temperature
        Rcal = 1.987;  #Universal gas constant in Cal/mol/K
        #3.49 ? 10^12 ? T^?0.39 ? exp(-2.440 kcal/mol/(RT))
        return 3.49e12*T**-0.39*numpy.exp(-2440./(Rcal*T));
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1-',siteID[0]) and 
            (re.match(r'S1S1S1S1A1H',siteID[1]) or re.match(r'S1S1S1S1A1H',siteID[2]))):
            return True
        return False 

class pyGR_CsR5BtCsR5S(pyGReaction):
    #Dissertation Reaction 41
    isStructural = True
    number = 41
    name = 'CsR5BtCsR5S'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'S1S1S1A1H',siteID[1]): #First side is bay
            side = siteDir * -1        
        elif re.match(r'S1S1S1A1H',siteID[2]): #Second side is bay
            side = siteDir * +1
        else:
            raise ValueError('Invalid type for 6-member ring bay closure.')
        if (kEdge < 4 and side < 0): nTurned = 4
        elif (kEdge > (pyGNe.edge.length-5) and side > 0): nTurned = -4
        else: nTurned = 0
        pyGNe.edge.turn(nTurned)
        kEdge += nTurned
        i_0 = pyGNe.edge.nodes[kEdge]
        ip4 = pyGNe.edge.nodes[kEdge+side*4]
        pyGNe.nodes[i_0].cns.append(ip4)
        pyGNe.nodes[ip4].cns.append(i_0)
        pyGNe.nodes[i_0].ncn = ''
        pyGNe.nodes[ip4].ncn = ''
        pyGNe.nodes[i_0].type = 'S2'
        pyGNe.nodes[ip4].type = 'S2'
        if side < 0:
            i_pop = kEdge - 3
        else:
            i_pop = kEdge + 1
        for k in range(3):
            pyGNe.edge.nodes.pop(i_pop)
        #pyGNe.edge.turn(-nTurned)
        pyGNe.nR5 += 1
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        #%Rate for R5 Bay Closure
        #%Taken from Violi Embedded Migration Paper
        #%A. Violi, J. Phys. Chem. A 109 (2005) 7781-7787.
        T = pyGEn.temperature
        Rcal = 1.987  #Universal gas constant in Cal/mol/K
        #3.86 ? 10^11 ? T^-0.21 ? exp(-17.70 kcal/mol/(RT))
        return 3.86e11*T**-0.21*numpy.exp(-17700./(Rcal*T))
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A1-',siteID[0]) and 
            (re.match(r'S1S1S1A1H',siteID[1]) or re.match(r'S1S1S1A1H',siteID[2]))):
            return True
        return False 

class pyGR_CsR5mtCsR6(pyGReaction):
    #Dissertation Reaction 42
    isStructural = True
    number = 42
    name = 'CsR5mtCsR6'
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        siteID = pyGNe.edge.idStrings[kEdge].split('_')
        siteDir  = pyGNe.edge.directions[kEdge]
        if re.match(r'A2H',siteID[1]): #First side is bay
            side = siteDir * -1        
            SS = False
        elif re.match(r'A2H',siteID[2]): #Second side is bay
            side = siteDir * +1
            SS = False
        elif (re.match(r'S2',siteID[1]) and re.match(r'S2',siteID[2])):
            side = +1   #If at a corner then side doesn't matter
            SS = True
        else:
            raise ValueError('Invalid type for 6-member ring bay closure.')
        if (kEdge < 2 and side < 0): nTurned = 2
        elif (kEdge > (pyGNe.edge.length-3) and side > 0): nTurned = -2
        else: nTurned = 0
        pyGNe.edge.turn(nTurned)
        kEdge += nTurned
        im1 = pyGNe.edge.nodes[kEdge-side]
        i_0 = pyGNe.edge.nodes[kEdge]
        ip1 = pyGNe.edge.nodes[kEdge+side]
        ip2 = pyGNe.edge.nodes[kEdge+side*2]
        if SS:
            i_U1 = numpy.setdiff1d(pyGNe.nodes[ip1].cns,pyGNe.edge.nodes)
            i_U2 = numpy.setdiff1d(pyGNe.nodes[im1].cns,pyGNe.edge.nodes)
            (posA,posB) = pyGNe.nodes[ip1].setR6pos(pyGNe.nodes[i_U1],
                                                pyGNe.nodes[i_U2],
                                                pyGNe.nodes[im1])
        else:
            i_Under=numpy.intersect1d(pyGNe.nodes[im1].cns,pyGNe.nodes[ip2].cns)
            (posA,posB) = pyGNe.nodes[ip1].setR6pos(pyGNe.nodes[ip2],
                                                pyGNe.nodes[i_Under],
                                                pyGNe.nodes[im1])
        pyGNe.nodes[i_0].cns.remove(ip1)
        pyGNe.nodes[ip1].cns.remove(i_0)
        pyGNe.nodes[i_0].cns.append(pyGNe.nC)
        pyGNe.nodes[ip1].cns.append(pyGNe.nC)
        pyGNe.nodes[i_0].type = 'A1'
        pyGNe.nodes[im1].type = 'S1'
        if SS:
            pyGNe.nodes[ip1].type = 'S1'
        else:
            pyGNe.nodes[ip1].type = 'A1'
            pyGNe.nodes[ip2].type = 'S1'
        pyGNe.nodes[i_0].ncn = 'H'
        pyGNe.nodes[i_0].pos = posB
        pyGNe.nodes.append(pyGNode(pos=posA,cns=[i_0,ip1],ncn='H',type='A1'))
        if side < 0:
            pyGNe.edge.nodes.insert(kEdge,pyGNe.nC)
        else:
            pyGNe.edge.nodes.insert(kEdge+1,pyGNe.nC)
        #pyGNe.edge.turn(-nTurned)
        pyGNe.nR6 += 1
        pyGNe.nR5 -= 1
        pyGNe.nC += 1
        pyGNe.edge.length += 1
        pyGNe.classifyEdge()
        return kEdge
    def rate(self,pyGEn):
        #Rate for converting five-member ring to six member ring
        #Assumed to be high to test affect.
        cOfCH3 = pyGEn.concOf('CH3')
        return 1.0e13*cOfCH3
    
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if (re.match(r'A2-',siteID[0]) and 
            (re.match(r'[AS]2',siteID[1]) or re.match(r'[AS]2',siteID[2]))):
            return True
        return False

class pyGR_CsmtCsC2H3(pyGReaction):
    #Dissertation Reaction 43
    isStructural = True
    number = 43
    name = 'CsmtCsC2H3'
    count = 0
    def __init__(self):
        pyGReaction.__init__(self)
    def apply(self,pyGNe,kEdge):
        pyGR_CsmtCsC2H3.count += 1
        print 'haha', pyGR_CsmtCsC2H3.count
        kEdge = self.CAHM(pyGNe,kEdge,'S1')
        # kNode = pyGNe.edge.nodes[kEdge]     
        pyGNe.classifyEdge()
        # pyGNe.nCCAHM += 2
        return kEdge
    def rate(self,pyGEn):
        print 'hahahaha'
        T = pyGEn.temperature
        rateDict = {1500:5.54E+05, 2000:2.78E+07}
        cOfC2H2 = pyGEn.concOf('C2H2')
        if T in rateDict:
            rate = rateDict[T]
        else:
            #Interpolate (See JPC KMC Paper).
            rate = 4.07*T**3.25*numpy.exp(-17865./T)*cOfC2H2
        if T < 1500 or T > 2500:
            print 'Using Arhennius fit outside fitted temperature range.'
        return rate 
    def isApplicable(self,siteType):
        siteID = siteType.split('_')
        if re.match(r'A1H',siteID[0]) and re.match(r'S1A1HS1',siteID[1]) and re.match(r'S1A1HS1',siteID[2]):
        #and m1 == 1:
            return True
        return False
# class pyGR_CsmtCsC2H3(pyGReaction):
#     #Dissertation Reaction 43
#     isStructural = False
#     number = 43
#     name = 'CsmtCsC2H3'
#     count = 0
#     def __init__(self):
#         pyGReaction.__init__(self)
#     def apply(self,pyGNe,kEdge):
#         pyGR_CsmtCsC2H3.count += 1
#         print 'haha', pyGR_CsmtCsC2H3.count
#         kNode = pyGNe.edge.nodes[kEdge]
#         pyGNe.nodes[kNode].ncn = 'C2H3'
#         pyGNe.classifyEdge(kEdge=kEdge)
#         pyGNe.nCCAHM += 2
#         return kEdge
#     def rate(self,pyGEn):
#         print 'hahahaha'
#         T = pyGEn.temperature
#         rateDict = {1500:5.54E+05, 2000:2.78E+07}
#         cOfC2H2 = pyGEn.concOf('C2H2')
#         if T in rateDict:
#             rate = rateDict[T]
#         else:
#             #Interpolate (See JPC KMC Paper).
#             rate = 4.07*T**3.25*numpy.exp(-17865./T)*cOfC2H2
#         if T < 1500 or T > 2500:
#             print 'Using Arhennius fit outside fitted temperature range.'
#         return rate 
#     def isApplicable(self,siteType):
#         siteID = siteType.split('_')
#         if re.match(r'A1H',siteID[0]) and re.match(r'S1A1HS1',siteID[1]) and re.match(r'S1A1HS1',siteID[2]):
#         #and m1 == 1:
#             return True
#         return False

def pyGRLoad(excludedReactions=[],objList=dir()):
        reactions = []
        for m in range(len(objList)):
            try:
                if ((not objList[m] == 'pyGReaction') and
                     issubclass(eval(objList[m]),pyGReaction)):
                    tmpReaction = eval(objList[m]+'()')
                    if not (tmpReaction.number in excludedReactions):
                          reactions.append(tmpReaction)
            except:
                pass
        return reactions

if __name__ == '__main__':
    from pyGEdge import pyGEdge
    from pyGEnvironment import pyGEnvironment
    from pyGNetwork import pyGNetwork
    import matplotlib.pyplot as plt
    env = pyGEnvironment(temp=1500)
#    print pyGR_CstCsm().rate(env)
    pyGNe = pyGNetwork(baseType='linear',nRings=7)
    rx1 = pyGR_CstCsm()
    rx2 = pyGR_CsmtCs()
    rx3 = pyGR_CsmtCsR5H()
    rx4 = pyGR_CsmtCsC2H()
    rx5 = pyGR_CsC2HtCsm()
    rx6 = pyGR_CsC2HtCsR5H()
    rx7 = pyGR_CsR5HtCsm()
    rx8 = pyGR_CsR5HtCsC2H()
    rx9 = pyGR_CsR5HtR5HCs()
    rx10 = pyGR_CsR5tCsR5m()
    rx11 = pyGR_CsR5mtCsR5()
    rx12 = pyGR_CsR5tCsR5H()
    rx13 = pyGR_CsR5HtCsR5()
    rx14 = pyGR_CsR5HpR5tCsR6R5H()
    rx15 = pyGR_CsR6R5HtCsR5HpR5()
    rx16 = pyGR_CsR6R5HtCsC2HpR5()
    rx17 = pyGR_CsR6R5HtCsmpR5()
    rx18 = pyGR_CsR6R5HpR5tCsR5HR6pR5()
    rx19 = pyGR_CsR6R5HpR5R6tCsR5SR6pR5()
    rx20 = pyGR_CsR6R5HpR6tCsR6R6pR5()
    rx21 = pyGR_CsR5HpC2HtCsR6m()
    rx22 = pyGR_CsR6R5HpC2HtCsR6mpR5()
    rx23 = pyGR_CsR6mtCsR5HpC2H()
    rx24 = pyGR_CsR6mtCsC2HpC2H()
    rx25 = pyGR_CsR6mtCsmpC2H()
    rx26 = pyGR_CsmpC2HtCsR6m()
    rx27 = pyGR_CsmtCsR6()
    rx28 = pyGR_CsR5HpR6tCsR6R6()
    rx29 = pyGR_CsR5R6mtCsR6mR5E()
    rx30 = pyGR_CsR5R6mtCsR6mR5I()
    rx31 = pyGR_CsR5R6mtCsR6mR5EI()
    rx32 = pyGR_CsR5R6mtCsR6mR5IE()
    rx33 = pyGR_CsmpR5tCsR6R5H()
    rx34 = pyGR_CsR5mtCsR6R5H()
    rx35 = pyGR_CsmpR5StCsR6()
    rx36 = pyGR_CsC2HpR5StCsR6()
    rx37 = pyGR_CsR5HpR5StCsR6() 
    rx38 = pyGR_CsR6R5HpR5StCsR6pR5()
    rx39 = pyGR_CsR5mpR5tCsR6()
    rx40 = pyGR_CsR6BtCsR6S()
    rx41 = pyGR_CsR5BtCsR5S()
    rx42 = pyGR_CsR5mtCsR6()
    rx43 = pyGR_CsmtCsC2H3()

    reactions = pyGRLoad()
    rx_0 = reactions[0]

#    for k in range(len(reactions)):
#        print "{0} {name}  {1} {2}".format(reactions[k].number,reactions[k].rate(env),env.temperature,name=reactions[k].name)

    #rx1.apply(pyGNe,7)
    #rx1.apply(pyGNe,10)
    #rx4.apply(pyGNe,7)
    #print rx26.isApplicable(pyGNe.edge.idStrings[9])
    rx43.apply(pyGNe,9)
    # print pyGNe.edge.idStrings[7]

    #print rx23.isApplicable(pyGNe.edge.idStrings[11])
    #rx10.apply(pyGNe,8)
    #rx34.apply(pyGNe,8)
    #rx10.apply(pyGNe,9)
    #print rx42.isApplicable(pyGNe.edge.idStrings[9])
    #rx42.apply(pyGNe,9)
    # pyGNe.plot()

    # plt.show()
    print 'Done.'

    
