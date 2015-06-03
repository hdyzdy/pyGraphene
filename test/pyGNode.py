

import numpy

class pyGNode:
    def __init__(self,cns=[],ncn='',pos=numpy.array([0,0,0]),type=''):
        self.cns = cns
        self.ncn = ncn
        self.pos = pos
        self.type = type
    def isActive(self):
        bool = False
        if (len(self.cns) < 3) and (self.ncn != 'solid'):
            bool = True
        return bool

    def setCAHMpos(self,node2,node3):
        p1 = self.pos
        p2 = node2.pos
        p3 = node3.pos
        v2 = p1 - p2
        v3 = p1 - p3
        v4 = v2 + v3
        v4u = v4/numpy.linalg.norm(v4)
        p4 = p1 + 1.45*v4u
        p5 = p4 + 1.3*v4u
        # print v4,v4u
        # print p4,p5
        return (p4, p5)

    def setR6pos(self,node2,node3,node4):
        #This function sets the positions of a the acetylene atoms in a HACA step
        #by mirroring the existing middle nodes (2,3) across the axis defined by
        #the outer nodes (1,4)
        p1 = self.pos
        p2 = node2.pos
        p3 = node3.pos
        p4 = node4.pos
        v1 = p4 - p1
        v1u = v1/numpy.linalg.norm(v1)
        v2 = p2 - p1
        v3 = p3 - p1
        v5a = p1 + numpy.dot(v1u,v2)*v1u
        v6a = p1 + numpy.dot(v1u,v3)*v1u
        p5 = v5a + v5a - p2
        p6 = v6a + v6a - p3
        return (p5,p6)
    def setR5pos(self,node2,node3):
        p1 = self.pos
        p2 = node2.pos
        p3 = node3.pos
        vX = (p3 - p2)/numpy.linalg.norm(p3 - p2)
        vY = (p1 - p2)/numpy.linalg.norm(p1 - p2) * numpy.sin(numpy.pi*2/5)**-1 + vX*numpy.tan(numpy.pi*2/5)**-1
        D = 2.*numpy.linalg.norm(p3-p2)*numpy.cos(numpy.pi/5)
        p4 = p2 + D*(numpy.cos(numpy.pi*2/5)*vX + numpy.sin(numpy.pi*2/5)*vY)
        p5 = p2 + D*(numpy.cos(numpy.pi/5)*vX + numpy.sin(numpy.pi/5)*vY)
        return (p4,p5)

if __name__ == '__main__':
    # node1 = pyGNode(pos=numpy.array([-numpy.cos(numpy.pi*2./5),numpy.sin(numpy.pi*2./5),0.]))
    # node2 = pyGNode(pos=numpy.array([0.,0.,0.]))
    # node3 = pyGNode(pos=numpy.array([1.,0.,0.]))
    # (p4,p5) = node1.setR5pos(node2,node3)

    # print p4, p5

    node1 = pyGNode(pos=numpy.array([0.,0.,0.]))
    node2 = pyGNode(pos=numpy.array([-1.42*numpy.sin(numpy.pi/6.),-1.42*numpy.cos(numpy.pi/6.),0.]))
    node3 = pyGNode(pos=numpy.array([1.42*numpy.sin(numpy.pi/6.),-1.42*numpy.cos(numpy.pi/6.),0.]))
    (p4,p5) = node1.setCAHMpos(node2,node3)
    print p4, p5



