
class pyGEdge:
    def __init__(self,nodes=[],idStrings=[],directions=[]):
        self.nodes = pyGCircList(nodes)
        self.idStrings = pyGCircList(idStrings)
        self.directions = pyGCircList(directions)
        self.str = pyGCircList('')
        self.actives = pyGCircList([])
        self.length = 0
    def turn(self,n=1):
        self.nodes.turn(n)
        self.idStrings.turn(n)
        self.directions.turn(n)
        #Last two are for classifying edge and will be reinitialized after turning
        self.str = pyGCircList('')
        self.actives = pyGCircList([])

class pyGCircList(list):
    def __init__(self,data):
        super(pyGCircList, self).__init__(data) #Call parent init

    def __getitem__(self,index):
        try:
            val = super(pyGCircList,self).__getitem__(index)
        except:
            if index > len(self) - 1:
                nTurns = index-(len(self)-1)
                self.turn(-nTurns)
                index -= nTurns
                val = super(pyGCircList,self).__getitem__(index)
                self.turn(nTurns)
            elif index < 0:
                nTurns = 0 - index
                self.turn(nTurns)
                index += nTurns
                val = super(pyGCircList,self).__getitem__(index)
                self.turn(-nTurns)
            else:
                raise
        return val

    def __getitem_safe__(self,index):
        #Used when caller knows that index is in the right range
        return super(pyGCircList,self).__getitem__(index)

    def __setitem__(self,index,value):
        if index > len(self) - 1:
            nTurns = index-(len(self)-1)
            self.turn(-nTurns)
            index -= nTurns
            super(pyGCircList,self).__setitem__(index,value)
            self.turn(nTurns)
        elif index < 0:
            nTurns = 0 - index
            self.turn(nTurns)
            index += nTurns
            super(pyGCircList,self).__setitem__(index,value)
            self.turn(-nTurns)
        else:
            super(pyGCircList,self).__setitem__(index,value)

    def circgetslice(self,i,j):
        if (abs(j-i) > len(self)): raise ValueError("Slice larger than array")
        if i > j : j += len(self)
        vals = []
        if j > len(self) - 1:
            nTurns = -(j-(len(self)-1))
            self.turn(nTurns)
            i += nTurns
            j += nTurns
        elif i < 0:
            nTurns = 0 - i
            self.turn(nTurns)
            i += nTurns
            j += nTurns
        else:
            nTurns = 0
        vals = super(pyGCircList,self).__getitem__(slice(i,j))
        self.turn(-nTurns)
        return vals

    def circsetslice(self,i,j,sequence):
        if i > j : j += len(self)
        if ((j-i) > len(self)): raise ValueError("Slice larger than array")        
        if not ((j-i) == len(sequence)): raise ValueError("Slice and sequence sizes don't match")
        for m in range(i,j):
            self.__setitem__(m,sequence[m-i])

    def __getslice__(self,i,j):
        raise IndexError("Use circgetslice(i,j) instead.")
        #Have to do this because getslice is no longer supported for user classes
        #and fiddles with i's and j's before passing them to overloaded function.

    def __setslice__(self,i,j,sequence):
        raise IndexError("Use circsetslice(i,j,sequence) instead.")
        #See __getslice__() explanation above.

    def pop(self,index):
        if index > len(self) - 1:
            nTurns = index-(len(self)-1)
            self.turn(-nTurns)
            index -= nTurns
            val = super(pyGCircList,self).pop(index)
            self.turn(nTurns-1)
        elif index < 0:
            nTurns = 0 - index
            self.turn(nTurns)
            index += nTurns
            val = super(pyGCircList,self).pop(index)
            self.turn(-nTurns)
        else:
            val = super(pyGCircList,self).pop(index)
        return val

    def insert(self,index,value):
        if index > len(self) - 1:
            nTurns = index-(len(self)-1)
            self.turn(-nTurns)
            index -= nTurns
            val = super(pyGCircList,self).insert(index,value)
            self.turn(nTurns+1)
        elif index < 0:
            nTurns = 0 - index
            self.turn(nTurns)
            index += nTurns
            val = super(pyGCircList,self).insert(index,value)
            self.turn(-nTurns-1)
        else:
            val = super(pyGCircList,self).insert(index,value)
        return val

    def turn(self,n=1):
        super(pyGCircList, self).__setitem__(slice(0,None),
                                             super(pyGCircList,self).__getitem__(slice(-n,None)) + 
                                             super(pyGCircList,self).__getitem__(slice(0,-n)))    

if __name__ == '__main__':
    cl = pyGCircList(('A','B','C','D','E','F'))
    x = -2
    y = 1
    print cl.circgetslice(x,y)
    print cl
    cl.turn(-1)
    print cl
    #cl.circsetslice(x,y,['G','H','I','J'])
    cl.pop(8)
    print cl
    cl.append('G')
    print cl
    cl.insert(6,'H')
    print cl
    cl.turn(-3)
    print cl
    cl.turn(3)
    print cl 
    print cl[-2]
    print cl[8]
    print 'Done.' 
