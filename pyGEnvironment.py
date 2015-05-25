
class pyGEnvironment:
    def __init__(self,temp=1500,press=101325,xOfH=0.00001,xOfH2=0.1,xOfC2H2=0.1,xOfCH3=0.0):
        self.temperature = temp;
        self.pressure = press;
        self.moleFrac = { 'H' : xOfH, 'H2' : xOfH2, 'C2H2' : xOfC2H2, 'CH3' : xOfCH3 }
    def concOf(self,species):
        totalC = self.pressure/(self.temperature * 8.314) / 1.0e6 #Concentration is in mol/cc
        return self.moleFrac[species]*totalC


if __name__ == '__main__':
    pge = pyGEnvironment()
    print pge.concOf('C2H2')

