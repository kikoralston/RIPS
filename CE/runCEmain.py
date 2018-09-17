from RIPSMasterScript17Nov2017 import *
from Parameters import *

print('Loading parameters and setting up initial data')
# Load parameters
genparam = Generalparameters.Generalparameters()
genparam.load(fname='./generalparameters.txt')

reserveparam = Reserveparameters.Reserveparameters()
reserveparam.load(fname='./reserveparameters.txt')

curtailparam = Curtailmentparameters.Curtailmentparameters()
curtailparam.load(fname='./curtailmentparameters.txt')

print()
print('Parameters loaded! Calling CE masterFunction...')
print()
print()
masterFunction(genparam, reserveparam, curtailparam)
