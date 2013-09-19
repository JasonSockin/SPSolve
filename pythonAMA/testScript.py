import numpy
from Amalg import Amalg

h = numpy.matrix([[0,0,1.2,0,-1,-1,0,-0.01],[0,-.7,0,1,0,0,0,0]])
lead, lag, neq, condn, upbnd = 2, 1, 2, 0.00000001, 1.00000001

Amalg(h,neq,lag,lead,condn,upbnd)
