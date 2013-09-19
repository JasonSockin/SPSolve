# This file contains the following functions:
#          Exact_shift, Numeric_shift, Shiftright

# Original author: Gary Anderson
# Original file downloaded from:
# http://www.federalreserve.gov/Pubs/oss/oss4/code.html
# Adapted for Dynare by Dynare Team.

# This code is in the public domain and may be used freely.
# However the authors would appreciate acknowledgement of the source by
# citation of any of the following papers:

# Anderson, G. and Moore, G.
# "A Linear Algebraic Procedure for Solving Linear Perfect Foresight
# Models."
# Economics Letters, 17, 1985.

# Anderson, G.
# "Solving Linear Rational Expectations Models: A Horse Race"
# Computational Economics, 2008, vol. 31, issue 2, pages 95-113

# Anderson, G.
# "A Reliable and Computationally Efficient Algorithm for Imposing the
# Saddle Point Property in Dynamic Models"
# Journal of Economic Dynamics and Control, 2010, vol. 34, issue 3,
# pages 472-489


# Import the numpy and scipy packages
from numpy import *
from scipy import *
from scipy.sparse import *
from scipy.linalg import *

##########################################################################


def Shiftright(x,n):

    #  Shift the rows of x to the right by n columns, leaving zeros in the
    #  first n columns. 
    
    rows, cols = x.shape
    left  = range(0,cols-n)  
    right = range(n,cols) 
    
    y = matrix(zeros(shape=((rows,cols))))
    y[ix_(range(0,rows),right)] = x[range(0,rows),:][:,left].toarray()
    
    return y

##########################################################################


def exactShift(h,q,iq,qrows,qcols,neq):

    # Compute the exact shiftrights and store them in q.

    hs = lil_matrix(h)
    
    nexact = 0
    left   = range(0,qcols)
    right  = range(qcols,qcols+neq)
    zerorows = list()
    sumVector = abs(hs[:,right].T).sum(axis=1)
    for i in range(0,len(sumVector)):
        if sumVector[i] == 0:
            zerorows.append(i)

    while len(zerorows) > 0  and iq <= qrows:
        nz = len(zerorows)
        hsRows, hsCols = hs.get_shape()
        q[ix_(range(iq,iq+nz),range(0,qcols))] = hs[zerorows,:][:,left].toarray()
        hs[zerorows,:] = Shiftright( hs[zerorows,:], neq )
        iq = iq + nz
        nexact = nexact + nz
        zerorows = list()
        newSumVector = abs(hs[:,right].T).sum(axis=1)
        for i in range(0,len(newSumVector)):
            if newSumVector[i] == 0:
                zerorows.append(i)
        h = lil_matrix(hs).todense()

    return h, q, iq, nexact

#########################################################################

def numericShift(h,q,iq,qrows,qcols,neq,condn):

    # Compute the numeric shiftrights and store them in q.
    
    nnumeric = 0
    left     = range(0,qcols)
    right    = range(qcols,qcols+neq)
    
    Q, R  = qr(h[:,right])
    zerorows = list()
    testVector = abs(diagonal(R))
    for i in range(0,len(testVector)):
        if testVector[i] <= condn:
            zerorows.append(i)
    
    while len(zerorows) > 0 and iq <= qrows:
        h = lil_matrix(h)
        Q = lil_matrix(Q)
        h = lil_matrix(Q.T * h)
        nz = len(zerorows)
        q[ix_(range(iq,iq+nz),range(0,qcols))] = h[zerorows,:][:,left].toarray()
        h[zerorows,:] = Shiftright( h[zerorows,:], neq )
        iq = iq + nz
        nnumeric = nnumeric + nz
        Q, R  = qr(lil_matrix(h[:,right]).todense())
        zerorows = list()
        testVector = abs(diagonal(R))
        for i in range(0,len(testVector)):
            if testVector[i] <= condn:
                zerorows.append(i)
                
    return h, q, iq, nnumeric
