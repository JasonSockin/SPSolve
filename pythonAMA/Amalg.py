# Import the numpy package
from numpy import *

# Import necessary user-defined functions
from Shifts import exactShift, numericShift
from buildA import buildA
from Eigensystem import Eigensystem
from copyW import copyW
from reducedForm import reducedForm
from checkA import existsNaN, existsInf
from makePhi import makePhi
from makeF import makeF

def Amalg(h,neq,nlag,nlead,condn,uprbnd):

    #  Solve a linear perfect foresight model using the matlab eig
    #  function to find the invariant subspace associated with the big
    #  roots.  This procedure will fail if the companion matrix is
    #  defective and does not have a linearly independent set of
    #  eigenvectors associated with the big roots.
    
    #  Input arguments:
    # 
    #    h         Structural coefficient matrix (neq,neq*(nlag+1+nlead)).
    #    neq       Number of equations.
    #    nlag      Number of lags.
    #    nlead     Number of leads.
    #    condn     Zero tolerance used as a condition number test
    #              by numeric_shift and reduced_form.
    #    uprbnd    Inclusive upper bound for the modulus of roots
    #              allowed in the reduced form.
    
    #  Output arguments:
    # 
    #    b         Reduced form coefficient matrix (neq,neq*nlag).
    #    rts       Roots returned by eig.
    #    ia        Dimension of companion matrix (number of non-trivial
    #              elements in rts).
    #    nexact    Number of exact shiftrights.
    #    nnumeric  Number of numeric shiftrights.
    #    lgroots   Number of roots greater in modulus than uprbnd.
    #    aimcode     Return code: see function aimerr.
    
    # Original author: Gary Anderson
    # Original file downloaded from:
    # http://www.federalreserve.gov/Pubs/oss/oss4/code.html
    # Adapted for Dynare by Dynare Team, in order to deal
    # with infinite or nan values.
    
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


    if nlag < 1 or nlead < 1: 
        error('Aim_eig: model must have at least one lag and one lead.')

    # Initialization.
    nexact = 0
    nnumeric = 0
    lgroots = 0
    iq = 0
    aimcode = 0
    b = 0
    qrows = neq*nlead
    qcols = neq*(nlag+nlead)
    bcols = neq*nlag
    q = matrix(zeros(shape=((qrows,qcols))))
    rts = matrix(zeros(shape=((qcols,1))))
    originalH = h

    # Compute the auxiliary initial conditions and store them in q.

    h, q, iq, nexact = exactShift(h,q,iq,qrows,qcols,neq)
    if iq > qrows: 
        aimcode = 61
        return
    
    h, q, iq, nnumeric = numericShift(h,q,iq,qrows,qcols,neq,condn)
    if iq > qrows: 
        aimcode = 62
        return

    #  Build the companion matrix.  Compute the stability conditions, and
    #  combine them with the auxiliary initial conditions in q.  

    a, ia, js = buildA(h,qcols,neq)

    if ia != 0:
        if existsNaN(a) or existsInf(a) : 
            print("A is NAN or INF")
            aimcode=63 
            return 
        w, rts, lgroots = Eigensystem(a,uprbnd,min(len(js),qrows-iq+1))
        q = copyW(q,w,js,iq,qrows)

    test = nexact + nnumeric + lgroots
    if test > qrows:
        aimcode = 3
    elif test < qrows:
        aimcode = 4

    phi = makePhi(q,originalH,nlag,nlead,neq)
    F = makeF(phi,originalH,q,nlag,nlead,neq)
    print()
    print("phi =")
    print()
    print(phi)
    print()
    print()
    print("F = ")
    print()
    print(F)
    print()

    # If the right-hand block of q is invertible, compute the reduced form.

    if aimcode == 0:
        nonsing, b = reducedForm(q,qrows,qcols,bcols,neq,condn)
    
    if nonsing and aimcode == 0:
        aimcode =  1
    elif not(nonsing) and aimcode == 0:
        aimcode =  5
    elif not(nonsing) and aimcode == 3:
        aimcode = 35
    elif not(nonsing) and aimcode == 4:
        aimcode = 45

    return (b,rts,ia,nexact,nnumeric,lgroots,aimcode)
