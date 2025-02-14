"""Principal Convex Hull Analysis (PCHA) / Archetypal Analysis."""

# Modified from the original implementation by:
# Authors: Morten Mørup, Technical University of Denmark (Matlab)
# Converted to Python by: Ulf Aslak, Technical University of Denmark
### Orginiated from the paper:
# Mørup, Morten, and Lars Kai Hansen. "Archetypal analysis for machine learning and data mining.".

# The original implementation can be found at:
#  https://mortenmorup.dk/?page_id=2 
# Orginal Python implementation:
# https://github.com/ulfaslak/py_pcha (Ulf Aslak)

# Note that this implementation has been modified significantly from the original implementation and is not optimized for binary data.



import numpy as np
from datetime import datetime as dt
import time



def PCHABer(X, C,S, noc, I=None, U=None, delta=0, verbose=False, conv_crit=1E-6, maxiter=500):
    """Return archetypes of dataset.

    Note: Commonly data is formatted to have shape (examples, dimensions).
    This function takes input and returns output of the transposed shape,
    (dimensions, examples).

    Parameters
    ----------
    X : numpy.2darray
        Data matrix in which to find archetypes

    noc : int
        Number of archetypes to find

    I : 1d-array
        Entries of X to use for dictionary in C (optional)

    U : 1d-array
        Entries of X to model in S (optional)


    Output
    ------
    XC : numpy.2darray
        I x noc feature matrix (i.e. XC=X[:,I]*C forming the archetypes)

    S : numpy.2darray
        noc x length(U) matrix, S>=0 |S_j|_1=1

    C : numpy.2darray
        noc x length(U) matrix, S>=0 |S_j|_1=1

    Ber : float
        Sum of Squared Errors


    """
    def S_update(S,PC,X,R,Ber, muS, niter):
        """Update S for one iteration of the algorithm."""
        noc, J = S.shape
        e = np.ones((noc, 1))
        for k in range(niter):
            Ber_old = Ber
            g = PC.T@(- X/R + (1-X)/(1-R))
            g = g - e * np.sum(g.A * S.A, axis=0)


            S_old = S
            while True:
                S = (S_old - g * muS).clip(min=0)
                S = S / np.dot(e, np.sum(S, axis=0))
                R = PC@S
                Ber = - (np.sum(np.multiply(X,np.log(R)) + np.multiply((1-X),np.log(1-R))))
                if Ber <= Ber_old * (1 + np.sign(Ber_old)*1e-9):
                    muS = muS * 1.2

                    break
                else:
                    muS = muS / 2

        return S,R, Ber, muS

    def C_update(X,P,R,S, C,Ber, delta, muC, niter=1):
        """Update C for one iteration of the algorithm."""
        J, nos = C.shape

        if delta != 0:
            alphaC = np.sum(C, axis=0).A[0]
            C = np.dot(C, np.diag(1 / alphaC))

        e = np.ones((J, 1))
        XtXSt = np.dot(X.T, XSt)

        for k in range(niter):

            # Update C
            Ber_old = Ber
            g = P.T@(-X/R + (1-X)/(1-R))@S.T

            g = g.A - e * np.sum(g.A * C.A, axis=0)

            C_old = C
            while True:
                C = (C_old - muC * g).clip(min=0)
                nC = np.sum(C, axis=0) + np.finfo(float).eps
                C = np.dot(C, np.diag(1 / nC.A[0]))


                PC = np.dot(P, C)
                R = PC@S
                Ber = - (np.sum(np.multiply(X,np.log(R)) + np.multiply((1-X),np.log(1-R))))

                if Ber <= Ber_old * (1 + np.sign(Ber_old)*1e-9):
                    muC = muC * 1.2
                    break
                else:
                    muC = muC / 2

        return C, Ber,R, muC, PC

    N, M = X.shape
    

    if I is None:
        I = range(M)
    if U is None:
        U = range(M)

    eps = 1e-3
    P = X+eps-2*eps*X
    PC = P@C
    R = PC@S

    # Initialize C

    muS, muC = 1, 1

    # Initialise S
    
    Ber = - (np.sum(np.multiply(X,np.log(R)) + np.multiply((1-X),np.log(1-R))))
    S,R, Ber, muS = S_update(S,PC,X,R,Ber, muS, 25)

    # Set PCHA parameters
    iter_ = 0
    dBer = np.inf
    t1 = dt.now()
    Ber_list = []

    if verbose:
        print('\nPrincipal Convex Hull Analysis / Archetypal Analysis')
        print('A ' + str(noc) + ' component model will be fitted')
        print('To stop algorithm press control C\n')

    dheader = '%10s | %10s | %10s | %10s | %10s | %10s' % ('Iteration', 'Cost func.', 'Delta Berf.', 'muC', 'muS', ' Time(s)   ')
    dline = '-----------+------------+------------+-------------+------------+------------+------------+------------+'

    while np.abs(dBer) >= conv_crit * np.abs(Ber) and iter_ < maxiter:
        if verbose and iter_ % 100 == 0:
            print(dline)
            print(dheader)
            print(dline)
        told = t1
        iter_ += 1
        Ber_old = Ber

        # C (and alpha) update
        XSt = np.dot(X[:, U], S.T)
        C, Ber,R, muC, PC = C_update(X,P,R,S, C,Ber, delta, muC, 10)
        # S update
        S,R, Ber, muS = S_update(S,PC,X,R,Ber, muS, 10)
        Ber_list.append(Ber)


        # Evaluate and display iteration
        dBer = Ber_old - Ber
        t1 = dt.now()
        if iter_ % 1 == 0:
            time.sleep(0.000001)
            if verbose:
                print('%10.0f | %10.4e | %10.4e | %10.4e | %10.4e | %10.4f \n' % (iter_, Ber, dBer/np.abs(Ber), muC, muS, (t1-told).seconds))

    # Display final iteration

    if verbose:
        print(dline)
        print(dline)
        print('%10.0f | %10.4f | %10.4e | %10.4e | %10.4e | %10.4f \n' % (iter_, Ber, dBer/np.abs(Ber), muC, muS, (t1-told).seconds))

    # Sort components according to importance
    ind, vals = zip(
        *sorted(enumerate(np.sum(S, axis=1)), key=lambda x: x[0], reverse=1)
    )
    S = S[ind, :]
    C = C[:, ind]
    PC = PC[:, ind]

    
    return PC, S, C, Ber_list