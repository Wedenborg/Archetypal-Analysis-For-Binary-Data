import torch
import itertools
import numpy as np
from src.methods.fnnls import fastnnls
torch.set_printoptions(precision=8)


torch.set_default_dtype(torch.double)
## S update!! 
def S_updateTorch(S,CtXtXC,XCtX,SSt,n_arc, n_samples,gridS = False,device = 'cpu'):
    """
    Updates the S parameter in archetypal analysis using sequential minimal optimization.

    Args:
        :S (ndarray): Matrix of archetypes.
        :type S: ndarray
        CtXtXC (ndarray): Gram matrix of the transposed input data X.
        XCtX (ndarray): Cross-correlation matrix between the input data X and the transposed input data X.
        CtXtX (ndarray): Cross-correlation matrix between the transposed input data X and the target variable y.
        SST (float): Total sum of squares of the target variable y.
        SSt (ndarray): Sum of squares of the archetypes.
        n_arc (int): Number of archetypes.
        n_samples (int): Number of samples in the input data X.

    Returns:
        S (ndarray): Updated matrix of archetypes.
        loss (list): List of loss values.
        SSt (ndarray): Updated S@S.T matrix.
    """    
    #alpha = torch.zeros(n_samples,device = device)


    tol = torch.tensor(1e-3,dtype=torch.double,device = device)
    eta = torch.tensor(1e-16,dtype=torch.double,device = device)

    # Create a grid of all possible pairs of archetypes.
    grid = torch.tensor(list(itertools.combinations(range(n_arc),2)))
    grid = grid[torch.randperm(grid.shape[0]),:]

    if gridS:
        grid = np.vstack((grid,grid))

    alpha = torch.zeros(n_samples,dtype=torch.double,device = device)

    for j in range(grid.shape[0]):
        SS = S[grid[j],:].sum(axis=0)
        S0 = S[grid[j],:]

        

        hess = CtXtXC[grid[j],:][:,grid[j]] 

        h_lin = CtXtXC[grid[j],:]
        
        d1 = -2*XCtX[grid[j],:]+2*(h_lin@S)
        d2 = 2*(hess@S0)
        d = d1-d2


        
        denominator = 2*SS[SS>tol]*(hess[0,0]-hess[0,1]-hess[1,0]+hess[1,1])

        nominator = SS[SS>tol] *(hess[0,1]+hess[1,0]-2*hess[1,1]) + d[0][SS>tol] - d[1][SS>tol]

        # Compute the alpha values.
        alpha[SS>tol]  = - torch.divide(nominator,denominator+eta)

        # Clip the values of alpha to the interval [0,1].
        alpha[alpha<0] = 0
        alpha[alpha>1] = 1

        # Update the archetypes using the computed alpha values.
        S[grid[j,0],SS>tol] = SS[SS>tol] * alpha[SS>tol]
        S[grid[j,1],SS>tol] = SS[SS>tol] * (1-alpha[SS>tol])

        
        SSt[grid[j],:] = S[grid[j],:]@S.T
        SSt[:,grid[j]] = SSt[grid[j],:].T


    return S,SSt





## C update!!

def C_updateTorch(X,C,S,SSt,XC,CtXtXC,n_arc,PP,device='cpu'):
    """
    Updates the C parameter in archetypal analysis using non-negative least squares.

    Args:
        C (ndarray): 
        S (ndarray): 
        SSt (ndarray): S@S.T
        XC (ndarray): X@C
        CtXtXC (ndarray): C.T@X.T@X@C
        SST (float): Total sum of squares of the full data matrix X.
        XSt (ndarray): X@S.T .
        n_arc (int): Number of archetypes.
        XtX (ndarray): Quadratic matrix.

    Returns:
        C (ndarray): Updated matrix.
        loss (list): List of loss values.
        CtXtXC (ndarray): Updated Gram matrix of the transposed data XC.
    """
    loss = []
    L = []


    for i in range(n_arc):


        Xtilde = torch.multiply(X,torch.sqrt(SSt[i,i]))

        Xty =  2*(X.T@(X@S[i,:])-X.T@(X@(C@(S@S[i,:])))) +Xtilde.T@(Xtilde@C[:,i])

        C_k,_,PP[i] = fastnnls(Xtilde,Xty,1e-9,C[:,i],PP[i],device)


        C[:,i] = C_k
        XC[:,i] = X@C_k

        CtXtXC[i,:] = XC[:,i]@XC
        CtXtXC[:,i] = CtXtXC[i,:]



    return C,XC, CtXtXC,PP


## AA function!!


def AALS(X,n_samples,n_arc,C=None,S=None,gridS = False,maxIter = 1000,device='cpu'):


    if C is None:
        C =torch.log(torch.rand(n_samples,n_arc,dtype=torch.double,device=device,requires_grad=False))
        C[torch.abs(C)<0.5] = 0
        C = C / C.sum(axis=0)[None,:]

    if S is None:
        S = torch.log(torch.rand(n_arc,n_samples,dtype=torch.double,device=device,requires_grad=False))
        S = S / S.sum(axis=0)[None,:] 

    L = []
    EV = []
    iter = 0 
    loss = torch.tensor(0)
    loss_old = torch.tensor(1e10)
    SSt = S@S.T
    

    SST = (X**2).sum()
    
    XC = X@C

    CtXtXC = XC.T@XC
    XCtX = XC.T@X

    PPP = [[] for i in range(n_arc)]

    #torch.random.shuffle(grid) # Do I really need to shuffle the grid?

    while (torch.abs(loss_old-loss)>1e-6*torch.abs(loss_old) and iter<maxIter):
        loss_old = loss

        
        S,SSt = S_updateTorch(S,CtXtXC,XCtX,SSt,n_arc, n_samples,gridS = False,device = device)

        

        C ,XC, CtXtXC,PPP = C_updateTorch(X,C,S,SSt,XC,CtXtXC,n_arc,PPP,device = device)

        XC = X@C
        XCtX = XC.T@X


        loss = SST-2*torch.sum(XCtX*S)+torch.sum(CtXtXC*SSt)
 
        L.append(loss.item())

        varExpl = (SST-loss)/SST

        EV.append(varExpl.item())

        iter = iter + 1
    return C,S,L, EV