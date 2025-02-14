import numpy as np 

import itertools
import torch

from src.methods.fnnls import fastnnls
from src.methods.AALS import AALS





def SUpdate(X,PC,S,SSt,n_arc,n_samples,device):

    """
    Update function for the S parameter in archetypal analysis.

    This function updates the S parameter in archerypal analysis using SMO (sequential minimal optimization) updates for faster convergence. 

    The method if developed for bernoulli and Poisson distributed data, but can be extended to other distributions.
    Both torch and numpy are supported, but the format needs to be specified in the base argument.

    Args:
        X (ndarray): Data matrix.
        P (ndarray): Numerically stable version of the data matrix. Implemented as P = X+eps-2*eps*X
        PC (ndarray): Matrix product of the P and C parameters.
        S (ndarray): Matrix the data in archetypal space.
        SSt (ndarray): S@S.T
        n_arc (int): Number of archetypes.
        n_samples (int): Number of samples.
        method (str): Method for the loss function. Choose between 'Poisson' and 'Bernoulli'.
        base (str): Choose between 'numpy' and 'torch'.
    
    """

    alpha = torch.zeros(n_samples,dtype = torch.double, device = device)

    tol = 1e-12

    # Create a grid of all possible combinations of archetypes
    grid = torch.tensor(list(itertools.combinations(range(n_arc),2)))

    grid = torch.vstack((grid,grid))
    grid = grid[torch.randperm(grid.shape[0]),:]

    for j in range(grid.shape[0]):

        SS = S[grid[j],:].sum(axis=0)
        S0 = S[grid[j],:]

        K =  -(PC[:,grid[j]]).T@(X/(PC@S)) + (PC[:,grid[j]]).T@((1-X)/(1-(PC@S)))

        QuadS = torch.einsum('ji,jk,jl->ikl',PC[:,grid[j]],PC[:,grid[j]],X/((PC@S)**2)+(1-X)/((1-(PC@S))**2))
        
        linS = K - 2*torch.einsum('kl,ikl->il',S0,QuadS)
    

        denominator = 2*SS[SS>tol]*(QuadS[0,0,:][SS>tol] - QuadS[0,1,:][SS>tol] - QuadS[1,0,:][SS>tol] + QuadS[1,1,:][SS>tol])
        nominator = SS[SS>tol] *(QuadS[0,1,:][SS>tol] + QuadS[1,0,:][SS>tol] - 2*QuadS[1,1,:][SS>tol]) + linS[0][SS>tol] - linS[1][SS>tol]

        # Compute the alpha values.
        alpha[SS>tol]  = - torch.divide(nominator,denominator)

        alpha[alpha<0] = 0
        alpha[alpha>1] = 1

        S[grid[j,0],SS>tol] = SS[SS>tol] * alpha[SS>tol]
        S[grid[j,1],SS>tol] = SS[SS>tol] * (1-alpha[SS>tol])

        
        SSt[grid[j],:] = S[grid[j],:]@S.T
        SSt[:,grid[j]] = SSt[grid[j],:].T


    return S,SSt


def CUpdate(X,P,C,S,PC,n_arc,PPP,device):

    """

    Update function for the C parameter in archetypal analysis using fnnls and an active set strategy.

    This function updates the C parameter in archerypal analysis using fnnls and an active set strategy for faster convergence.
    This utilises the sparsity of the C matrix, as there are only few samples in each archetype.

    The method if developed for bernoulli data, but can be extended to other distributions.
    
    Args:
        X (ndarray): Data matrix.
        P (ndarray): Numerically stable version of the data matrix. Implemented as P = X+eps-2*eps*X
        C (ndarray): Matrix to be updated
        S (ndarray): Representation of the data in archetypal space.
        PC (ndarray): Matrix product of the P and C parameters.
        n_arc (int): Number of archetypes.
        PPP (list): List of active set indices.
    """

    for i in range(n_arc):

        v = ((X/((P@C@S)**2))@S[i,:]**2) + ((((1-X)/((1-(P@C@S))**2))@S[i,:]**2))

        grad = (P.T@((X/(P@C@S))@S[i,:]) -P.T@(((1-X)/(1-P@C@S))@S[i,:]))

        xtilde  = (P*torch.sqrt(v[:,None]))
        Xty =  grad + xtilde.T@(xtilde@C[:,i])

        tol = 1e-9
        C_k,_,PPP[i] = fastnnls(xtilde,Xty,tol,C[:,i],PPP[i])

        C[:,i] = C_k
        PC[:,i] = P@C_k

    return C,PC,PPP



def LossFunc(X,PC,S,device):

    """
    This function computes the loss function for archetypal analysis and returns the loss value.

    The method if developed for bernoulli distributed data, but can be extended to other distributions. 
   
    Args:
        X (ndarray): Data matrix.
        PC (ndarray): Matrix product of the P and C parameters.
        S (ndarray): Matrix of archetypes.


    Returns:
        loss (float): Loss value.

    """

    loss = torch.multiply(-X,torch.log(PC@S)).sum().sum() - torch.multiply((1-X),torch.log(1-(PC@S))).sum().sum() 
        

    loss = loss.item()


    return loss


def Bernoulli_Archetypal_Analysis(X,n_arc, maxIter = 500, convCriteria = 1e-6,  C = None, S = None,init = False):

    """
    An iterative update function for archetypal analysis. 

    This function calls the SUpdate and CUpdate functions to update the S and C parameters in an iterative manner.
    Both numpy and torch are supported, but the format needs to be specified in the base argument.

    The method if developed for bernoulli and Poisson distributed data, but can be extended to other distributions.


    Future work:
        Implement method for other distributions.

    Args:
        X (ndarray): Data matrix.
        idxC (ndarray): Relevant indices for X to form the archetypes
        n_arc (int): Number of archetypes.
        maxIter (int): Maximum number of iterations.
        convCriteria (float): Convergence criteria for the stopping criterion.
        C (ndarray): Initial guess for the C parameter.
        S (ndarray): Initial guess for the S parameter.

    Returns:
        C (ndarray): Final C parameter.
        S (ndarray): Final S parameter.
        L (list): List of loss values.
    
    
    """

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    eps = 1e-3
    _ , n_samples = X.shape

    P = X+eps-2*eps*X
    P = P.to(device)

    iter = 0 
    loss = 9
    lossOld = 10
    L = []

    if C is None:


        C =torch.abs(torch.log(torch.rand(n_samples,n_arc,dtype=torch.double, device = device)))
        C = C / C.sum(axis=0)[None,:]
        
    if S is None:

        S = torch.abs(torch.log(torch.rand(n_arc,n_samples,dtype=torch.double, device = device)))
        S = S / S.sum(axis=0)[None,:]

    if init:
        C,S,_,_ = AALS(X,n_samples,n_arc,C,S,gridS = True,maxIter = 50)

    PC = P@C
    SSt = S@S.T
    PPP = [[] for i in range(n_arc)]

    while (np.abs(loss - lossOld)/np.abs(loss) > convCriteria and iter<maxIter):
        
        
        #np.abs(loss - lossOld)/loss < eps
        lossOld = loss

        # Update S
        S,SSt = SUpdate(X,PC,S,SSt,n_arc,n_samples,device)

        # Update C
        C, PC,PPP = CUpdate(X,P,C,S,PC,n_arc,PPP,device) 

        # Update loss

        loss = LossFunc(X,PC,S,device)


        L.append(loss)

        iter += 1

    return C,S,L


    