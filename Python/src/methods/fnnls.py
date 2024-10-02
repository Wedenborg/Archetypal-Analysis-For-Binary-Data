
## FNNLS update!!
import torch
import numpy as np
torch.set_printoptions(precision=8)

def fastnnls(X,Xty,tol,b,PP,device='cpu'):
    
    """
    Non-negative least squares solver using the fast NNLS algorithm.

    Args:
        XtX (ndarray): Gram matrix of the input data X.
        Xty (ndarray): Cross-correlation vector between the input data X and the target variable y.
        tol (float): Tolerance for the stopping criterion.
        b (ndarray): Initial guess for the solution.

    Returns:
        x (ndarray): Solution to the non-negative least squares problem.
        P (ndarray): Set of active indices in the solution.
    """


    #device = "cpu"
    n = X.shape[1] # Number of columns in the input data X.
    # Tag det gamle aktive sæt, som er de værdier, der ikke er 0. 

    if not len(PP):
        P = (torch.empty((n,1), dtype=torch.long,device=device)*np.nan).flatten()
        PP = torch.argwhere(~torch.isnan(P)).flatten()

        Z = torch.arange(1,n+1,dtype=torch.float,device=device) # Initialize the set of candidate indices.
        ZZ = (Z-1).type(torch.long)

    else:
        P = (torch.empty((n,1), dtype=torch.float,device=device)*np.nan).flatten()
        P[PP] = PP.float()

        Z = torch.arange(1,n+1,dtype=torch.float,device=device) # Initialize the set of candidate indices.
        Z[PP] = np.nan
        ZZ = torch.argwhere(torch.isnan(P)).flatten()


    
    x = b.clone().type(torch.double)   #torch.zeros(n,dtype = torch.double,device=device) # Initialize the solution vector.
    z = b.clone().type(torch.double)   # Initialize the search direction.
   
    ij = None
    iterOuter = 0

    scale = torch.mean(torch.sum(X*X,axis =1))
    lambda1 = 1e9*scale
    lambda2 = 1e-6*scale
    t = X[:,PP]@x[PP]
    w = Xty +lambda1 - X.T@t - torch.sum(x)*lambda1 - lambda2*x # Compute the gradient of the objective function.
    iter = 0 
    # ændre tol til tol/np.mean(np.diag(XtX))

    if len(PP) == 0:
        XtX_PP = None 
    else:
        XtX_PP = X[:,PP].T@X[:,PP]

    
    while(len(ZZ)>0 and torch.sum(w[ZZ]>tol)>0 and iterOuter<1000):

        iter = 0
        temp = w[ZZ] # Select the active gradient components.
        t = torch.argmax(temp)  # Select the index with the largest gradient.
        t = ZZ[t] # Convert the index to the original numbering.
        #P[t] = t # Add the index to the set of active indices.
        #Z[t] = np.nan # Remove the index from the set of candidate indices.


        #PP = PP.append(t) #torch.argwhere(~torch.isnan(P)).flatten() # Select the active indices.
        PP = torch.cat((PP,torch.tensor([t],dtype=torch.long,device=device)),0)
        ZZ = ZZ[ZZ!=t]  #torch.argwhere(~torch.isnan(Z)).flatten() # Select the candidate indices.

        if XtX_PP is None:
            XtX_PP = X[:,PP].T@X[:,PP]
        else:
            #XtX_PP = XtX_PP.append(X[:,t].T@X[:,PP[:-1]],axis = 1)
            #XtX_PP = XtX_PP.append(X[:,PP].T@X[:,t],axis = 0)

            if len(PP[:-1]) == 1:
                XtX_PP = torch.hstack((XtX_PP[0],X[:,t].T@X[:,PP[:-1]]))
                XtX_PP = torch.vstack((XtX_PP,X[:,PP].T@X[:,t]))

            else:

                XtX_PP = torch.hstack((XtX_PP,(X[:,t].T@X[:,PP[:-1]]).unsqueeze(1)))
                XtX_PP = torch.vstack((XtX_PP,X[:,PP].T@X[:,t]))


        #lamb_inner =  torch.mean(torch.diag(XtX_PP)) 
        #lamb = lamb_inner*10e9 # ændret fra 1

        activeSetChange = True
        z = z*0
        #z[PP] = torch.linalg.solve(X[:,PP].T@X[:,PP]*SSt+lamb,Xty[PP]+lamb) # Compute the search direction.
        z[PP] = torch.linalg.solve(XtX_PP+lambda1+torch.eye(len(PP))*lambda2,Xty[PP]+lambda1) # Compute the search direction.


        while (torch.sum(z[PP]<=tol)>0 and iter<500): 


            temp = torch.argwhere(z[PP] <= tol).T[0] # Select the indices with small search direction.

            QQ =PP[temp] # ]temp[(temp.view(1, -1) == temp2.view(-1, 1)).any(dim=0)]


            #QQ = np.intersect1d(temp,temp2) # Select the indices with small search direction and active indices.

            a = x[QQ]/((x[QQ]-z[QQ]))
            a[torch.isnan(a)] = np.inf
            alpha = torch.min(a) # Compute the step size.


            x = x + alpha*(z-x) # Update the solution vector.


            t1 = torch.argwhere(x[PP]<tol).flatten() # Select the indices with small solution.
            ij = PP[t1] # Select the indices with small solution and active indices.

            idx_PP = torch.ones_like(PP).bool() if len(ij)==0 else ~torch.eq(PP[:, None], ij).T[0] 

            PP = PP[idx_PP]

            
            ZZ = torch.cat((ZZ,ij),0)



            if XtX_PP is None:
                XtX_PP = X[:,PP].T@X[:,PP]
            else:
                XtX_PP = XtX_PP[idx_PP][:,idx_PP]




            #lamb = torch.mean(torch.diag(XtX_PP))*10e9

            #z[PP] = torch.linalg.solve(X[:,PP].T@X[:,PP]*SSt+lamb,Xty[PP]+lamb) # Compute the search direction.
            z = z*0
            z[PP] = torch.linalg.solve(XtX_PP+lambda1+torch.eye(len(PP))*lambda2,Xty[PP]+lambda1) # Compute the search direction.



            iter = iter + 1

        
        x = z.clone()
        
        t = X[:,PP]@x[PP]
        w = Xty +lambda1 - X.T@t - torch.sum(x)*lambda1 - lambda2*x 
    
        iterOuter = iterOuter + 1

    
    out = z.clone()
        
    return out,w,PP
