import numpy as np
import torch
import torch.distributions as d

class DataGenerator():
    """
    Class that generates fake data.
    Implemented from Probabilistic Archetypal Analysis by Seth et al as described in section 5.

    Args:
        n_samples: int, number of subjects
        n_features: int, number of features
        n_arc: int, number of archetypes
        distribution: str, distribution of the data
        sparsity: float [0,1], low sparsity --> more sparse data
        alpha: float [0,1], low alpha --> more data around the archetypes
        seed: int, random seed
        noise: float [0,1], low noise --> more noise
        base: str, base of the data
    """

    def __init__(self, n_samples, n_features, n_arc, distribution, sparsity, alpha, seed, noise,base):
        self.n_samples = n_samples
        self.n_features = n_features
        self.n_arc = n_arc
        self.distribution = distribution
        self.sparsity = sparsity
        self.alpha = alpha
        self.seed = seed
        self.noise = noise
        self.base = base

    def BernoulliData(self):
        """
        Function that generates Bernoulli distributed data.
        Implemented from Probabilistic Archetypal Analysis by Seth et al as described in section 5.

        """

        eta = torch.ones([self.n_features,self.n_arc]) 
        x = torch.ones([self.n_features,self.n_samples])
        h = torch.ones([self.n_arc,self.n_samples])
        tmp = torch.ones([self.n_features,self.n_samples])

        while eta.shape[1]!=torch.unique(eta,dim=1).shape[1]:
            eta = d.Bernoulli(self.sparsity).sample([self.n_features,self.n_arc])

        for j in range(self.n_samples):
            ind = torch.randint(self.n_features, (1,)) 
            #exp = Exponential(torch.ones([n_arc])*(alpha)).sample()
            h[:,j] = d.Dirichlet(torch.ones([self.n_arc])*(self.alpha)).sample()
            tmp[:,j]=(eta@h[:,j]).clip(max=0.99)
            x[:,j] = d.Bernoulli(tmp[:,j]).sample()


        #x = (x+1*10**(-6)).clip(min=0)

        return x

    def PoissonData(self):
        Bernoulli = self.BernoulliData()

        data = d.Poisson(5*Bernoulli).sample()

        return data
    
    def GaussianData(self):
        eta = torch.ones([self.n_features,self.n_arc]) 
        x = torch.ones([self.n_features,self.n_samples])
        h = torch.ones([self.n_arc,self.n_samples])
        tmp = torch.ones([self.n_features,self.n_samples])

        while eta.shape[1]!=torch.unique(eta,dim=1).shape[1]:
            eta = d.Bernoulli(self.sparsity).sample([self.n_features,self.n_arc])

        for j in range(self.n_samples):
            ind = torch.randint(self.n_features, (1,)) 
            #exp = Exponential(torch.ones([n_arc])*(alpha)).sample()
            h[:,j] = d.Dirichlet(torch.ones([self.n_arc])*(self.alpha)).sample()
            tmp[:,j]=(eta@h[:,j]) #.clip(max=0.99)
            data = tmp + d.Normal(0,self.noise).sample([self.n_features,self.n_samples])



        return data
    
    def CreateData(self):

        if self.distribution == "Bernoulli":
            data=self.BernoulliData()
        elif self.distribution == "Poisson":
            data = self.PoissonData()
        elif self.distribution == "Gaussian":
            data= self.GaussianData()
        else:
            raise ValueError("Distribution not implemented")
        
        if self.base == "numpy":
            data = data.numpy()

        return data


