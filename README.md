# Introduction
This repository contains the code introduced in the paper Archetypal Analysis for Binary Data, accepted at 2025 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)
If used please cite [Archetypal Analysis for Binary Data](https://arxiv.org/abs/2502.04172)

Archetypal analysis (AA) is a matrix decomposition method that identifies distinct patterns using convex combinations of the data points denoted archetypes with each data point in turn reconstructed as convex combinations of the archetypes. AA thereby forms a polytope representing trade-offs of the distinct aspects in the data. Most existing methods for AA are designed for continuous data and do not exploit the structure of the data distribution. In this paper, we propose two new optimization frameworks for archetypal analysis for binary data. i) A second order approximation of the AA likelihood based on the Bernoulli distribution with efficient closed-form updates using an active set procedure for learning the convex combinations defining the archetypes, and a sequential minimal optimization strategy for learning the observation specific reconstructions. ii) A Bernoulli likelihood based version of the principal convex hull analysis (PCHA) algorithm originally developed for least squares optimization. We compare these approaches with the only existing binary AA procedure relying on multiplicative updates and demonstrate their superiority on both synthetic and real binary data. Notably, the proposed optimization frameworks for AA can easily be extended to other data distributions providing generic efficient optimization frameworks for AA based on tailored likelihood functions reflecting the underlying data distribution. 

## :zap: Archetypal Analysis for Binary Data 
*Important!:* R is still not fully supported! Please use the python code instead!

The example notebook will take you through a small example with synthetic data for both Bernoulli and Gaussian data 
## :file_folder: File Structure
```
.
├── src
│   ├── methods
│   │   ├── __init__.py
│   │   └── AABernoulli.py
│   │   └── AALS.py
│   │   └── fnnls.py
│   │   └── PCHABer.py
│   ├── utils
│   │   ├── __init__.py
│   │   └── DataGenerator.py
|   |   └── NMI.py
│   └── visualizations
│       ├── plot_loss.py
│       └── PlotLossArc.py
│       └── PlotNMIStability.py
├── ExampleNB.ipynb
└── README.md
```
## :star2: Credit/Acknowledgment
This work is authored by Anna Emilie J. Wedenborg and Morten Mørup.


##  :lock: License
This work has been given a MIT license.
