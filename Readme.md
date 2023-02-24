# SdeEstimation
SdeEstimation toolbox is a set of procedures designed for the state and parameters estimation of stochastic differential equation from discrete observations. The toolbox provides various implementations of Local Linearization filters for the state estimation and, consequently, of the Innovation Estimators for the parameters. The users can choose the implementation that better fits their needs in correspondence with the  particular model and data under consideration, and of the required accuracy for the estimation. This includes deterministic and stochastic filters with fixed step sizes and number of samples, with adaptive time stepping algorithms, with adaptive sampling algorithms, as well as local and global optimization algorithms for computing the innovation estimators.  Parallel computations are automatically set when more than one core are available, plus Krylov-subspace approximations for computing the predictions when the model has more than six equations. Confidence intervals are provided for the estimated parameters as well as the AIC and BIC values of the fitted  model. For the estimated parameters, the filtering algorithms return the predictions and filters of the states, the innovation process, the values of the Kolmogorov-Smirnov, Jarque-Bera, Ljung-Box, and Engle test for the innovation process, and the log likelihood value. The provided measures of goodness of fit of the model to the data make the toolbox also useful for the model selection, for designing of models and for optimal experimental design. The toolbox is intended for the recurrent practical situation in which a diffusion process should be identified from a reduced number of possible noisy and partial observations of the state variables, distant in time and with missing data. A number of illustrative test models and demos are included. 

## SdeEstimation Toolbox: Instructions for the users
- [InstructionsForTheUsers.pdf](https://github.com/locallinearization/SdeEstimation/raw/main/InstructionsForTheUsers.pdf)

## Functions for parameter and state estimation
- InnovEstimator: computes innovation estimators of the unknown parameters and unobserved states 
- EstimSet: sets the options for the parameter estimation
- DisplayParameterEstimates: summarizes the results of the parameter estimation 
- PlotStatTestFitInn: plots the statistics of the standarized fitting-innovation

## Functions for filtering
- LLfilter: computes Local Linearization filters 
- PlotFilterEstimates: plots the observations and their filter estimates
- PlotModelEstimates: plots exact and approximate conditional moments

## Functions for model simulation
- Euler: computes the Euler-Maruyama scheme for SDEs 
- Observations: generates time series of noisy observations of the model state variables
- PlotModelRealization: plots the simulated state variables and their observations

## Demos
- [Demo_CIR](./Demos/Demo_CIR.m)
- [Demo_ExactAdd](./Demos/Demo_ExactAdd.m.m)
- [Demo_ExactMul](./Demos/Demo_ExactMul.m.m)
- [Demo_FitzhughNagumo](./Demos/Demo_FitzhughNagumo.m.m)
- [Demo_LinearOscillator](./Demos/Demo_LinearOscillator.m.m)
- [Demo_SIR](./Demos/Demo_SIR.m.m)
- [Demo_VanderPolAdd](./Demos/Demo_VanderPolAdd.m.m)
- [Demo_VanderPolMul](./Demos/Demo_VanderPolMul.m.m)

Toolbox's initialization paths
- SdeEstimationPath: Always run first this script to access all the toolbox functions

