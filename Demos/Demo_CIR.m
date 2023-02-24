close all
disp('   ')
disp('SDE: CIR model, finance')
disp('  dx = (alpha+beta*x)*dt + sigma*sqrt(x)*dw')
disp('Parameters')
disp('  alpha=0.0189,  beta=-0.2339,  sigma=0.0073')
disp('Data')
disp('  307 equidistant noisy observations of x on [0.01,25.398]')
disp('Estimation method')
disp('  Exact innovation estimator')
disp('Optimization method')
disp('  UMDA+fmincon')

Like_CIR_Demo

disp('Warning: The fitting-innovation is not a white noise')
disp('         The provided data are insufficient for the estimation')
disp('   ')

