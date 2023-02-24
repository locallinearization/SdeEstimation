close all
disp('   ')
disp('SDE')
disp('  dx = alpha*t*x*dt + sigma*sqrt(x)*dw')
disp('Parameters')
disp('  alpha=-0.1,  sigma=0.1')
disp('Data')
disp('  11 equidistant noisy observations of x on [0.5,10.5]')
disp('Estimation method')
disp('  Innovation estimator with adaptive determinist LL filter')
disp('Optimization method')
disp('  fmincon')

Like_ExactMul_Demo