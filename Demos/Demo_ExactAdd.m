close all
disp('   ')
disp('SDE')
disp('  dx = alpha*t*x*dt + sigma1*t^2*exp(0.5*alpha*t^2)*dw1 + sigma2*sqrt(t)*dw2')
disp('Parameters')
disp('  alpha=0.25,  sigma1=5,   sigma2=0.1')
disp('Data')
disp('  11 equidistant noisy observations of x on [0.01,10.01]')
disp('Estimation method')
disp('  Exact Innovation estimator')
disp('Optimization method')
disp('  fmincon')

Like_ExactAdd_Demo

disp('Warning: The fitting-innovation is not a Gaussian white noise')
disp('         The provided data are insufficient for the estimation')
disp('   ')
