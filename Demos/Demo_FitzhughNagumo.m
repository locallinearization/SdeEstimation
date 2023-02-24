close all
disp('   ')
disp('SDE: stochastic Fitzhugh-Nagumo equation')
disp('  dx1 = alpha1*(x1-x1^3/3+x2)*dt')
disp('  dx2 = -(1/alpha1)*(x1-alpha2)*dt + alpha3*dw')
disp('Parameters')
disp('  alpha1=1,  alpha2=1,   alpha3=0.1')
disp('Data')
disp('  101 equidistant noisy observations of (x1,x2) on [0,50]')
disp('Estimation method')
disp('  Innovation estimator with adaptive determinist LL filter')
disp('Optimization method')
disp('  UMDA+fmincon')

Like_FitzhughNagumo_Demo
