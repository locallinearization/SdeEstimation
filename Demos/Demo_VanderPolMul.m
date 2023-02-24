close all
disp('   ')
disp('SDE: Van der Pool oscillator with random frequency')
disp('  dx1 = x2*dt')
disp('  dx2 = (-(x1^2-1)*x2-alpha*x1)*dt + sigma*x1*dw')
disp('Parameters')
disp('  alpha=1,  sigma=1')
disp('Data')
disp('  31 equidistant noisy observations of x1 on [0,30]')
disp('Estimation method')
disp('  Innovation estimator with adaptive determinist LL filter')
disp('Optimization method')
disp('  fmincon')

Like_VanderPolMul_Demo
