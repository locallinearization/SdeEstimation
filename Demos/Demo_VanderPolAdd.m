close all
disp('   ')
disp('SDE: Van der Pool oscillator with random force')
disp('  dx1 = x2*dt')
disp('  dx2 = (alpha-(x1^2-1)*x2-x1)*dt + sigma*dw')
disp('Parameters')
disp('  alpha=0.5,  sigma=0.75')
disp('Data')
disp('  31 equidistant noisy observations of x1 on [0,30]')
disp('Estimation method')
disp('  Innovation estimator with adaptive determinist LL filter')
disp('Optimization method')
disp('  fmincon')

Like_VanderPolAdd_Demo
