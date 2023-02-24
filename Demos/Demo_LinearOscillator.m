close all
disp('   ')
disp('SDE: Stochastic harmonic oscillator with random frequency and force')
disp('  dx1 = x2*dt')
disp('  dx2 = (alpha-omega*x1)*dt - sigma1*x1*dw1 + sigma1*dw2')
disp('Parameters')
disp('  alpha=10,  omega=30,   sigmal=0.5,   sigma2=0.5')
disp('Data')
disp('  801 equidistant noisy observations of x2 on [0,80]')
disp('Estimation method')
disp('  Exact innovation estimator')
disp('Optimization method')
disp('  fmincon')

Like_LinearOscillator_Demo
