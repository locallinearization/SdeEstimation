close all
disp('   ')
disp('SDE: SIR model, epidemic')
disp('  dS = -beta*S*I*dt')
disp('  dI = (beta*S*I-gamma*I)*dt')
disp('  dR = gamma*I*dt')
disp('Parameters')
disp('  beta=0.0004,  gamma=0.04')
disp('Data')
disp('  100 equidistant noisy observations of (I,R) on [1,100] with random missing data')
disp('Estimation method')
disp('  Innovation estimator with adaptive determinist LL filter')
disp('Optimization method')
disp('  fmincon')

Like_SIR_Demo
