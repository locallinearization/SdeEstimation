function [c, ceq] = ConFunLike(ThetaVar,EstimOptions)

global gct    % variable defined in LogLKFunction.m

% nonlinear inequality constraints
c=gct';

% nonlinear equality constraints
ceq = [];

