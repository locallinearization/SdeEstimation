function ARCH = ARCHtest(res,varargin)
%ARCHTEST Engle test for residual heteroscedasticity
%
% Syntax:
%
%   [h,pValue,stat,cValue] = archtest(res)
%   [h,pValue,stat,cValue] = archtest(res,param1,val1,param2,val2,...)
%
% Description:
%
%   The ARCH test of Engle assesses the null hypothesis that a series of
%   residuals r(t) exhibits no conditional heteroscedasticity (ARCH
%   effects), against the alternative that an ARCH(L) model
%
%       r(t)^2 = a0 + a1*r(t-1)^2 + ... + aL*r(t-L)^2 + e(t)
%
%   describes the series with at least one nonzero ak for k = 0, ..., L.
%   ARCHTEST estimates the ARCH(L) model for a specified number of lags L
%   and then computes the Lagrange multiplier statistic T*R^2, where T is
%   the sample size and R^2 is the coefficient of determination of the
%   regression. Under the null, the asymptotic distribution of the test
%   statistic is chi-square with L degrees of freedom.
%
% Input Arguments:
%
%   res - Vector of residuals for which the test statistic is computed. The
%   	last element corresponds to the most recent observation. Typically,
%       res contains (standardized) residuals obtained by fitting a model
%       to an observed time series.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME        VALUE
%
%   'Lags'      Scalar or vector of positive integers indicating the number
%               of lags L used to compute the test statistic. Each element
%               must be less than length(res)-2. The default value is 1.
%
%   'Alpha'     Scalar or vector of nominal significance levels for the
%               tests. Elements must be greater than zero and less than
%               one. The default value is 0.05.
%
%   Scalar parameter values are expanded to the length of any vector value
%   (the number of tests). Vector values must have equal length. If any
%   value is a row vector, all outputs are row vectors.
%
% Output Arguments:
%
%   h - Vector of Boolean decisions for the tests, with length equal to the
%   	number of tests. Values of h equal to 1 indicate rejection of the
%       null of no ARCH effects in favor of the alternative. Values of h
%       equal to 0 indicate a failure to reject the null.
%
%   pValue - Vector of p-values of the test statistics, with length equal
%       to the number of tests.
%
%   stat - Vector of test statistics, with length equal to the number of
%       tests.
%
%   cValue - Vector of critical values for the tests, determined by Alpha,
%       with length equal to the number of tests.
%
% Notes:
%
%   o A suitable value for the number of lags must be determined in order
%     to draw valid inferences from the test. One method is to fit and 
%     obtain loglikelihoods for a sequence of progressively smaller ARCH(L) 
%     models, and then use LRATIOTEST to evaluate the significance of each 
%     restriction. Another method is to combine a measure of fit with the 
%     information criteria provided by AICBIC.
%
%   o Residuals in an ARCH process are dependent, but not correlated. Thus
%     ARCHTEST tests for heteroscedasticity without autocorrelation. To
%     test for autocorrelation, use LBQTEST.
%
%   o GARCH(P,Q) processes are locally equivalent to ARCH(P+Q) processes.
%     If ARCHTEST(res,'Lags',L) shows evidence of conditional 
%     heteroscedasticity in residuals from a mean model, then a GARCH(P,Q) 
%     model with P+Q = L may be appropriate.
%
% Example:
%
%   % Test exchange rates for ARCH effects:
%
%   load Data_MarkPound
%   returns = price2ret(Data);
%   residuals = returns-mean(returns);
%   h = archtest(residuals)
%
% References:
%
%   [1] Box, G.E.P., G.M. Jenkins, and G.C. Reinsel. Time Series Analysis:
%       Forecasting and Control. 3rd ed. Upper Saddle River, NJ:
%       Prentice-Hall, 1994.
%
%   [2] Engle, R. "Autoregressive Conditional Heteroscedasticity with
%       Estimates of the Variance of United Kingdom Inflation."
%       Econometrica. Vol. 96, 1988, pp. 893-920.
%
% See also LBQTEST, AICBIC.

% Copyright 2018 The MathWorks, Inc.   

% Parse inputs and set defaults. As of R2018a, this parse no longer supports 
% the deprecated syntax in which lags and alpha are input as the second and 
% third arguments, respectively. 

N = length(res); % Sample size
defaultLags = 1;
defaultAlpha = 0.05;

parseObj = inputParser;
parseObj.addRequired ('res'  ,               @(x) validateattributes(x, {'double'}, {'nonempty' 'vector'}, '', 'residual series'));
parseObj.addParameter('Lags' , defaultLags , @(x) validateattributes(x, {'double'}, {'nonempty' 'integer'  'vector' 'positive' '<=' (N-2)}, '', 'number of lags'));
parseObj.addParameter('Alpha', defaultAlpha, @(x) validateattributes(x, {'double'}, {'nonempty' 'vector' '>' 0 '<' 1}, '', 'significance levels'));

try
  parseObj.parse(res,varargin{:});
catch ME       
  throwAsCaller(ME)   
end

res   = parseObj.Results.res(:);
lags  = parseObj.Results.Lags;
alpha = parseObj.Results.Alpha;

% Check parameter values for commensurate lengths, expand scalars, and
% convert all variables to columns:

[numTests,rowOutput,lags,alpha] = sizeCheck(lags,alpha);
res2 = res.^2; % Squared residuals
   
% Perform the regression:

R2 = zeros(numTests,1);

for order = 1:numTests
    X         =  [ones(N,1),lagmatrix(res2,1:lags(order))];
    X         =  X(lags(order)+1:end,:);  % Regression design matrix
    y         =  res2(lags(order)+1:end); % Response variable
    yHat      =  X*(X\y);                 % Predicted response
    T         =  N-lags(order);           % Effective sample size
    yHat      =  yHat-sum(yHat)/T;        % De-meaned predicted response
    y         =  y-sum(y)/T;              % De-meaned observed response
    R2(order) =  (yHat'*yHat)/(y'*y);     % Centered R-squared
end

% Compute the test statistics:

stat = R2.*(N-lags);

% Compute p-values:

pValue = 1-chi2cdf(stat,lags);

% Compute critical values
  
cValue = chi2inv(1-alpha,lags);

% Perform the test:

h = (alpha >= pValue);

% Display outputs as row vectors if lags is a row vector:

if rowOutput
    
   h = h';
   pValue = pValue';
   stat = stat';

   if ~isempty(cValue)
       
      cValue = cValue';
      
   end
   
end

% output
ARCH.h=h;
ARCH.pValue=pValue;
ARCH.stat=stat;
ARCH.cValue=cValue;


%-------------------------------------------------------------------------
% Check parameter values for commensurate lengths, expand scalars, and
% convert all variables to columns
function [numTests,rowOutput,varargout] = sizeCheck(varargin)

% Initialize outputs:

numTests = 1;
rowOutput = false;

% Determine vector lengths, number of tests, row output flag:

for i = 1:nargin
        
    ivar = varargin{i};
    iname = inputname(i);
    
    paramLength.(iname) = length(ivar);
    numTests = max(numTests,paramLength.(iname));
    
    if ~isscalar(ivar)
        rowOutput = rowOutput || (size(ivar,1) == 1);
    end    
    
end

% Check for commensurate vector lengths:

for i = 1:(nargin-1)
    iname = inputname(i);
    for j = (i+1):nargin
        jname = inputname(j);
        if (paramLength.(iname) > 1) && (paramLength.(jname) > 1) ...
            && (paramLength.(iname) ~= paramLength.(jname))
        
            error(message('econ:archtest:ParameterSizeMismatch', iname, jname))
              
        end        
    end
end

% Expand scalars:

for i = 1:nargin
    
    ivar = varargin{i};
    if paramLength.(inputname(i)) == 1
        varargout{i} = ivar(ones(numTests,1)); %#ok
    else
        varargout{i} = ivar(:);  %#ok Column output
    end
    
end

end % sizeCheck

end % ARCHTEST