function LBQ = LBQtest(res,maxLag,numPar,alpha)
% Ljung-Box Q-test for residual autocorrelation
%
% Description:
%
%   The "portmanteau" test of Ljung and Box assesses the null hypothesis
%   that a series of residuals exhibits no autocorrelation for a fixed
%   number of lags L, against the alternative that some autocorrelation
%   coefficient rho(k), k = 1, ..., L, is nonzero. The test statistic is
%
%                  L
%   	Q = T(T+2)Sum(rho(k)^2/(T-k)),
%                 k=1
%
%   where T is the sample size, L is the number of autocorrelation lags,
%   and rho(k) is the sample autocorrelation at lag k. Under the null, the
%   asymptotic distribution of Q is chi-square with L degrees of freedom.
%
% Input Arguments:
%
%   res - Vector of residuals for which the test statistic is computed. The
%   	last element corresponds to the most recent observation. Typically,
%       res contains (standardized) residuals obtained by fitting a model
%       to an observed time series.
%
% Output Arguments:
%
%   LQB.h - Vector of Boolean decisions for the tests, with length equal 
%   	to the number of tests. Values of h equal to 1 indicate rejection of the
%       null of no autocorrelation in favor of the alternative. Values of h
%       equal to 0 indicate a failure to reject the null.
%
%   LQB.pValue - Vector of p-values of the test statistics, with length equal
%       to the number of tests.
%
%   LQB.stat - Vector of test statistics, with length equal to the number of
%       tests.
%
%   LQB.cValue - Vector of critical values for the tests, determined by alpha,
%       with length equal to the number of tests.
%
% Notes:
%
%   o The input lags affects the power of the test. If L is too small, the
%     test will not detect high-order autocorrelations; if it is too large,
%     the test will lose power when significant correlation at one lag is
%     washed out by insignificant correlations at other lags. The default
%     value of min[20,T-1] is suggested by Box, Jenkins, and Reinsel [1].
%     Tsay [4] cites simulation evidence that a value approximating log(T)
%     provides better power performance.
%
%   o When res is obtained by fitting a model to data, the degrees of
%     freedom are reduced by the number of estimated coefficients,
%     excluding constants. For example, if res is obtained by fitting an
%     ARMA(p,q) model, dof should be L-p-q.
%
%   o LBQTEST does not test directly for serial dependencies other than
%     autocorrelation, but it can be used to identify conditional
%     heteroscedasticity (ARCH effects) by testing squared residuals. See,
%     e.g., McLeod and Li [3]. Engle's test, implemented by ARCHTEST, tests
%     for ARCH effects directly.
%
% Example:
%
%   % Test exchange rates for autocorrelation, ARCH effects:
%
%   load Data_MarkPound
%   returns = price2ret(Data);
%   residuals = returns-mean(returns);
%   h1 = lbqtest(residuals)
%   h2 = lbqtest(residuals.^2)
%
% References:
%
%   [1] Box, G.E.P., G.M. Jenkins, and G.C. Reinsel. Time Series Analysis:
%       Forecasting and Control. 3rd ed. Upper Saddle River, NJ:
%       Prentice-Hall, 1994.
% 
%   [2] Gourieroux, C. ARCH Models and Financial Applications. New York:
%       Springer-Verlag, 1997.
%
%   [3] McLeod, A.I. and W.K. Li. "Diagnostic Checking ARMA Time Series
%       Models Using Squared-Residual Autocorrelations." Journal of Time
%       Series Analysis. Vol. 4, 1983, pp. 269-273.
%
%   [4] Tsay,R.S. Analysis of Financial Time Series. Hoboken, NJ: John
%       Wiley & Sons, Inc., 2005.
%

T = length(res);       
lags = (1:maxLag)';
dof = max([lags-numPar ones(maxLag,1)],[],2);

% Compute the sample ACF out to the largest lag:
ACF = AutoCorr(res,maxLag);        
ACF = ACF(2:end);           % Strip off ACF at lag 0

% Compute Q-statistics to the largest lag; keep only those requested:
idx = (T-(1:maxLag))';
stat = T*(T+2)*cumsum((ACF.^2)./idx);
stat = stat(lags);

% Compute p-values:
pValue = 1-chi2cdf(stat,dof);

% Compute critical values 
cValue = chi2inv(1-alpha,dof);

% Perform the test:
h = (alpha >= pValue);

% output
LBQ.h=h;
LBQ.pValue=pValue;
LBQ.stat=stat;
LBQ.cValue=cValue;


function acf = AutoCorr(y,maxLag)
%Sample autocorrelation function
%
% Description:
%
%   Compute the sample autocorrelation function (ACF) of a univariate, 
%   stochastic time series y. 
%
% Reference:
%
%   [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%       Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%       NJ: Prentice-Hall, 1994.

y = y(:);           

% The ACF computation is based on [1], pages 30-34, 188:

nFFT = 2^(nextpow2(length(y))+1);
F = fft(y-mean(y),nFFT);
F = F.*conj(F);
acf = ifft(F);
acf = acf(1:(maxLag+1)); % Retain non-negative lags
acf = acf./acf(1); % Normalize
acf = real(acf);

