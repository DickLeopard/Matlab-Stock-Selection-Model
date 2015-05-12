function stats = olswnw(y,X,robust,nlags)
% Multiple linear regression using ordinary least squares (OLS) with 
% either White (1980) heteroskedasticity consistent or Newey-West (1987) 
% heteroskedasticity and autocorrelation consistent (HAC) covariance 
% estimator. 
%
% stats = olswnw(y,X,robust,nlags) calculates OLS regression with either 
% White heteroskedastic consistent or Newey-West HAC covariance estimator. 
% It solves the OLS problem either through a QR matrix or Cholesky 
% decomposition, depending on the size of the data (borrowed from 
% LeSage, J.P., 2010).
%
% Inputs:
%   y       An n-by-1 vector of dependent variable observations
%   X       An n-by-k matrix of independent variables observations, with 
%           rows corresponding to observations and columns to independent 
%           variables
%   robust  Indicator for type of estimator used to obtain robust  
%           covariance matrix, select between: 
%           '0'     White (1980) heteroskedasticity consistent estimator
%           '1'     Newey-West (1987) heteroskedasticity and 
%                   autocorrelation consistent (HAC) estimator
%   nlags   Non-negative integer containing the lag length to be used in 
%           the HAC estimator, which represents the number of non-zero 
%           autocorrelations in the dependent data. When nlags is not 
%           supplied, nlags equals the common lag length as suggested by 
%           Newey-West (1987)
%  
% Outputs:
%   stats   An output structure containing the following statistics:
%           'beta'          Regression coefficients
%           'yhat'          Fitted values of the dependent variable
%           'r'             Residuals
%           'mse'           Mean squared error
%           'covb'          Covariance matrix of regression coefficients
%           'rsquare'       R-square statistic
%           'adjrsquare'    Adjusted R-square statistic
%           'tstat'         standard errors, t statistics and p-values 
%                           based on White or Newey-West robust covariance 
%                           estimator of regression coefficients
%           'fstat'         F statistic and p-value
%           'dwstat'        Durbin-Watson statistic 
% 
% Comments:
%   X should include a column of ones so that the model contains a constant
%   term. The R-square value is 1 minus the ratio of the error sum of 
%   squares to the total sum of squares. This value can be negative for 
%   models without a constant, which indicates that the model is not 
%   appropriate for the data. The F statistic and its p-value are 
%   calculated under the assumption that the model contains a constant 
%   term and they are not correct for models without a constant. 
%
% References:
%   LeSage, J.P. (2010) Econometrics Toolbox. Retrieved from 
%       http://www.spatial-econometrics.com/.
%   Newey, W. K. & West, K. D. (1987). A Simple, Positive Semi-Definite, 
%       Heteroskedasticity and Autocorrelation Consistent Covariance Matrix, 
%       Econometrica, 55, 703708.
%   White, H.L. (1980). A Heteroskedasticity-Consistent Covariance Matrix
%       Estimator and a Direct Test for Heteroskedasticity, Econometrica, 
%       48, 817-838.

[nobs, nvar] = size(X);

% Solution method for least squares problem in order to calculate (X'*X)^-1 
if nobs < 10000     % QR decomposition for problems with < than 10000 obs
    [q r] = qr(X,0);
    xpxi = (r'*r) \ eye(nvar);
else                % Cholesky decomposition for large problems
    xpxi = (X'*X) \ eye(nvar);
end

stats.beta = xpxi * (X'*y);
nvar = length(stats.beta);
dfe = nobs - nvar;
dft = nobs - 1;
dfr = nvar - 1;
stats.yhat = X * stats.beta;
stats.r = y - stats.yhat;
sse = stats.r' * stats.r;
stats.mse = sse / (nobs - nvar);
stats.covb = xpxi * stats.mse;
ydemean = y - mean(y); 
stats.rsquare = 1 -(sse / (ydemean' * ydemean));
adjrsquare1 = sse / (nobs - nvar);
adjrsquare2 = (ydemean' * ydemean) / (nobs - 1);
stats.adjrsquare = 1 - (adjrsquare1 / adjrsquare2);

stats.tstat.beta = stats.beta;
stats.tstat.se = sqrt(diag(stats.covb));
stats.tstat.t = stats.beta ./ stats.tstat.se;
stats.tstat.pval = 2*(tcdf(-abs(stats.tstat.t),dfe));
stats.tstat.dfe = nobs - nvar;

% Covariance matrix robust to
if robust == 0      % White heteroskedasticity 
    hhat = X .* stats.r(:,ones(1,nvar));
    xuux = hhat' * hhat;
    stats.tstat.covb = xpxi * (xuux) * xpxi;
    stats.tstat.se = sqrt(diag(stats.tstat.covb));
    stats.tstat.t = stats.beta ./ stats.tstat.se;
    stats.tstat.pval = 2*(tcdf(-abs(stats.tstat.t),dfe));
    stats.tstat.dfe = nobs - nvar;
elseif robust == 1  % Newey-West heteroskedasticity and autocorrelation
    if nargin < 4 && nlags <= 0 
        nlags = floor(4 * ((nobs / 100)^(2 / 9)));  % Use common lag length
    end 
    emat = [];
    for i = 1:nvar
        emat = [emat
            stats.r'];
    end;       
    hhat = X' .* emat;
    G = zeros(nvar,nvar); 
    w = zeros(2 * nlags + 1,1);
    j = 0;
    while j ~= nlags + 1
        ga = zeros(nvar,nvar);
        w(nlags + 1 + j,1) = (nlags + 1 - j) / (nlags + 1);
        za = hhat(:,(j + 1):nobs) * hhat(:,1 : nobs - j)';
        if j == 0
            ga = ga + za;
        else
            ga = ga + za + za';
        end
        G = G + w(nlags + 1 + j,1) * ga;
        j = j + 1;
    end;
    stats.tstat.covb = xpxi * G * xpxi;
    stats.tstat.se = sqrt(diag(stats.tstat.covb));
    stats.tstat.t = stats.beta ./ stats.tstat.se;
    stats.tstat.pval = 2*(tcdf(-abs(stats.tstat.t),dfe));
    stats.tstat.dfe = nobs - nvar;
end

stats.fstat.sse = sse;
stats.fstat.dfe = dfe;
stats.fstat.dfr = dfr;
stats.fstat.ssr = sum((stats.yhat - mean(stats.yhat)).^2);
stats.fstat.f = (stats.rsquare/(stats.fstat.dfr)) / ((1 - stats.rsquare)/stats.fstat.dfe);
stats.fstat.pval = 1-fcdf(stats.fstat.f,stats.fstat.dfr,stats.fstat.dfe);

stats.loglike = -nobs/2 * (1 + log(2*pi) + log(stats.fstat.sse / nobs));
ediff = stats.r(2:nobs) - stats.r(1:nobs - 1);
stats.dwstat = (ediff' * ediff) / stats.fstat.sse;
end
