function stats = olswnwpanel(y,X,firmid,timeid,nlags)
% Multiple linear regression using ordinary least squares (OLS) with 
% Newey-West (1987) heteroskedasticity and autocorrelation consistent (HAC)
% covariance estimator modified for panel data sets in order to estimate 
% only correlations between lagged residuals in the same cluster (borrowed 
% from Gow, I.D., Ormazabal, G., & Taylor, D.J., 2010). It can be used for 
% a single time-series or panel data set. 
%
% stats = olswnwpanel(y,X,firmid,timeid,nlags) estimates OLS regression with 
% Newey-West HAC covariance estimator modified for panel data set. It 
% solves the OLS problem either through a QR matrix or Cholesky 
% decomposition, depending on the size of the data (borrowed from 
% LeSage, J.P., 2010).
%
% Inputs:
%   y       An n-by-1 vector of dependent variable observations
%   X       An n-by-k matrix of independent variables observations, with 
%           rows corresponding to observations and columns to independent 
%           variables
%	firmid  An n-by-1 vector for the firm identifier, which is the variable 
%           that denotes each firm (e.g., cusip, permno, or gvkey)
%   timeid  An n-by-1 vector for the time identifier, which is the variable 
%           that identifies the time dimension (e.g., year)
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
%           'tstat'         standard errors, t statistics, and p-values 
%                           based on Newey-West heteroskedasticity and 
%                           autocorrelation consistent covariance estimator 
%                           of regression coefficients
%           'fstat'         F statistic and p-value
%           'dwstat'        Durbin-Watson statistic 
% 
% Comments:
%   X should include a column of ones so that the model contains a constant
%   term. When using this function for a single time-series, firmid and 
%   timeid are n-by-1 vectors of ones. When this function is used for a 
%   panel data set (so that only observations within a cluster may be 
%   correlated), include firmid and timeid. This specification allows for 
%   observations on the same firm in different years to be correlated, 
%   i.e., a firm effect. If you want to allow for observations on different 
%   firms but in the same year to be correlated, reverse the firm and time 
%   identifiers. If you are clustering on some other dimension besides firm 
%   or time (e.g., industry or country), use that variable instead. You can
%   specify any lag length up to T - 1, where T is the number of time 
%   periods per firm. The R-square value is 1 minus the ratio of the error 
%   sum of squares to the total sum of squares. This value can be negative 
%   for models without a constant, which indicates that the model is not 
%   appropriate for the data. The F statistic and its p-value are 
%   calculated under the assumption that the model contains a constant term 
%   and they are not correct for models without a constant. 
%
% References:
%   Gow, I.D., Ormazabal, G., & Taylor, D.J. (2009). Correcting for 
%       Cross-Sectional and Time-Series Dependence in Accounting Research, 
%       The Accounting Review, 85, 483-512.
%   Greene, W.H. (2002). Econometric Analysis, Englewood Cliffs, 
%       NJ: Prentice-Hall.
%   LeSage, J.P. (2010) Econometrics Toolbox. Retrieved from 
%       http://www.spatial-econometrics.com/.
%   Newey, W. K. & West, K. D. (1987). A Simple, Positive Semi-Definite, 
%       Heteroskedasticity and Autocorrelation Consistent Covariance Matrix, 
%       Econometrica, 55, 703708.
%   Petersen, M.A. (2009). Estimating Standard Errors in Finance Panel Data 
%       Sets: Comparing Approaches, Review of Financial Studies, 22, 435-480.

[nobs, nvar] = size(X);

% Sort the data by firmid and then by timeid
firmidlist = unique(firmid);
i = 0;
for i = 1:length(firmidlist);
    rows = find(firmid == firmidlist(i));
    [tmp, idx] = sort(timeid(rows,:));
    rows = rows(idx,:);
    if i == 1
        y2 = y(rows,:);
        X2 = X(rows,:);
        F2 = firmid(rows,:);
        T2 = timeid(rows,:);
    else
        y2 = [y2; y(rows,:)];
        X2 = [X2; X(rows,:)];
        F2 = [F2; firmid(rows,:)];
        T2 = [T2; timeid(rows,:)];
    end
end

% Reassign the variables after sorting
y = y2;
X = X2;
firmid = F2;
timeid = T2;

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
r = y - stats.yhat;
stats.r = r; 
sse = stats.r' * stats.r;
stats.mse = sse / (nobs - nvar);
stats.covb = xpxi * stats.mse;
ydemean = y - mean(y);
stats.rsquare = 1 - sse / (ydemean' * ydemean);
adjrsquare1 = sse / (nobs - nvar);
adjrsquare2 = (ydemean' * ydemean) / (nobs - 1);
stats.adjrsquare = 1 - (adjrsquare1 / adjrsquare2);

stats.tstat.beta = stats.beta;
stats.tstat.se = sqrt(diag(stats.covb));

% Covariance matrix robust to Newey-West heteroskedasticity and 
% autocorrelation (see Greene, 2010)
if nargin < 4 && nlags <= 0 
    nlags = floor(4 * ((nobs / 100)^(2 / 9)));  % Use common lag length
end
Q = 0;
for l = 0:nlags
    wl = 1 - l / (nlags + 1);
    for t = l + 1:nobs
        if (l == 0)     % This calculates the S_0 portion
            Q = Q  + r(t)^2 * X(t, :)' * X(t,:);
        else            % This calculates the off-diagonal terms
            if firmid(t,1) == firmid(t - l,1)
                Q = Q + wl * r(t) * r(t - l) * ...
                    (X(t, :)' * X(t - l,:) + X(t - l, :)' * X(t,:));
            end
        end
    end
end
Q = 1 / (nobs - nvar) * Q;

stats.tstat.covb = nobs * xpxi * Q * xpxi;
stats.tstat.se = sqrt(diag(stats.tstat.covb));
stats.tstat.t = stats.beta ./ stats.tstat.se;
stats.tstat.pval = 2*(tcdf(-abs(stats.tstat.t),dfe));
stats.tstat.dfe = nobs - nvar;

stats.fstat.sse = sse;
stats.fstat.dfe = dfe;
stats.fstat.dfr = dfr;
stats.fstat.ssr = sum((stats.yhat - mean(stats.yhat)).^2);
stats.fstat.f = (stats.rsquare/(stats.fstat.dfr)) / ((1 - stats.rsquare)/stats.fstat.dfe);
stats.fstat.pval = 1 - fcdf(stats.fstat.f,stats.fstat.dfr,stats.fstat.dfe);

stats.loglike = -nobs/2 * (1 + log(2*pi) + log(stats.fstat.sse / nobs));
ediff = stats.r(2:nobs) - stats.r(1:nobs - 1);
stats.dwstat = (ediff' * ediff) / stats.fstat.sse;
end
