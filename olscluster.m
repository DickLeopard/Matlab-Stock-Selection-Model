function stats = olscluster(y,X,g,varargin)
% Multiple linear regression using ordinary least squares (OLS) with 
% 1-way or 2-way cluster-robust covariance estimator (borrowed from 
% Gow, I.D., Ormazabal, G., & Taylor, D.J., 2010).
%
% stats = olscluster(y,X,g,varargin) calculates OLS regression with 1-way 
% or 2-way cluster-robust covariance estimator. It solves the OLS problem 
% either through a QR matrix decomposition or Cholesky decomposition, 
% depending on the size of the data (borrowed from LeSage, J.P., 2010).
%
% Inputs:
%  y        An n-by-1 vector of dependent variable observations
%  X        An n-by-k matrix of independent variables observations, with 
%           rows corresponding to observations and columns to independent 
%           variables
%  g        An n-by-1 vector of 1st cluster variable observations
%  varargin An n-by-1 vector of (optional) 2nd cluster variable
%           observations
% 
% Outputs:
%  stats    An output structure containing the following statistics:
%           'beta'          Regression coefficients
%           'yhat'          Fitted values of the dependent variable
%           'r'             Residuals
%           'mse'           Mean squared error
%           'covb'          Covariance matrix of regression coefficients
%           'rsquare'       R-square statistic
%           'adjrsquare'    Adjusted R-square statistic
%           'tstat'         standard errors, t statistics, and p-values 
%                           based on 1-way or 2-way cluster-robust 
%                           covariance estimator of regression coefficients
%           'fstat'         F statistic and p-value
%           'dwstat'        Durbin-Watson statistic 
% 
% Comments:
%   y, g, and varargin should be vectors of equal length. X should have the 
%   same number of observations as y and it should include a column of ones
%   so that the model contains a constant term. The R-square value is 1
%   minus the ratio of the error sum of squares to the total sum of 
%   squares. This value can be negative for models without a constant, 
%   which indicates that the model is not appropriate for the data. The F 
%   statistic and its p-value are calculated under the assumption that the 
%   model contains a constant term and they are not correct for models 
%   without a constant. 
%
% References:
%   Cameron, A.C., Gelbach, J.B., & Miller, D.L. (2011). Robust Inference
%       with Multi-way Clustering, Journal of Business and Economic 
%       Statistics, 29, 238-249.
%   Gow, I.D., Ormazabal, G., & Taylor, D.J. (2009). Correcting for 
%       Cross-Sectional and Time-Series Dependence in Accounting Research, 
%       The Accounting Review, 85, 483-512.
%   LeSage, J.P. (2010) Econometrics Toolbox. Retrieved from 
%       http://www.spatial-econometrics.com/.

[nobs, nvar] = size(X);

if nargin < 3 || nargin > 4
    error('Either 1 or 2 cluster variables are required');
end

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

% Covariance matrix robust to clustering on the 1st variable
covb = onecluster(xpxi, X, r, g);
% If 2 cluster variables are supplied, adjust the covariance matrix for 
% the 2nd cluster variable
if nargin == 4
    h = varargin{1};
    gh = [g h];
    % Add the component attributable to the 2nd cluster
    tmp = covb + onecluster(xpxi, X, r, h);
    % Subtract the variance attributable to the intersection cluster (see 
    % Cameron, A.C., Gelbach, J.B., & Miller, D.L., 2011)
    covb = tmp - onecluster(xpxi, X, r, gh);
end
    
stats.tstat.covb = covb;
stats.tstat.se = sqrt(diag(covb));
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

% Nested function that handles a single cluster
    function varb = onecluster(xpxi, X, r, g)
        % Identify unique clusters
        G = unique(g, 'rows');
        M = size(G,1);
        % Calculate the central matrix for the sandwich variance estimator 
        % (see Cameron, A.C., Gelbach, J.B., & Miller, D.L., 2011)
        mid = 0;
        % This code handles clusters defined by more than one variable, 
        % treating observations having the same value for all cluster 
        % variables as members of the same cluster
        for i = 1:size(G,1)
            test = [1:size(g,1)]';
            for j = 1:size(G,2)
                test2 = find(g(:,j) == G(i,j));
                test = intersect(test,test2);
            end
            Xg = X(test,:);
            rg = r(test,:);
            mid = mid + Xg' * rg * rg' * Xg;
        end;
        % Calculate the cluster-robust covariance matrix estimate
        qc = (nobs - 1) / (nobs - nvar) * M / (M - 1);
        varb = qc * xpxi * mid * xpxi;
    end
end
