clear all; close all; clc

% 4.1 import data
exlServerObj = actxserver('Excel.Application');
exlWkbkObj = exlServerObj.Workbooks;
exlPath = 'D:\UMD\764 Quantitative Investment Strategy\Matlab Exercises';
exlFileList = 'lab4_data.xlsx';
exlFile = exlWkbkObj.Open([exlPath '\' exlFileList]);
exlSheet = exlFile.Sheets.Item('New Search Criteria');
exlRngObj = exlSheet.Range('A2:I8871');
exlData = exlRngObj.Value;
exlServerObj.Quit;
clear exlServerObj exlPath exlWkbkObj exlFileList exlFile exlSheet exlRngObj;
% exlData: 8870¡Á9 cell

% delete data which has missing values
rowsIndx = find(prod(double(~cellfun('isempty', exlData)), 2))'; %[1 2 3 ... 8870]
trimmedExlData = cell([numel(rowsIndx), size(exlData, 2)]); % 8870¡Á9 empty cells
count = 0;
for rowLoop = rowsIndx;  % rowLoop = 1 at the first time
    cellsWithNoTxt = ~cellfun('isclass', exlData(rowLoop, :), 'char'); % string:0, nonstring:1
    cellsWithNumsIndx = cellsWithNoTxt; 
    invalidNums = isnan([exlData{rowLoop, cellsWithNumsIndx}]); % at the first time, for the 1st row of the number's part, NaN:1, nonNaN:0
    if ~any(invalidNums);% at the 1st time, for the 1st row, all numbers with no NaN: 1, has NaN: 0
        count = count + 1;
        for columnLoop  = 1: size(exlData, 2);
            trimmedExlData{count, columnLoop} = exlData{rowLoop, columnLoop};
        end
    end
end
trimmedExlData = trimmedExlData([1: count], :);
clear rowsIndx count rowLoop cellsWithNoTxt cellsWithNumsIndx invalidNums 
clear columnLoop exlData % only trimmedExlData left
% trimmedExlData: 4986¡Á9 cell

% Create the corresponding variables
name = {trimmedExlData{:,1}}';
ticker = {trimmedExlData{:,2}}';
exchange = {trimmedExlData{:,3}}';
eP = [trimmedExlData{:,4}]';
clsgPrice = [trimmedExlData{:,5}]';
size = [trimmedExlData{:,6}]';
beMe = [trimmedExlData{:,7}]';
m12M = [trimmedExlData{:,8}]';
yr1Ret = [trimmedExlData{:,9}]';

stocks = dataset({nominal(name),'name'},...
    {nominal(ticker),'ticker'},{nominal(exchange),'exchange'},...
    {eP,'eP'},{clsgPrice,'clsgPrice'},{size,'size'},...
    {beMe,'beMe'},{m12M,'m12M'},{yr1Ret,'yr1Ret'});
% stocks: 4986¡Á9 dataset

% 4.2 Find  and delete outliers and extreme values
% Define outliers as three standard deviation from the mean
means = datasetfun(@mean,stocks,'DataVars',{'size','beMe','m12M','yr1Ret'}); 
stds = datasetfun(@std,stocks,'DataVars',{'size','beMe','m12M','yr1Ret'});
outliersIndx = abs(stocks.size - means(1)) > (3*stds(1)) |...
    abs(stocks.beMe - means(2)) > (3*stds(2)) |...
    abs(stocks.m12M - means(3)) > (3*stds(3)) |...
    abs(stocks.yr1Ret - means(4)) > (3*stds(4)) |...
    stocks.eP < 0 |...
    stocks.beMe < 0 |...
    stocks.clsgPrice < 5 |...
    stocks.size < 100000000;    
nOut = sum(outliersIndx); % nOut = 1986
stocks(outliersIndx,:) = []; % 3000 stocks left
% stocks: 3000¡Á9 dataset

% 4.3 lnSize
stocks.lnSize = log(stocks.size);
clear beMe clsgPrice eP exchange m12M means nOut name...
    outliersIndx size stds ticker trimmedExlData yr1Ret
% stocks: 3000¡Á10 dataset
% 4.4
% 4.4.1 scatter plot between explained variable and explanatory variables
figure, scatter(stocks.lnSize, stocks.yr1Ret, 'gx'), axis tight,  lsline;
title('scatter plot(lnsize-return)')
xlabel('lnSize'), ylabel('return')
figure, scatter(stocks.beMe, stocks.yr1Ret, 'gx'), axis tight,  lsline;
title('scatter plot(beMe-return)')
xlabel('beMe'), ylabel('return')
figure, scatter(stocks.m12M, stocks.yr1Ret, 'gx'), axis tight,  lsline;
title('scatter plot(m12M-return)')
xlabel('m12M'), ylabel('return')


figure, plotregression(stocks.lnSize, stocks.yr1Ret), axis auto;
axis([18 25 -2 4]);
% figure, plotregression(stocks.yr1Ret, stocks.lnSize), axis auto;
% axis([-2 4 18 25]);
figure, plotregression(stocks.yr1Ret, stocks.beMe), axis tight;
figure, plotregression(stocks.yr1Ret, stocks.m12M), axis tight;

% 4.4.2 correlation coefficient
[r1, p1] = corrcoef([stocks.lnSize stocks.beMe stocks.m12M]); 
% r1 = 
%     1.0000   -0.1408    0.0417
%    -0.1408    1.0000   -0.1597
%     0.0417   -0.1597    1.0000
% p1 = 
%     1.0000    0.0000    0.0224
%     0.0000    1.0000    0.0000
%     0.0224    0.0000    1.0000
[i, j] = find(p1 > 0.05); % [i j] = [1 1;2 2;3 3], the p-value matrix shows that the correlation coefficients between explanatory variables under 5% level are significant
[i, j] = find(p1 > 0.01); % while under 1% level, [i j] = [1 1;3 1;2 2;1 3;3 3], which shows that the correlation coefficients between lnSize and m12M is not significant
[r2, p2] = corrcoef([stocks.yr1Ret stocks.lnSize stocks.beMe stocks.m12M]);
r2 = r2(1, 2:end);
p2 = p2(1, 2:end);
% r2 =  0.0733    0.0424    0.1114, p2 = 0.0001    0.0201    0.0000
[i] = find(p2 > 0.05); % empty, which shows that all the correlation coefficients between explained variable and explanatory variables are significant under 5% level
[i] = find(p2 > 0.01); % i = 2, which shows that the correlation coefficients between yr1Ret and beMe is not significant under 1% level

% 4.5 Standardize and regression
stocks.zLnSize = zscore(stocks.lnSize);
stocks.zBeMe = zscore(stocks.beMe);
stocks.zM12m = zscore(stocks.m12M);
stocks.zYr1Ret = zscore(stocks.yr1Ret);
% stocks: 3000¡Á14 database
stats = regstats(stocks.zYr1Ret, [stocks.zLnSize stocks.zBeMe stocks.zM12m],...
    'linear', {'beta','tstat', 'fstat','mse','rsquare','adjrsquare','r'});
regression_table(stats);

% 4.6
%4.6.1 Fitness of the model - Adjusted R square
adj_rsquare = stats.adjrsquare;
% adj_rsquare = 0.0212 
% the adjusted rsquare is low, which means that this model does not fits
% the data well

% 4.6.2 significance of this model
q = 3; % there¡¯re 3 restrictions in the null hypothesis H0
           % H0: ŠÂ1=0, ŠÂ2=0, ŠÂ3=0 
n = length(stocks.zYr1Ret);
k = 3; % there's 3 regressors
fCrit_1 = finv(0.95, q, n-k-1); 
% fCrit_1 = 2.6079
%F_stat_1 = 22.6527, F_stat_1 > fCrit_1, so I can reject 
% the null hypothesis at level of 5%, which means
% that the regression is significant. 

% 4.6.3
% a. unrestricted sum of the squared residuals(USSR)
ussr = nansum(stats.r.^2);

% b. restricted regression
rStats = regstats(stocks.zYr1Ret, stocks.zM12m,...
    'linear', {'beta','tstat', 'fstat','mse','rsquare','adjrsquare','r'});

% c. restricted sum of the squared residuals (RSSR)
rssr = nansum(rStats.r.^2);

% d. F-test, H0: ŠÂ1=0, ŠÂ2=0
q = 2;
n = length(stocks.zYr1Ret);
k = 3;
F_stat_2 = ((rssr-ussr)/q)/(ussr/(n-k-1));

% e. critical value and P-value
fCrit_2 = finv(0.95, q, n-k-1); 
p_value = 1-fcdf(F_stat_2, q, n-k-1);
% p_value = 3.4251e-07, close to 0

% f. whether zLnSize and zBeMe are significant or not
% F_stat_2 = 14.9612, fCrit_2 = 2.9987
% F_stat_2 > fCrit_2, so I can reject the null 
% hypothesis at level of 5%, which means
% that zLnSize and zBeMe are significant in this model 
% and we should not drop these two regressors. 

% 4.7 
%              Coefficient  Std. Err.        t-stat          p-value    [95% Conf. Interval]
% ------------------------------------------------------------------------------
%       x1   0.078547   0.018251   4.303667   0.000017   0.042775   0.114319 
%       x2   0.072606   0.018472   3.930533   0.000087   0.036401   0.108812 
%       x3   0.119737   0.018304   6.541495   0.000000   0.083861   0.155613 

% From the table above, we can get that the slop coefficient of these 3
% factors are: 0.078547, 0.072606, 0.119737. The coefficients are not very
% big. Whild we can also see that p-values are all less than 5%, which
% means that though the coefficients are not big, these 3 factors do have a
% significant relationship with the explained varialbe. And the standard
% errors are also not big, and we can get the changing range of these 3
% factors in the 95% confidence interval. 

% Further more, use t-test to see whether these coefficients are significant for
% each individual factor
p = [0.025 0.975];
tCrit = tinv(p,n-1);
% tCrit = [-1.9608 1.9608]
% t_stat_zLnSize = 4.303667, t_stat_zBeMe = 3.930533, t_stat_zM12m = 6.541495, 
% these 3 t-stat are all greater than 1.96, so they are all significant
% under 5% level. So we can reject the null hypothesis and we can conclude
% that the slope coefficients are different from zero. 

% 4.8
% zYr1Ret = 0.0785zLnSize + 0.0726zBeMe + 0.1197zM12m
% from the resulting model we can see that size and BE/ME has a slight
% positive relationship with next year's return. From a view of all the
% companies listed in these 3 exchanges, an increase in the size can lead
% to a slight increase in next year's return. But the research of 3 factor
% model shows that frims with a small size have a higher return than firms
% with a largee size. In my point of view, I thought that the coefficient
% of the size should be nengative, but it is positive. Maybe the model is
% too simple or there's something wrong in selecting data. And the positive
% coefficient of Be/Me shows matches the research of 3 factor model at least 
% in direction, although the coefficient is not large. The coefficient of 12
% month momentum is also positive and is larger than the coefficient of size 
% and Be/Me, which means that firms with a higher return in 2004 will
% have a higher return in 2005. This matches the study of momentum which
% says that firms with a positive momentun in last 12 months and firms with
% a negative momentum last month tend to have a higher return in the
% future. 
