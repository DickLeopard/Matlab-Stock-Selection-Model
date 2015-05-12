
% PART A | Initiation 
clear all; close all; clc

exlPath = 'D:\UMD\764 Quantitative Investment Strategy\Matlab Exercises'; 
exlFileListInSample =...
    {'lab5_data_2003_2004.xlsx', 'lab5_data_2004_2005.xlsx',...
    'lab5_data_2005_2006.xlsx', 'lab5_data_2006_2007.xlsx','lab5_data_2007_2008.xlsx'};
exlFileListOutOfSample =...
    {'lab5_data_2008_2009.xlsx', 'lab5_data_2009_2010.xlsx'}; 
exlSheet = 'New Search Criteria'; 
exlRange = 'A2:I8871'; 

%% PART B | In-Sample Estimation
% Beginning of loop over list of in-sample Excel files    
coefficients = zeros(5);
for fileCount = 1:numel(exlFileListInSample)
    exlFile = exlFileListInSample{fileCount};
    
    [stocks.name stocks.ticker stocks.exchange stocks.clsgPrice...
        stocks.size stocks.beMe stocks.eP stocks.m12M stocks.yr1Ret]...
        = exlDatafromMorningstar([exlPath '\' exlFile],...
        exlSheet,exlRange);
 % 5.2.3 filters the data
    filteredIndex = stocks.beMe > 0 & stocks.eP > 0 &...
        stocks.clsgPrice > 5 & stocks.size >= 100000000;
    stocks.name = {stocks.name{filteredIndex}}';
    stocks.ticker = {stocks.ticker{filteredIndex}}';
    stocks.exchange = {stocks.exchange{filteredIndex}}';
    stocks.clsgPrice = stocks.clsgPrice(filteredIndex);
    stocks.size = stocks.size(filteredIndex);
    stocks.beMe = stocks.beMe(filteredIndex);
    stocks.eP = stocks.eP(filteredIndex);
    stocks.m12M = stocks.m12M(filteredIndex);
    stocks.yr1Ret = stocks.yr1Ret(filteredIndex);
    
    % 5.2.2 remove outliers
    stocks = struct2dataset(stocks);
    means = datasetfun(@mean,stocks,'DataVars',{'clsgPrice', 'size','beMe','eP', 'm12M','yr1Ret'}); 
    stds = datasetfun(@std,stocks,'DataVars',{'clsgPrice', 'size','beMe','eP', 'm12M','yr1Ret'}); 
    outliersIndx = abs(stocks.clsgPrice - means(1)) > (3*stds(1)) |...
    abs(stocks.size - means(2)) > (3*stds(2)) |...
    abs(stocks.beMe - means(3)) > (3*stds(3)) |...
    abs(stocks.eP - means(4)) > (3*stds(4)) |...
    abs(stocks.m12M - means(5)) > (3*stds(5)) |...
    abs(stocks.yr1Ret - means(6)) > (3*stds(6));
    stocks(outliersIndx,:) = []; %  ____ stocks left
    
    % 5.2.4 log of the market return of equity and zscores
    stocks.lnSize = log(stocks.size);
    stocks.zLnSize = zscore(stocks.lnSize);
    stocks.zBeMe = zscore(stocks.beMe);
    stocks.zEP = zscore(stocks.eP);
    stocks.zM12M = zscore(stocks.m12M);
    stocks.zyr1Ret = zscore(stocks.yr1Ret);
    
    % 5.2.5&5.2.6 regression and the output table
    disp(['Regression for file ' exlFile])
    stats = regstats(stocks.zyr1Ret,...
        [stocks.zLnSize stocks.zBeMe stocks.zEP stocks.zM12M],'linear',...
        {'beta','tstat','fstat','mse','rsquare','adjrsquare','r'});
    regression_table(stats)
    disp(' ');
     
    % 5.2.7 save all the beta into a matrix
        coefficients(:,fileCount) = stats.beta;
    
    clear stocks
    % End of loop over list of in-sample Excel files    
end
clear fileCount exlFile noOutliersIndx filteredIndex...
    fileCount exlFileListInSample stats means stds outliersIndx

%% 5.3 Fama-Macbeth procedure
% 5.3.1 
coefficientsMean = mean(coefficients, 2);
% 5.3.2
coefficientsStd = std(coefficients, 1, 2);
% 5.3.3
t_stats = coefficientsMean./coefficientsStd*sqrt(5);
% 5.3.4
% a. the estimated factor return premia is the amount paid by the market
% for each standard deviation of exposure to each factor during the
% estimation window and we want to use them as estimated amount paid 
% by the market for each standard deviation of exposure during the
% forecasting window. And the FM t-statistics for each factor premia are
% used to see whether these estimation are statistically significant or not

% b. when the windows are non-overlapping and independent, we can derive
% the FM t-statistics easily. If there are overlappings, we need to do some
% adjustment. 
% clear coefficientsMean coefficientsStd
%% 5.4 select factors for forecasting windows
% t-statistics =    -0.5319    1.7989    1.4100    0.2035    0.0041; 
% the 2nd one > 1.5 and the 3rd one is close to 1.5. 
% the 2nd one is the premia of lnSize, according to 3 factor model, this
% one should be negative, but there are some exceptions. For example,
% professor mensioned on Monday that during some years in the early 2000's,
% the returns of big firms exceed the returns of small firms. So I will
% keep the 2nd one though the sigh should be negative according to 3 factor
% model. And the 3rd one is for Be/Me. According to 3 factor model, companies 
% with hign Be/Me ratios will have a higher return than companies with low
% Be/Me ratios. So the sign is correct. Though the t-statistic is less than
% 1.5, it is 1.41, close to 1.5. I will keep this one to test in the
% forecasting windows. 

%% PART C | Out-of-Sample Test
% 5.5.1 imports and cleans for missing and erroneous entries
portValue = zeros(2);
portRet = zeros(2);
for fileCount = 1:numel(exlFileListOutOfSample)
    exlFile = exlFileListOutOfSample{fileCount};
    [stocks.name stocks.ticker stocks.exchange stocks.clsgPrice...
    stocks.size stocks.beMe stocks.eP stocks.m12M stocks.yr1Ret] =...
    exlDatafromMorningstar([exlPath '\' exlFile],...
    exlSheet,exlRange);
% 5.5.2 % 5.5.3 filter the data and remove the outliers
    filteredIndex = stocks.beMe > 0 & stocks.eP > 0 &... 
        stocks.clsgPrice > 5 & stocks.size >= 100000000;
    stocks.name = {stocks.name{filteredIndex}}';
    stocks.ticker = {stocks.ticker{filteredIndex}}';
    stocks.exchange = {stocks.exchange{filteredIndex}}';
    stocks.clsgPrice = stocks.clsgPrice(filteredIndex);
    stocks.size = stocks.size(filteredIndex);
    stocks.beMe = stocks.beMe(filteredIndex);
    stocks.eP = stocks.eP(filteredIndex);
    stocks.m12M = stocks.m12M(filteredIndex);
    stocks.yr1Ret = stocks.yr1Ret(filteredIndex);
    
    stocks = struct2dataset(stocks); % 2395¡Á9; 2174¡Á9
    means = datasetfun(@mean,stocks,'DataVars',{'clsgPrice', 'size','beMe','eP', 'm12M','yr1Ret'}); 
    stds = datasetfun(@std,stocks,'DataVars',{'clsgPrice', 'size','beMe','eP', 'm12M','yr1Ret'}); 
    outliersIndx = abs(stocks.clsgPrice - means(1)) > (3*stds(1)) |...
    abs(stocks.size - means(2)) > (3*stds(2)) |...
    abs(stocks.beMe - means(3)) > (3*stds(3)) |...
    abs(stocks.eP - means(4)) > (3*stds(4)) |...
    abs(stocks.m12M - means(4)) > (3*stds(5)) |...
    abs(stocks.yr1Ret - means(4)) > (3*stds(6));
    stocks(outliersIndx,:) = []; %  2127¡Á9; 2054¡Á9
    
    % 5.5.4 create lnSize and compute zscores 
   stocks.lnSize = log(stocks.size);
   stocks.zLnSize = zscore(stocks.lnSize);
   stocks.zBeMe = zscore(stocks.beMe);
   stocks.zEP = zscore(stocks.eP);
   stocks.zM12M = zscore(stocks.m12M);
   % stocks: 2027¡Á14; 2054¡Á14
    
    scores = t_stats(2) * stocks.zLnSize + t_stats(3) * stocks.zBeMe;
%     scores = t_stats(3) * stocks.zBeMe + t_stats(4) * stocks.eP;
    % scores: 2127¡Á1; 2054¡Á1
    % each stock has a score in the list 'scores'
    % stocks: dataset; scores: vector

    dollarsInvested = 10000000;
    numberStocksLong = round(numel(scores)*0.10);
    numberStocksShort = numberStocksLong;
    % long top 10% and short bottom 10%

    [sortedScores, scoreSortIndx] = sort(scores);
    hiScoreIndx = [numel(scores):-1:(numel(scores) - numberStocksLong + 1)];
    lowScoreIndx = [numberStocksShort:-1:1];
    disp(' ')
    disp('Portfolio stocks, sorted by score')
   
    for kk = hiScoreIndx
        disp(['Long (order by score): '  stocks.name{scoreSortIndx(kk)}...
            ' (' stocks.ticker{scoreSortIndx(kk)}...
            '): Score = ' num2str(scores(scoreSortIndx(kk)))])
    end
    for kk = lowScoreIndx
      disp(['Short (order by score): '  stocks.name{scoreSortIndx(kk)}...
           ' (' stocks.ticker{scoreSortIndx(kk)}...
         '): Score = ' num2str(scores(scoreSortIndx(kk)))])
    end

    disp(' ')
    disp('Portfolio stocks, sorted by performance')
    longScoresIndx = scoreSortIndx(hiScoreIndx);
    longRets = 100 * stocks.yr1Ret(longScoresIndx);
    longNames = stocks.name{longScoresIndx};
    longTickers = stocks.ticker{longScoresIndx};
    longScores = scores(longScoresIndx);
    [sortedLongRets sortedLongIndx] = sort(longRets);
    for kk = numberStocksLong:-1:1
        disp(['Long (order by performance): '  ...
            stocks.name{longScoresIndx(sortedLongIndx(kk))}...
            ' (' stocks.ticker{longScoresIndx(sortedLongIndx(kk))}...
            '): Score = ' num2str(longScores(sortedLongIndx(kk)))...
           ', PercentReturn = ' num2str(longRets(sortedLongIndx(kk)))])
    end
    shortScoresIndx = scoreSortIndx(lowScoreIndx);
    shortRets = 100 * stocks.yr1Ret(shortScoresIndx);
    shortNames = stocks.name{shortScoresIndx};
    shortTickers = stocks.ticker{shortScoresIndx};
    shortScores = scores(shortScoresIndx);
    [sortedShortRets sortedShortIndx] = sort(shortRets);
    for kk = 1:numberStocksShort
      disp(['Short (order by performance): '  ...
            stocks.name{shortScoresIndx(sortedShortIndx(kk))}...
            ' (' stocks.ticker{shortScoresIndx(sortedShortIndx(kk))}...
           '): Score = ' num2str(shortScores(sortedShortIndx(kk)))...
            ', PercentReturn = ' num2str(-shortRets(sortedShortIndx(kk)))])
    end

    longPortValue = 0;
    for kk = 1:numberStocksLong
        longPortValue = longPortValue +...
            dollarsInvested / numberStocksLong *...
            (1 + stocks.yr1Ret(scoreSortIndx(end - kk + 1)));
    end
    shortPortValue = 0;
    for kk = 1:numberStocksShort
        shortPortValue = shortPortValue +...
           dollarsInvested / numberStocksShort *...
           (1 - stocks.yr1Ret(scoreSortIndx(kk)));
    end
    portValue(1, fileCount) = longPortValue;
    portValue(2, fileCount) = shortPortValue;
    longPortRetPer = (longPortValue/dollarsInvested-1)*100;
    shortPortRetPer = (shortPortValue/dollarsInvested-1)*100;
    portRet(1, fileCount) = longPortRetPer;
    portRet(2, fileCount) = shortPortRetPer;
    clear stocks scores
end
clear kk fileCount exlFile filteredIndex means stds outliersIndx dollarsInvested ...
    numberStocksLong numberStocksShort sortedScores scoreSortIndx hiScoreIndx...
    lowScoreIndx longScoresIndx longRets longNames longTickers longScores...
    sortedLongRets sortedLongIndx shortScoresIndx shortRets shortNames...
    shortTickers shortScores sortedShortRets sortedShortIndx longPortValue...
    shortPortValue longPortRetPer shortPortRetPer
 
% the selected stock list is too long, I put it in another file

% 2008-2009: 
% longPortValue = 1.4466e+07
% shortPortValue =7.5801e+06
% long-return = ((1.4466e+07/10000000)-1)*100 = 44.67%
% short-return = ((7.5801e+06/10000000)-1)*100 = -24.19%
%    
% 2009-2010: 
% longPortValue = 1.1659e+07
% shortPortValue = 8.7180e+06
% long-return = ((1.1659e+07/10000000)-1)*100 = 16.59%
% short-return = ((8.7180e+06/10000000)-1)*100 = -12.82%

% 5.6 8 estimation windows and 2 forecasting windows
exlFileListInSample =...
    {'lab5_data_2000_2001.xlsx', 'lab5_data_2001_2002.xlsx', 'lab5_data_2002_2003.xlsx', 'lab5_data_2003_2004.xlsx', 'lab5_data_2004_2005.xlsx',...
    'lab5_data_2005_2006.xlsx', 'lab5_data_2006_2007.xlsx','lab5_data_2007_2008.xlsx'};
coefficients = zeros(5, 8);
t_stats = coefficientsMean./coefficientsStd*sqrt(8);
% modify the codes above and run from 5.1-5.3
% t-statistics' =   -1.1012   -0.3148    2.3755    2.0447    0.0425
% the 3rd beta which is for BeMe and 4th beta which is for eP are
% significant. And the beta for BeMe is positive which matches the study of
% 3 factor model, and the positive beta for eP ratio also make sence for a
% higher earning will lead to a higher return the next year. So this time I
% will choose BeMe and eP. 

% only modify one line and run part 5.5: 

% scores = t_stats(3) * stocks.zBeMe + t_stats(4) * stocks.eP;

% 2008-2009: 
% longPortValue = 1.4438e+07
% shortPortValue =7.5155e+06
% long-return = ((1.4466e+07/10000000)-1)*100 = 44.38%
% short-return = ((7.5801e+06/10000000)-1)*100 = -24.85%
%    
% 2009-2010: 
% longPortValue = 1.2286e+07
% shortPortValue = 8.0174e+06
% long-return = ((1.1659e+07/10000000)-1)*100 = 22.86%
% short-return = ((8.7180e+06/10000000)-1)*100 = -19.83%

% 5.7 difference 
% in the 5 estimation windows situation, I chose lnSize and BeMe, 
% while in the 8 estimation windows situation, I chose BeMe and 
% eP. This change lead to the change in the list of scores. So that the 
% stocks selected according to the sorted scores were different. And as a
% result the long-short portfolio values and returns are different in two
% different situations. 












