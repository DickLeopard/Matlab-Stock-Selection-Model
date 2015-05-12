%3.1 Import Data

%%
exlServerObj = actxserver('Excel.Application');
exlWkbkObj = exlServerObj.Workbooks;
exlPath = 'D:\UMD\764 Quantitative Investment Strategy\Matlab Exercises';
exlFileList = 'lab3_data.xlsx';
exlFile = exlWkbkObj.Open([exlPath '\' exlFileList]);
exlSheet = exlFile.Sheets.Item('New Search Criteria');
exlRngObj = exlSheet.Range('A2:L3984');
exlData = exlRngObj.Value;
exlServerObj.Quit;
clear exlServerObj exlPath exlWkbkObj exlFileList exlFile exlSheet exlRngObj %only exlData left
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% delete data which has missing values
rowsIndx = find(prod(double(~cellfun('isempty',exlData)),2))';%1¡Á3983 row vector from 1 to 3983[1, 2, ¡­, 3983]
trimmedExlData = cell([numel(rowsIndx),size(exlData,2)]); %3983¡Á12 empty cells
count = 0;
for rowLoop = rowsIndx % rowLoop = 1 at the first time
    cellsWithNoTxt = ~cellfun('isclass',exlData(rowLoop,:),'char');% string:0, nonstring:1
    cellsWithNumsIndx = cellsWithNoTxt;%[0 0 0 1 1 1 1 1 1 1 1 1]
    invalidNums = isnan([exlData{rowLoop,cellsWithNumsIndx}]);% at the first time, for the 1st row of the number's part, NaN:1, nonNaN:0
    if ~any(invalidNums)% at the 1st time, for the 1st row, all numbers with no NaN: 1, has NaN: 0
        count = count + 1;
        for columnLoop = 1:size(exlData,2)% [1 2 ¡­ 12]
            trimmedExlData{count,columnLoop} = exlData{rowLoop,columnLoop};
        end
    end
end

trimmedExlData = trimmedExlData([1:count],:);
clear rowsIndx count rowLoop cellsWithNoTxt cellsWithNumsIndx invalidNums 
clear columnLoop exlData % only trimmedExlData left

%%
%Create the corresponding variables
name = {trimmedExlData{:, 1}}';
ticker = {trimmedExlData{:, 2}}';
exchange = {trimmedExlData{:, 3}}';
be = [trimmedExlData{:,4}]';
clsgPrice = [trimmedExlData{:, 5}]';
size = [trimmedExlData{:, 6}]';
cP = [trimmedExlData{:, 7}]';
eP = [trimmedExlData{:, 8}]';
ebitdaP = [trimmedExlData{:, 9}]';
fcfP = [trimmedExlData{:, 10}]';
sP = [trimmedExlData{:, 11}]';
yr1Ret = [trimmedExlData{:,12}]';

nyse = dataset({nominal(name),'name'},...   
    {nominal(ticker),'ticker'},{nominal(exchange),'exchange'},...    
    {be,'be'},{clsgPrice,'clsgPrice'},{size,'size'},...  
    {cP,'cP'},{eP,'eP'},{ebitdaP,'ebitdaP'},...   
    {fcfP,'fcfP'},{sP,'sP'},{yr1Ret,'yr1Ret'}); 
% dataset({size, 'SIZE'}): get a dataset that has a title of SIZE and follows the numbers in variable size. 
% name is string, if do not use nominal, there will be a '' around the names. For the number, do not need to use nominal function
varNames = get(nyse,'VarNames');

% 3.2 Find  and delete outliers
% remove outliers 
means = datasetfun(@mean,nyse,'DataVars',{'cP','eP','ebitdaP',... % calculate mean in nyse specified by Vars
    'fcfP','sP','yr1Ret'});                                                                            % means: 1¡Á6 row vector
stds = datasetfun(@std,nyse,'DataVars',{'cP','eP','ebitdaP',...
    'fcfP','sP','yr1Ret'});                                                                            %stds: 1¡Á6 row vector
outliersIndex = abs(nyse.cP - means(1)) > (3*stds(1))|...
    abs(nyse.eP - means(2)) > (3*stds(2))|...
    abs(nyse.ebitdaP - means(3)) > (3*stds(3))|...
    abs(nyse.fcfP - means(4)) > (3*stds(4))|...
    abs(nyse.sP - means(5)) > (3*stds(5))|...
    abs(nyse.yr1Ret - means(6)) > (3*stds(6));

nOut = sum(outliersIndex); % calculate the number of outliers
nyse(outliersIndex, :) = []; % delete the outliers

%% 3.3 Delete extreme values
ExtremeValue = nyse.be < 0|nyse.eP < 0|nyse.clsgPrice < 5|...
    nyse.size < 100000000;                                                                 % Find the extreme values
nyse(ExtremeValue, :) = [];                                                                 % Delete the extreme values

%% 3.4 Correlation Analysis
% 3.4.1 Scatter plot of 5 explanatory variables
figure, plotmatrix([nyse.cP nyse.eP nyse.ebitdaP nyse.fcfP nyse.sP]),...
    title({'scatter plotmatrix of explanatory variables' ; 'cP                    eP                 ebitdaP                  fcfP                   sP'});
ylabel('sP            fcfP            ebitdaP           eP              cP');


figure, subplot(5, 5, 1); hist(nyse.cP),ylabel('cP'), title('cP');                                                         % Use subplot function to simulate adding lsline to scatter plots got from plotmatrix
subplot(5, 5, 2); scatter(nyse.cP, nyse.eP, 1, 'filled'), lsline, title('eP');
subplot(5, 5, 3); scatter(nyse.cP, nyse.ebitdaP, 1, 'filled'), lsline, title('ebitdaP');
subplot(5, 5, 4); scatter(nyse.cP, nyse.fcfP, 1, 'filled'), lsline, title('fcfP');
subplot(5, 5, 5); scatter(nyse.cP, nyse.sP, 1, 'filled'), lsline, title('sP');
subplot(5, 5, 6); scatter(nyse.eP, nyse.cP, 1, 'filled'), lsline, ylabel('eP');
subplot(5, 5, 7); hist(nyse.eP);
subplot(5, 5, 8); scatter(nyse.eP, nyse.ebitdaP, 1, 'filled'), lsline;
subplot(5, 5, 9); scatter(nyse.eP, nyse.fcfP, 1, 'filled'), lsline;
subplot(5, 5, 10); scatter(nyse.eP, nyse.sP, 1, 'filled'), lsline;
subplot(5, 5, 11); scatter(nyse.ebitdaP, nyse.cP, 1, 'filled'), lsline, ylabel('ebitdaP');
subplot(5, 5, 12); scatter(nyse.ebitdaP, nyse.eP, 1, 'filled'), lsline;
subplot(5, 5, 13); hist(nyse.ebitdaP);
subplot(5, 5, 14); scatter(nyse.ebitdaP, nyse.fcfP, 1, 'filled'), lsline;
subplot(5, 5, 15); scatter(nyse.ebitdaP, nyse.sP, 1, 'filled'), lsline;
subplot(5, 5, 16); scatter(nyse.fcfP, nyse.cP, 1, 'filled'), lsline, ylabel('fcfP');
subplot(5, 5, 17); scatter(nyse.fcfP, nyse.eP, 1, 'filled'), lsline;
subplot(5, 5, 18); scatter(nyse.fcfP, nyse.ebitdaP, 1, 'filled'), lsline;
subplot(5, 5, 19); hist(nyse.fcfP);
subplot(5, 5, 20); scatter(nyse.fcfP, nyse.sP, 1, 'filled'), lsline;
subplot(5, 5, 21); scatter(nyse.sP, nyse.cP, 1, 'filled'), lsline, ylabel('sP');
subplot(5, 5, 22); scatter(nyse.sP, nyse.eP, 1, 'filled'), lsline;
subplot(5, 5, 23); scatter(nyse.sP, nyse.ebitdaP, 1, 'filled'), lsline;
subplot(5, 5, 24); scatter(nyse.sP, nyse.fcfP, 1, 'filled'), lsline;
subplot(5, 5, 25); hist(nyse.sP);

% figure, plotregression(nyse.cP,nyse.eP,'cP-eP',nyse.cP,nyse.ebitdaP,'cP-ebitdaP',nyse.cP,nyse.fcfP,'cP-fcfP',nyse.cP,nyse.sP,'cP-sP',...
%     nyse.eP, nyse.cP, 'eP-cP', nyse.eP, nyse.ebitdaP, 'eP-ebitdaP', nyse.eP, nyse.fcfP, 'eP-fcfP', nyse.eP,nyse.sP, 'eP-sP',...
%     nyse.ebitdaP, nyse.cP, 'ebitdaP-cP', nyse.ebitdaP, nyse.eP, 'ebitdaP-eP', nyse.ebitdaP, nyse.fcfP, 'ebitdaP-fcfP',...
%     nyse.ebitdaP,nyse.sP, 'ebitdaP-sP', nyse.fcfP, nyse.cP, 'fcfP-cP', nyse.fcfP, nyse.eP, 'fcfP-eP', nyse.fcfP, nyse.ebitdaP, 'fcfP-ebitdaP', ...
%     nyse.fcfP, nyse.sP, 'fcfP-sP', nyse.sP, nyse.cP, 'sP-cP', nyse.sP, nyse.eP, 'sP-eP', nyse.sP, nyse.ebitdaP, 'sP-ebitdaP',nyse.sP,nyse.fcfP, 'sP-fcfP' );
% 
% figure, plotregression(nyse.cP,nyse.eP,'cP-eP',nyse.cP,nyse.ebitdaP,'cP-ebitdaP',nyse.cP,nyse.fcfP,'cP-fcfP',nyse.cP,nyse.sP,'cP-sP',...
%     nyse.eP, nyse.ebitdaP, 'eP-ebitdaP', nyse.eP, nyse.fcfP, 'eP-fcfP', nyse.eP,nyse.sP, 'eP-sP',...
%     nyse.ebitdaP, nyse.fcfP, 'ebitdaP-fcfP',nyse.ebitdaP,nyse.sP,
%     'ebitdaP-sP',nyse.fcfP, nyse.sP, 'fcfP-sP');
% 
% Tried to use plotregression function, but can't read the figure. All the figures are crumpled up. 

%3.4.2 Correlation coefficients and 3.4.3 matrix of p-values
[r,p]=corrcoef([nyse.cP nyse.eP nyse.ebitdaP nyse.fcfP nyse.sP]); % r is the correlation coefficients matrix, and p is p-value matrix
[i,j]=find(p> 0.05); % [i, j] = [1 1; 2 2; 3 3; 4 4; 5 5], all on the diagnal. The p-value shows that the correlation coefficients are significant. 
[i,j]=find(p> 0.01); % result doesn't change, which shows that even under 1% level, the correlation coefficients are significant. 
% correlation coefficients matrix: 
% r = 
%     1.0000    0.2467    0.5546    0.4538    0.2193
%     0.2467    1.0000    0.4155    0.1930    0.2038
%     0.5546    0.4155    1.0000    0.0923    0.2941
%     0.4538    0.1930    0.0923    1.0000    0.1241
%     0.2193    0.2038    0.2941    0.1241    1.0000 
% From the correlation coefficient matrix r, I can get that corr(cP ebitdaP) = 0.5546, 
% corr(cP fcfP) = 0.4538, corr(eP ebitdaP) = 0.4155, these are 3 big corr, others are relatively small

% Matrix of p-values: 
% p = 
%     1.0000    0.0000    0.0000    0.0000    0.0000
%     0.0000    1.0000    0.0000    0.0000    0.0000
%     0.0000    0.0000    1.0000    0.0014    0.0000
%     0.0000    0.0000    0.0014    1.0000    0.0000
%     0.0000    0.0000    0.0000    0.0000    1.0000
% From the matrix of p-values I can get that all the p-values of regression
% between two different variables are smaller than 1%, which shows that the
% coefficient got above are all significant. 

%% 3.5 
% 3.5.1 scatter plot of explained variable and each of the explanatory variables
figure, subplot(3, 2, 1); hist(nyse.yr1Ret), ylabel('yr1Ret'); 
subplot(3, 2, 2); scatter(nyse.yr1Ret, nyse.cP, 1, 'filled'), lsline, ylabel('cP');
subplot(3, 2, 3); scatter(nyse.yr1Ret, nyse.eP, 1, 'filled'), lsline, ylabel('eP');
subplot(3, 2, 4); scatter(nyse.yr1Ret, nyse.ebitdaP, 1, 'filled'), lsline, ylabel('ebitdaP');
subplot(3, 2, 5); scatter(nyse.yr1Ret, nyse.fcfP, 1, 'filled'), lsline, ylabel('fcfP');
subplot(3, 2, 6); scatter(nyse.yr1Ret, nyse.sP, 1, 'filled'), lsline, ylabel('sP');

% 3.5.2  Correlation coefficients and 3.5.3 matrix of p-values
[r,p]=corrcoef([nyse.yr1Ret nyse.cP nyse.eP nyse.ebitdaP nyse.fcfP nyse.sP]);
r = r(1, :);
% correlation coefficients matrix: 
% r =  1.0000    0.0808    0.0121    0.0830   -0.0111    0.0500
% From the correlation coefficients matrix, I can get that the yearly
% returns are not highly correlated with the explanatory variables, and
% among the 5 explanatory variables, the yearly return has a negative
% relationship with fcfP, and has positive relationships with other 4 variables. 
p = p(1, :);
%  Matrix of p-values: 
% p = 1.0000    0.0052    0.6759    0.0040    0.7011    0.0836
i=find(p< 0.05); % i = 2, 4, which shows that correlation coefficient between yr1Ret and cP, 
                            % and between yr1Ret and ebitdaP are significant, others are not significant. 

%% 3.6 Standardize the variables using z-score
% 3.6.1 compute the z-scores of the explanatory variables and the correlation coefficients
nyse.zCP = zscore(nyse.cP);
nyse.zEP = zscore(nyse.eP);
nyse.zEbitdaP = zscore(nyse.ebitdaP);
nyse.zFcfP = zscore(nyse.fcfP);
nyse.zSP = zscore(nyse.sP);
[r,p]=corrcoef([nyse.zCP nyse.zEP nyse.zEbitdaP nyse.zFcfP nyse.zSP]);
% correlation coefficients matrix: 
% r =  
%     1.0000    0.2467    0.5546    0.4538    0.2193
%     0.2467    1.0000    0.4155    0.1930    0.2038
%     0.5546    0.4155    1.0000    0.0923    0.2941
%     0.4538    0.1930    0.0923    1.0000    0.1241
%     0.2193    0.2038    0.2941    0.1241    1.0000
%  Matrix of p-values: 
% p = 
%     1.0000    0.0000    0.0000    0.0000    0.0000
%     0.0000    1.0000    0.0000    0.0000    0.0000
%     0.0000    0.0000    1.0000    0.0014    0.0000
%     0.0000    0.0000    0.0014    1.0000    0.0000
%     0.0000    0.0000    0.0000    0.0000    1.0000
% The correlation coefficients matrix and matrix of p-values are the same
% with the results in 3.4. It can be easily proven that the correlation
% coefficients between z-scores of two variables are the same with the
% correlation coefficients between these two variables. 
% 3.6.2 Aggregate explanatory variables
nyse.zAggr = mean([nyse.zCP nyse.zEP nyse.zEbitdaP nyse.zFcfP nyse.zSP],2);
[r, p] = corrcoef([nyse.zAggr nyse.yr1Ret]);
r = r(2);
p = p(2);
% r = 0.0660, p = 0.0224, which shows that the explained variable yearly
% return has a low level of correlation between the aggregated explanatory
% variable, and the p-value of this test is 0.0224 < 5%, showing that it is
% significant under 5% level. 

%% 3.7 Problem of multicollinearity
% correlation coefficients matrix in 3.4: 
% r =                cP             eP        ebitdaP     fcfP           sP                          
%  cP            1.0000    0.2467    0.5546    0.4538    0.2193
%  eP            0.2467    1.0000    0.4155    0.1930    0.2038
%  ebitdaP   0.5546    0.4155    1.0000    0.0923    0.2941
%  fcfP         0.4538    0.1930    0.0923    1.0000    0.1241
%  sP            0.2193    0.2038    0.2941    0.1241    1.0000 
% Matrix of p-values: 
% p = 
%     1.0000    0.0000    0.0000    0.0000    0.0000
%     0.0000    1.0000    0.0000    0.0000    0.0000
%     0.0000    0.0000    1.0000    0.0014    0.0000
%     0.0000    0.0000    0.0014    1.0000    0.0000
%     0.0000    0.0000    0.0000    0.0000    1.0000
% From the matrix of p-values, I can get that all the correlation coefficients 
% are significant. And from correlation coefficients matrix we can see that 
% there are some pairs of variables that have big correlation coefficients:
% corr(cP ebitdaP) = 0.5546, corr(cP fcfP) = 0.4538, corr(eP ebitdaP) = 0.4155.
% These pairs of variables which has a big correlation coefficients can lead 
% to a problem of multicollinearity. 

% We can delete some variables to get rid of multicollinearity. For
% example, if a pair of variables have a big correlation coefficient, we
% can delete one of them. And there's other more sophisticated ways to deal
% with this problem, like stepwise regression. We can conduct regressions
% of explained variables with each of the explanatory variables and choose
% the one with the highest R-square at the 1st step. Then add other
% explanatory variables one by one into the regression model and keep the
% ones that will increase R-square. This may be a better way to deal with
% multicollinearity. 

%% 3.8 
% we would run a risk of overfitting a regression model if we incorporate
% all of the 5 explanatory variables to explain average return for the year
% 2005. Some of the variables are highly correlated, we don't need all of
% them. 
% step 1:
[r,p, rLo, rUp]=corrcoef([nyse.yr1Ret nyse.cP nyse.eP nyse.ebitdaP nyse.fcfP nyse.sP]);
r = r(1, 2:end);
p = p(1, 2:end);
rLo = rLo(1, 2:end);
rUp = rUp(1, 2:end);
% r = 0.0808    0.0121    0.0830   -0.0111    0.0500
% p = 0.0052    0.6759    0.0040    0.7011    0.0836
% the 2nd and 4th p-value are too big, so I will delete 2nd and 4th
% explanatory variables at the 1st step, that is delete eP and fcfP, and we
% have cP, ebitdaP and sP left. 

% step 2: 
[r, p] = corrcoef([nyse.cP nyse.ebitdaP nyse.sP]);
% r = 
%     1.0000    0.5546    0.2193
%     0.5546    1.0000    0.2941
%     0.2193    0.2941    1.0000
% p = 
%     1.0000    0.0000    0.0000
%     0.0000    1.0000    0.0000
%     0.0000    0.0000    1.0000
% we can see that corr(nyse.cP nyse.ebitdaP) = 0.5546, so I will delete one
% of them. In order to decide which one to delete, I conduct:
[r1,p1,rLo1,rUp1] = corrcoef([nyse.cP nyse.ebitdaP nyse.yr1Ret],'alpha',0.005);
% r1 = 
%     1.0000    0.5546    0.0808
%     0.5546    1.0000    0.0830
%     0.0808    0.0830    1.0000
% p1 = 
%     1.0000    0.0000    0.0052
%     0.0000    1.0000    0.0040
%     0.0052    0.0040    1.0000
% From r1 and p1 above, we can see that ebitdaP has a higher correlation
% coefficient than cP and the p-value of ebitdaP and less than the p-value
% of cP, so I will delete cP. 
% As a result, the model only has 2 variables left: ebitdaP and sP. 
