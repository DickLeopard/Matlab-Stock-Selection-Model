%I can import the Excel data by using "Import Data" on the Toolbar, or I
%can import by using xlsread
clear all; close all; clc

data = xlsread('lab2_data.xlsx');

%2.2
data(any(isnan(data),2),:) = [];

%2.3
shsOut = data(: , 1);
clsgPrice = data(: , 2);
pE = data(: , 3);
Insize = log(shsOut.*clsgPrice);
eP = 1./pE;

%Plot histogram
nbins = 30;
figure, hist(Insize, nbins);
title('hist(Insize)');
figure, hist(eP, nbins);
title('hist(eP)');

%2.4
%2.4.1 Calculate mean and median
InsizeMean = mean(Insize);
InsizeMean;  %InsizeMean = 17.8259
InsizeMedian = median(Insize);
InsizeMedian; %InsizeMedian = 17.7952

ePMean = mean(eP);
ePMean; %ePMean = -0.1843
ePMedian = median(eP);
ePMedian; %ePMedian = 0.0052

%2.4.2 Calculate variance, standard deviation, and interquartile range
InsizeVar = var(Insize);
InsizeVar;
InsizeStd = std(Insize);
InsizeStd;
InsizeIqr = iqr(Insize);
InsizeIqr;
ePVar = var(eP);
ePVar;
ePStd = std(eP);
ePStd;
ePIqr = iqr(eP);
ePIqr;

%2.4.3 Calculate the percentiles
InsizePercentile = prctile(Insize, [5 25 50 75 95]);
InsizePercentile;
ePPercentile = prctile(eP, [5 25 50 75 95]);
ePPercentile;

%2.5
%2.5.1 Box plot
figure, boxplot(Insize, 'label', {'boxplot(Insize)'}, 'notch', 'on', 'whisker', 2);
figure, boxplot(eP, 'label', {'boxplot(eP)'}, 'notch', 'on', 'whisker', 2);

%2.5.2 Calculate and remove outliers
outliers_Insize = abs(zscore(Insize))>3;
Num_Outl_Insize = sum(outliers_Insize); %number of outliers
Insize(outliers_Insize) = []; %remove outtliers
outliers_eP = abs(zscore(eP))>3;
Num_Outl_eP = sum(outliers_eP); %number of outliers
eP(outliers_eP) = []; %remove outtliers

%2.6 Calculate and plot z-score
zscore_Insize = zscore(Insize);
figure, plot(zscore_Insize);
zscore_eP = zscore(eP);
figure, plot(zscore_eP);
%We standardize factor exposures in order to make the data
%comparable in magnitude and dispersion before further analysis.
%z-score tells us how many standard deviation units away from the
%average a particular stock’s exposure is, which facilitates stock
%comparison across multiple characteristics.
