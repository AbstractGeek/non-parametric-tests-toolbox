function [h,P] = kWallis(rankedData,tiedNums,alpha)
% Kruskal Wallis test (non parametric one way analysis of variance)
% [h,P] = kWallis(rankedData,tiedNums)
% 
% Inputs - 1)Cell array with ranks
%        - 2)tiedNums (number of ties per rank) 
% (Obtained from the rankData function)
% 
% Outputs - 1) h=0 implies failure to reject null hypothesis
%             h=1 implies that the null hypothesis is rejected
%           2) p-value for the null hypothesis
% 
% The test is described in Section 11.6 (nonparametric multiple
% comparisions) of Biostatistical analysis by Jerrold H. Zar
% 
% Author: Dinesh Natesan 
% 

numOfInputs = length(rankedData);
N = NaN(numOfInputs,1);
R = NaN(numOfInputs,1);
tempsum = 0;
for i=1:numOfInputs
    N(i) = length(rankedData{i,1});
    R(i) = sum(rankedData{i,1});
    tempsum = tempsum + (R(i)*R(i)/N(i));
end
% Total N
sumN = sum(N);
% Kruskal-Wallis Statistic, H
H = (12/(sumN*(sumN+1)))*tempsum - 3*(sumN+1);
% Correction Factor
t = sum(tiedNums.^3 - tiedNums);
C = 1-(t/(sumN.^3-sumN));
% Corrected H
Hcor = H/C;
v = numOfInputs-1;
P = 1-chi2cdf(Hcor,v);

if (P<alpha)
    h=1;            % Reject the Null hypothesis H0 (i.e. distributions are not the same)
else
    h=0;            % Failure to reject Null hypothesis H0 (i.e. distributions are same)
end

end