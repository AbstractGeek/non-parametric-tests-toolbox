function [h,P,stats] = PerformStats(X,Group,alpha)
% Perform Kruskal Wallis and Nemenyi test on the dataset with the provided
% alpha
% [h,p,stats] = PerformStats(X,Group,alpha)
% 
% Inputs: 1) X - observations [vector]
%         2) Group - group number for each observation [vector]
%         3) alpha. Valid alpha values are:
%            [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, .001]
%
% Outputs: 1) h from Kruskal Wallis test
%             h=0 implies failure to reject null hypothesis
%             h=1 implies that the null hypothesis is rejected
%          2) p-value from Kruskal Wallis test for the null hypothesis
%          3) Multiple comparison table from Nemenyi test (cell array)
% 
% Author: Dinesh Natesan 
% 


if nargin == 2
    alpha = 0.05;
end

% Split vector into subvectors and remove NaN's if there are any.
data = cell(max(Group),1);
for i = 1:size(X,1)
    if (isnan(X(i,1))~= 1)
        data{Group(i,1),1} = [data{Group(i,1),1};X(i,1)];
    end
end

% Rank the obtained data.
[rankedData,tiedNums] = rankData(data);
% Analysis of variance
[h,P] = kWallis(rankedData,tiedNums,alpha);
fprintf('\nHypothesis = %d; p = %f \n',h,P);
% Post-hoc analysis - Tukey test
[stats] = nemenyi(rankedData,tiedNums,alpha);
stats.X = X;
stats.Group = Group;

end
