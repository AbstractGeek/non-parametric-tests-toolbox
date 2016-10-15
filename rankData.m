function [rankedData,tiedNums] = rankData(data)
% [rankedData,tiedNums] = rankData(data)
% This is a precursor function for the non parametric tests - Mann-Whitney and Kruskal Wallis test
% It outputs grouped data with its ranks
% 
% Input - Cell array with observations
% Output - 1)Cell array with ranks
%        - 2)tiedNums (number of ties per rank) 
% 
% Author: Dinesh Natesan 
% 


numOfInputs = length(data);
R = tiedrank(cell2mat(data));
rankedData = cell(numOfInputs,1);
% Break the tied ranks into groups and return the cell
temp = R;
for i = 1:numOfInputs
    cellend = size(data{i,1},1);
    rankedData{i} = temp(1:cellend);
    if i~=numOfInputs
        temp = temp(cellend+1:end);
    end
end

tiedNums = hist(R,(1:max(R))')';

end