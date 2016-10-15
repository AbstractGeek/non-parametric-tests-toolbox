function [stats] = nemenyi(rankedData,tiedNums,alpha)
% Nemenyi test - nonparametric post-hoc test for multiple comparisions
% [stats] = nemenyi(rankedData,tiedNums,alpha)
%
% This test is supposed to be used after Kruskal Wallis test rejects the
% null hypothesis with a p-value of alpha
% 
% Inputs - 1)Cell array with ranks
%        - 2)tiedNums (number of ties per rank) 
%            (Obtained from the rankData function)
%        - 3)alpha. Valid alpha values are:
%            [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, .001]
% 
% Outputs - A table of multiple comparisons between different groups. The
% result of each comparisons is neatly provided in as table (cell array)
% 
% The test is described in Section 11.6 (nonparametric multiple
% comparisions) of Biostatistical analysis by Jerrold H. Zar
% 
% Author: Dinesh Natesan
%

load 'Nemenyi Test Critical Values.mat';
numOfInputs = length(rankedData);
N = NaN(numOfInputs,1);
R = NaN(numOfInputs,1);
for i=1:numOfInputs
    N(i) = length(rankedData{i,1});
    R(i) = sum(rankedData{i,1});
end
% Total N
sumN = sum(N);
t = sum(tiedNums.^3 - tiedNums);
meanRanks = R./N;
stats.meanRanksOrdered = meanRanks;
% Arrange mean ranks and preserve group number
stats.meanRanks = NaN(length(meanRanks),2);
[stats.meanRanks(:,1),stats.meanRanks(:,2)] = sort(meanRanks,'ascend');
meanRanks = stats.meanRanks;
% Get Q value and get the number of comparisions
Qexp=NemenyiCriticalValues(numOfInputs,NemenyiCriticalValues(1,:)==alpha);          %#ok<NODEF>
numOfComp = (numOfInputs*(numOfInputs-1)/2);         % Number of comparisions

% Define stats - output variable, and make it a proper table
stats.table = cell(numOfComp+1,7);
stats.table{1,1} = 'Comparision (B vs. A)';
stats.table{1,2} = 'Difference (mean(RA) - mean(RB))';
stats.table{1,3} = 'SE';
stats.table{1,4} = 'Q';
stats.table{1,5} = 'Q (obtained from the table))';
stats.table{1,6} = 'Null Hypothesis - H0';
stats.table{1,7} = 'Conclusion';

% Loop variables
currentB = numOfInputs;     % Ranks
currentA = 1;

flag = 0;
doNotCompare = {};

for i=1:numOfComp
    % Get necessary parameters for B
    RB = meanRanks(currentB,1);
    Bindex = meanRanks(currentB,2);
    nB = N(Bindex);
    % Get necessary parameters for A
    RA = meanRanks(currentA,1);
    Aindex = meanRanks(currentA,2);
    nA = N(Aindex);
    
    for j=1:length(doNotCompare)
        tmp = doNotCompare{j};
        if (sum(tmp==currentB)&&sum(tmp==currentA))
            flag=1;
        end
    end
    
    if (flag==1)
        stats.table{i+1,1} = strcat(int2str(Bindex),'vs. ',int2str(Aindex));
        stats.table{i+1,7} = 'Do Not Test';
        flag = 0;
    else
        % Get SE,Q and h
        SE = sqrt(((sumN*(sumN+1)/12) - (t/(12*(sumN-1))))*(1/nA+1/nB));
        Q = (RB - RA)/SE;
        h = (Q>Qexp);        
        
        if h
            conclusion = 'Reject H0';
        else
            conclusion = 'Accept H0';
            tmp = currentA:currentB;
            if(length(tmp)>2 && ~isempty(doNotCompare))
                dNClength = length( doNotCompare );
                doNotCompare{ dNClength + 1 } = tmp;
%                 doNotCompare = {doNotCompare{1:end};tmp};
            elseif (length(tmp)>2 && isempty(doNotCompare))
                doNotCompare = {tmp};
            end
        end
        
        stats.table{i+1,1} = strcat(int2str(Bindex),'vs. ',int2str(Aindex));
        stats.table{i+1,2} = RB-RA;
        stats.table{i+1,3} = SE;
        stats.table{i+1,4} = Q;
        stats.table{i+1,5} = Qexp;
        stats.table{i+1,6} = h;
        stats.table{i+1,7} = conclusion;
        
    end
    
    % Propogate the loop
    currentA = currentA+1;
    if (currentA == currentB && currentB==1)
        % loop done
        disp('Something wrong, your loop is not ending properly! Check!!');
    elseif (currentA == currentB)
        currentB = currentB-1;
        currentA = 1;
    end
    
end

end
