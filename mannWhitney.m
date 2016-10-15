function [h,Uobs,Uexp] = mannWhitney(N1,N2)
% Mann-Whitney U Test (nonparametric t-test)
% [H,Uobs,Uexp] = mannWhitney(N1,N2)
% 
% Inputs - N1 (group 1)
%        - N2 (group 2)
% Ensure they are column vectors.
% 
% If significant (h=1), p value = 0.05
% 
% The test is described in Biostatistical analysis by Jerrold H. Zar.
% Alpha for this test is fixed at 0.05
% 
% Author: Dinesh Natesan
%

data{1,1} = N1;
data{2,1} = N2;
rankedData = rankData(data);
R1 = rankedData{1,1};
R2 = rankedData{2,1};
n1 = length(R1);
n2 = length(R2);
if (min(n1,n2)<3 && max(n1,n2)>5)
   disp('two small N. smallest of n1,n2 should be alteast 3 and largest of n1,n2 should be atleast 5');
   h = NaN;
   Uobs = NaN;
   Uexp = NaN;
   return;
end
U = (n1*n2) + (n1*(n1+1))/2 - sum(R1);
Udash = (n1*n2) - U;
Uobs = min(U,Udash);

load('Mann Whitney Critical Values.mat');
% row is smallest n, should be greater than 2
% column is largest n, should be greater than 4
% H0 - both the samples are from the same distribution

Uexp = mannWhitneyCriticalValues(min(n1,n2),max(n1,n2));
if (Uobs>Uexp)
    h = 0;          % Failure to reject Null hypothesis H0 (i.e. distributions are same)
else
    h = 1;          % Reject the Null hypothesis H0 (i.e. distributions are not the same)
    
end

end