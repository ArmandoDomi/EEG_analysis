function [pcM,ppcM] = mypartialcorr(xM,tau)
% [pcM,ppcM] = mypartialcorr(xM,tau)
% It computes the partial correlation for all time series in 'xM' (columns)
% and for the given lag 'tau'. In the output it gives also the p-value for
% each partial correlation from a parametric student test for significance.
% INPUTS
% - xM      : the n x K matrix of K time series.
% - tau     : the given lag (default=0 -> the K x K matrices of partial
%             correlations and p-values are symmetric)
% OUTPUTS
% - pcM     : the K x K matrix of partial correlations
% - ppcM    : the K x K matrix of p-values from the parametric student test 
%             for significance of the partial correlation for each pair of 
%             variables 

[n,K]=size(xM);
if nargin==1
    tau = 0;
end    
if tau<0
    tau = -tau;
end
%% Declare variables to store results in
pcM = ones(K,K);
ppcM = zeros(K,K);
%% Run for each pair of variables
for i=1:K-1
    if std(xM(:,i))==0
        break;
    end
    for j=i+1:K
        if std(xM(:,j))==0
            break;
        end
        %% Compute original r 
        % fprintf('pair (%d,%d) ... \n',i,j);
        xindV = setdiff([1:K],[i j]);
        if isempty(xindV)
            e1V = xM(:,i);
            e2V = xM(:,j);
        else
            [tmp1,tmp2,e1V] = regress(xM(:,i),[ones(n,1) xM(:,xindV)]);
            [tmp1,tmp2,e2V] = regress(xM(:,j),[ones(n,1) xM(:,xindV)]);
        end
        if tau==0
            tmpM = corrcoef(e1V,e2V);
            pcM(i,j)=tmpM(1,2);       
            %% Make a student parametric test of significance for phi 
            statnow = sqrt(n-2)*pcM(i,j)/sqrt(1-pcM(i,j).^2);
            ppcM(i,j) = 2*(1-tcdf(abs(statnow),n-2));
            pcM(j,i)=pcM(i,j);
            ppcM(j,i)=ppcM(i,j);
        else
            % if tau>0 then the partial correlation matrix is not symmetric        
            tmpM = corrcoef(e1V(1:n-tau),e2V(1+tau:n));
            pcM(i,j)=tmpM(1,2);       
            statnow = sqrt(n-tau-2)*pcM(i,j)/sqrt(1-pcM(i,j).^2);
            ppcM(i,j) = 2*(1-tcdf(abs(statnow),n-tau-2));
            tmpM = corrcoef(e1V(1+tau:n),e2V(1:n-tau));
            pcM(j,i)=tmpM(1,2);       
            statnow = sqrt(n-tau-2)*pcM(j,i)/sqrt(1-pcM(j,i).^2);
            ppcM(j,i) = 2*(1-tcdf(abs(statnow),n-tau-2));
        end    
    end % for j=i+1:K
end % for i=1:K-1
