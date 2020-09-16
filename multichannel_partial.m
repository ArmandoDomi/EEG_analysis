clear all

warning('off','all')
warning

fprintf('Load the epochs with the rejected artifacts ...\n');
xM = load('test1.mat');
nameM=xM.EEG.chanlocs;
nameM=struct2cell(xM.EEG.chanlocs);
nameM=nameM(1,:)';
xM=xM.ans;

[d1,d2,d3]=size(xM);

bconnections=zeros(30,2);
wconnections=zeros(30,2);
% Script 1: cross-correlation and partial correlation in EEG
% and networks 

tau = 10; % The delay for which the cross-correlation and partial 
         % correlation to be computed 
alpha = 0.005; % The level of significance
rthresh = 0.2; % Arbitrary threshold for the significance of 
               % cross-correlation and partial correlation
samplefreq = 1450; % The sampling frequency
montagenr = 2; % -> 1, the old system, for EEG files E, F, 
               % -> 2, the new system for all others EEG file of 10-20 system 
taus = 1/samplefreq; % The sampling time

% select only the channel that are close to fc1
index=[5,6,7,8,9,13,14,21,22,23];

K=length(index);

nameM=nameM(index);

figno=1;

rng(1);

for i=1:d3
    
    preTMS=xM(index,1451:2611,i);
    postTMS=xM(index,3191:4351,i);
    
    %% for window[-2000,-200] [200,2000]
    %{
    preTMS=xM(index,1:2610,i);
    postTMS=xM(index,3191:5800,i);
    %}
    %%
    
    
    fprintf('Computes the partial correlation (tau=%d) for all %d variables...\n',tau,K);
    %% For each pair of channels compute the partial correlation and form the 
    % partial correlation matrix
    [pcM1,ppcM1] = mypartialcorr(preTMS',tau);
    [pcM2,ppcM2] = mypartialcorr(postTMS',tau);

    pcM1(1:K+1:end)=0;
    ppcM1(1:K+1:end)=0;
    pcM2(1:K+1:end)=0;
    ppcM2(1:K+1:end)=0;
    
%% Plot the estimated partial correlation network
    % The network of weighted connections given by |r_{XY|Z}(tau)|
   
    % The network of simple connections given by threshold |r_{XY|Z}(tau)|>thresh
    pcthreshM1 = abs(pcM1) > rthresh;
    tit8txt = sprintf('|R_{XY|Z}(%d)| > %1.2f',tau,rthresh);
    %pre
    %plotnetworktitle(pcthreshM1,[0 1],nameM,tit8txt,figno);
    %saveas(gca, fullfile('C:\Users\Armando\Documents\partial\', num2str(figno)), 'jpeg');
    
    %post
    pcthreshM2 = abs(pcM2) > rthresh;
    %plotnetworktitle(pcthreshM2,[0 1],nameM,tit8txt,figno+1);
    %saveas(gca, fullfile('C:\Users\Armando\Documents\partial\', num2str(figno+1)), 'jpeg');
    
    figno=figno+2;

    %pre
    fprintf("preTMS \n");
    fprintf('Average degree (for binary connections) or strength (for weighted connections): \n');
    fprintf('For the weighted connections:for ParCorr=%3.2f \n',sum(sum(abs(pcM1)))/K);
    fprintf('For binary connections (threshohd=%1.5f): for ParCorr=%3.2f \n',rthresh,sum(sum(pcthreshM1))/K)
    fprintf('\n');
    
    %post
    fprintf("postTMS \n");
    fprintf('Average degree (for binary connections) or strength (for weighted connections): \n');
    fprintf('For the weighted connections:for ParCorr=%3.2f \n', sum(sum(abs(pcM2)))/K);
    fprintf('For binary connections (threshohd=%1.5f): for ParCorr=%3.2f \n',rthresh,sum(sum(pcthreshM2))/K)
    fprintf('\n');
    
    bconnections(i,1)=sum(sum(pcthreshM1))/K;
    bconnections(i,2)=sum(sum(pcthreshM2))/K;
    
    
    wconnections(i,1)=sum(sum(abs(pcM1)))/K;
    wconnections(i,2)=sum(sum(abs(pcM2)))/K;
end

%average degree
H=ttest(bconnections(:,1),bconnections(:,2));
%average strength
H_=ttest(wconnections(:,1),wconnections(:,2));

figure(figno)
plot(1:length(bconnections(:,1)),bconnections(:,1),1:length(wconnections(:,1)),wconnections(:,1));
title('PreTMS average strength / degree');
xlabel('epochs');
ylabel('strength / degree');
legend('degree', 'strength');

figure(figno+1)
plot(1:length(bconnections(:,2)),bconnections(:,2),1:length(wconnections(:,2)),wconnections(:,2));
title('PostTMS Average strength / degree');
xlabel('epochs');
ylabel('strength / degree');
legend('degree', 'strength');


