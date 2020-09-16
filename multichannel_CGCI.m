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


m = 10; % The order of the VAR model used for the computation of the 
        % GCI and CGCI 
alpha = 0.005; % The level of significance

%index=[5,6,7,8,9,13,14,15,21,22,23,24]; % select only the channel that are close to fc1
index=[5,6,7,8,9,13,14,21,22,23];
%index=[7,8,14,22,23];

K = length(index);  % The number of EEG channels (network nodes) to use (subset of the 
        % whole set of channels)
gcithresh = 0.05; % Arbitrary threshold for the significance of GCI and CGCI
maketest = 1; % 1-> make parametric significance test
samplefreq = 1450; % The sampling frequency
montagenr = 2; % -> 1, the old system, for EEG files E, F, 
               % -> 2, the new system for all others EEG file of 10-20 system 
taus = 1/samplefreq; % The sampling time
nameM=nameM(index);

rng(1);

bconnections=zeros(30,2);
wconnections=zeros(30,2);
figno=1;


for i=1:d3
    
    preTMS=xM(index,1451:2611,i);
    postTMS=xM(index,3191:4351,i);
    
    
    %% for window[-2000,-200] [200,2000]
    %{
    preTMS=xM(index,1:2610,i);
    postTMS=xM(index,3191:5800,i);
    %}
    %% For each pair of channels compute the conditional Granger causality 
    % index (CGCI) and form the CGCI-causality matrix
    fprintf('Computes the CGCI (m=%d) for all %d variables...\n',m,K);
    [CGCIM1,pCGCIM1] = CGCI(preTMS',m,maketest);
    [CGCIM2,pCGCIM2] = CGCI(postTMS',m,maketest);



    CGCIM1(1:K+1:end)=0;
    CGCIM1(CGCIM1<0)=0;


    CGCIM2(1:K+1:end)=0;
    CGCIM2(CGCIM2<0)=0;

    if size(pCGCIM1,1) == K
        pCGCIM1(1:K+1:end)=0;
    end

    if size(pCGCIM2,1) == K
        pCGCIM2(1:K+1:end)=0;
    end


    % The network of simple connections given by threshold GCI_{X->Y|Z}>thresh

    tit8txt = sprintf('GCI_{X->Y|Z}(%d) > %1.2f',m,gcithresh);
    %pre
    cgcithreshM1 = CGCIM1 > gcithresh;
    %plotnetworktitle(cgcithreshM1,[0 1],nameM,tit8txt,figno);
    %saveas(gca, fullfile('C:\Users\Armando\Documents\cgci\', num2str(figno)), 'jpeg');

    %post
    cgcithreshM2 = CGCIM2 > gcithresh;
    %plotnetworktitle(cgcithreshM2,[0 1],nameM,tit8txt,figno+1);
    %saveas(gca, fullfile('C:\Users\Armando\Documents\cgci\', num2str(figno+1)), 'jpeg');


    figno=figno+2;

    %pre
    fprintf("preTMS \n");
    fprintf('Average degree (for binary connections) or strength (for weighted connections): \n');
    fprintf('For the weighted connections:for CGCI=%3.2f \n',sum(sum(CGCIM1))/K);
    fprintf('For binary connections (threshohd=%1.5f):for CGCI=%3.2f \n',gcithresh,sum(sum(cgcithreshM1))/K)
    fprintf('\n');

    %post
    fprintf("postTMS \n");
    fprintf('Average degree (for binary connections) or strength (for weighted connections): \n');
    fprintf('For the weighted connections:for CGCI=%3.2f \n',sum(sum(CGCIM2))/K);
    fprintf('For binary connections (threshohd=%1.5f):for CGCI=%3.2f \n',gcithresh,sum(sum(cgcithreshM2))/K)
    fprintf('\n');

    bconnections(i,1)=sum(sum(cgcithreshM1))/K;
    bconnections(i,2)=sum(sum(cgcithreshM2))/K;
    
    wconnections(i,1)=sum(sum(abs(CGCIM1)))/K;
    wconnections(i,2)=sum(sum(abs(CGCIM2)))/K;

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

