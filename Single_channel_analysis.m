clear all

fprintf('Load the epochs with the rejected artifacts ...\n');
xM = load('test1.mat');
nameM=xM.EEG.chanlocs;
nameM=struct2cell(xM.EEG.chanlocs);
nameM=nameM(1,:)';
xM= xM.ans;

% select only the channel that are close to fc1
%{
index=[5,6,7,8,9,13,14,21,22,23];

xM=xM(index,:,:);
nameM=nameM(index);
%}

[d1,d2,d3]=size(xM);
%% calculate the mean for each channel
ch_mean=zeros(d1,d2);

for i=1:d1
    for j=1:d2
        ch_mean(i,j)=mean(xM(i,j,:));
    end
end

plotmts(ch_mean',1,0,d3,1/1450,nameM,1);

%plotmts(ch_mean(:,2900:3050)',1,0,d1,1/1450,nameM,1);

%% autocorrelation

tau=10;
H=zeros(d1,1);

for c=1:d1

    xm=xM(c,:,:);
    xm=squeeze(xm);
    
    x=zeros(d3,1);
    y=zeros(d3,1);

    for i=1:d3
        
        preTMS=xm(1451:2610,i);
        postTMS=xm(3191:4350,i);
            
        temp1=autocorrelation(preTMS,tau);
        x(i)=temp1(end,end);

        temp2=autocorrelation(postTMS,tau);
        y(i)=temp2(end,end);
    end
    %% performs a paired t-test of the hypothesis that two matched samples, in the vectors X and Y, come from distributions with equal means
    H(c)=ttest(x,y);
end

disp(H);


