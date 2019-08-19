clc
close all
clear all

% dimension=2;
% 
% [data_matrix] = CS_data_generate_Punit(-9,1,25000,dimension);
% data_matrix_with_lables_1=[data_matrix zeros(25000,1)+1];
% 
% [data_matrix] = CS_data_generate_Punit(0,2,50000,dimension);
% data_matrix=[data_matrix zeros(50000,1)+2];
% data_matrix_with_lables_1=[data_matrix_with_lables_1;data_matrix];
% 
% [data_matrix] = CS_data_generate_Punit(9,1,25000,dimension);
% data_matrix=[data_matrix zeros(25000,1)+3];
% data_matrix_with_lables_1=[data_matrix_with_lables_1;data_matrix];
% 
% Data=data_matrix_with_lables_1(:,1:end-1);
% Labels=data_matrix_with_lables_1(:,end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FOR MNIST Image Data
% images = loadMNISTImages('train-images.idx3-ubyte');
% labelsss = loadMNISTLabels('train-labels.idx1-ubyte');
% ImageData=[images' labelsss];
% Data=ImageData(:,1:end-1);
% Labels=ImageData(:,end)+1;
% 


% load('GaussianMixtureData2WellSeperated.mat','Data','Labels'); 
% load('ForestCoverData.mat','Data','Labels')
% load('KDDCUPFullDATA.mat','Data','Labels')


%% FOR HAR DATA

% Data_Train= importdata('X_train.txt');
% Label_Train= importdata('Y_train.txt');
% Data_Test= importdata('X_test.txt');
% Label_Test= importdata('Y_test.txt');
% Data=[Data_Train;Data_Test];
% Labels=[Label_Train;Label_Test];
% % cp=205 %%2 of total data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('ACTIVITYDATASETLARGEDATA.mat','ActivitySet');
% 
% fulldata=ActivitySet(:,1:end-1);
% Labels=ActivitySet(:,end);
% min10=min(fulldata,[],1);
% max10=max(fulldata,[],1);
% mean10=mean(fulldata,1);
% std10=std(fulldata,1);
% 
% minusminX = bsxfun(@minus,fulldata,min10);
% Data = bsxfun(@rdivide,minusminX,(max10-min10));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=size(Data,1);
INDX=randperm(n);

N=1;
Trials=1;

DN_New_CP_NS=zeros(Trials,N);
DN_New_CP=zeros(Trials,N);
% DN_Old=zeros(Trials,N);
% Old_DI_Time=zeros(Trials,N);
New_DI_Time=zeros(Trials,N);

kk=1;
samples= n%50000;
tt=INDX(1:samples);
Upspace_Mat=Data(tt,:);
lables= Labels(tt,:);
NoofK=length(unique(lables));
colors = distinguishable_colors(NoofK);




% figure;
% for i=1:NoofK
%     cluster_index=find(lables==i);
%     %      plot(Upspace_Mat(cluster_index,1),Upspace_Mat(cluster_index,2),'.','color',colors(i,:));
%     plot(Upspace_Mat(cluster_index,1),Upspace_Mat(cluster_index,2),'.','color','g');
%     hold on;
% end


% tic
% DN_Old(kk)= DunnIndex(Upspace_Mat,lables,'euclidean','SL','MaxD','LargeData')
% Old_DI_Time(kk)=toc;


    
% DN_Old=0.0060
R=cell(NoofK,1);
DN_New_CP_NS_iMM=0;
DN_New_CP_NS=0;
DN_New_CP=0;
DN=Inf;
DN1=Inf;
ff=[];
D=cell(NoofK,1);
AlphaPts=cell(NoofK,1);
X=cell(NoofK,1);
ApproxData=[];
for i=1:NoofK
    indi=find(lables==i);
    X{i}= Upspace_Mat(indi,:);
    D{i}=distance2(X{i}(1,:),X{i})';
    R{i}(:,1)=D{i};
end
cp=0;
j=0;
   
numerator=Inf;
denominator=-Inf;

while cp<1000
    cp=cp+1
    ns=2*cp;
    j=j+1;
    smp=cell(NoofK,1);
    smp1=cell(NoofK,1);
    ApproxData1=[];
    ApproxData2=[];
   
     for i=1:NoofK
     D{i}=min(D{i},R{i}(:,end));
     [max_d,m]=max(D{i});
     R{i}(:,end+1)=distance2(X{i}(m,:),X{i})';
     AlphaPts{i}(j)=m;  
%      ApproxData=[ApproxData; [X{i}(m,:) i]];
     [dd,ii]=min(R{i},[],2);
     for t=1:cp
     s = find(ii==t+1);
     nt = ceil(ns*length(s)/(size(X{i},1)-length(find(ii==1))))-1; %% -1 written by Punit for Testing.
     ind = ceil(rand(nt,1)*length(s));     
%      [d_t,ind] = sort(dd(s));
     smp1{i}=[smp1{i};s(ind)];
     
%      nt = ceil(ns*length(s)/(size(X{i},1)-length(find(ii==1))))-1 %% -1 written by Punit for Testing.
     nt = ceil(ns*length(s)/(size(X{i},1)))-1; %% -1 written by Punit for Testing.
     ff=[ff;nt];
     ind = ceil(rand(nt,1)*length(s));
     smp{i}=[smp{i}; s(ind)];
     end   
     [cp length(smp{i}) length(AlphaPts{i})];
     ApproxData1=[ApproxData1; [X{i}(smp{i},:) i*ones(length(smp{i}),1)]; [X{i}(AlphaPts{i},:) i*ones(length(AlphaPts{i}),1)]];
     ApproxData2=[ApproxData2; [X{i}(smp1{i},:) i*ones(length(smp1{i}),1)]; [X{i}(AlphaPts{i},:) i*ones(length(AlphaPts{i}),1)]];
%      ApproxData2=[ApproxData2; [X{i}(smp1{i},:) i*ones(length(smp1{i}),1)]];
     end   
%    plot(ApproxData1(:,1),ApproxData1(:,2),'*','color','k','markers',15);hold on;

  DN_New_CP_NS_iMM(j)= min(DN1,DunnIndex(ApproxData2(:,1:end-1),ApproxData2(:,end),'euclidean','SL','MaxD','SmallData'));
  DN_New_CP_NS(j)= min(DN, DunnIndex(ApproxData1(:,1:end-1),ApproxData1(:,end),'euclidean','SL','MaxD','SmallData'));
  DN=DN_New_CP_NS(j);
  DN1=DN_New_CP_NS_iMM(j);

%    plot(DN_Old*ones(j,1),'-+','LineWidth',8,'color','g'); hold on;
%    plot(DN_New_CP(1:j),'-+','LineWidth',8,'color','k'); hold on;
  plot(DN_New_CP_NS_iMM(1:j),'-+','LineWidth',8,'color','r'); hold on;
  plot(DN_New_CP_NS(1:j),'-+','LineWidth',8,'color','b'); hold on;
 ylim([-0.5  0.5])
  xlim([0  500])
  drawnow;
end









