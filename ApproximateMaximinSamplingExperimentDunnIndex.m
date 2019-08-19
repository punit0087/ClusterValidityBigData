clc
close all
clear all


% load('data_matrix_with_lables');
% data_matrix_without_lables=data_matrix_with_lables(:,1:end-1);
%
% Data=data_matrix_without_lables;
% Labels=data_matrix_with_lables(:,end);
%
dimension=100;

[data_matrix] = CS_data_generate_Punit(-12,2,250000,dimension);
data_matrix_with_lables_1=[data_matrix zeros(250000,1)+1];

[data_matrix] = CS_data_generate_Punit(-3,1,250000,dimension);
data_matrix=[data_matrix zeros(250000,1)+2];
data_matrix_with_lables_1=[data_matrix_with_lables_1;data_matrix];

[data_matrix] = CS_data_generate_Punit(3,1,250000,dimension);
data_matrix=[data_matrix zeros(250000,1)+3];
data_matrix_with_lables_1=[data_matrix_with_lables_1;data_matrix];

[data_matrix] = CS_data_generate_Punit(12,2,250000,dimension);
data_matrix=[data_matrix zeros(250000,1)+4];
data_matrix_with_lables_1=[data_matrix_with_lables_1;data_matrix];

Data=data_matrix_with_lables_1(:,1:end-1);
Labels=data_matrix_with_lables_1(:,end);

clear data_matrix data_matrix_with_lables_1;

% load('GaussianMixtureData2WellSeperated.mat','Data','Labels');


% load('ForestCoverData.mat','Data','Labels')
% load('KDDCUPFullDATA.mat','Data','Labels')

n=size(Data,1);
INDX=randperm(n);

N=1;
Trials=11;

DN_New_CP_NS=zeros(Trials,N);
DN_New_CP=zeros(Trials,N);
DN_New_CP_NS_NMM=zeros(Trials,N);
NearlyMaximinTime=zeros(Trials,N);
MaximinTime=zeros(Trials,N);
DN_Old=zeros(1,N);
Old_DI_Time=zeros(1,N);
New_DI_Time=zeros(Trials,N);
New_DI_Time_NMM=zeros(Trials,N);


kk=1;
samples= n;
tt=INDX(1:samples);
Upspace_Mat=Data(tt,:);
lables= Labels(tt,:);
NoofK=length(unique(lables));
colors = distinguishable_colors(NoofK);


% tic
% DN_Old(kk)= DunnIndex(Upspace_Mat,lables,'euclidean','SL','MaxD','LargeData');
% Old_DI_Time(kk)=toc;


% figure;
% for i=1:NoofK
%     cluster_index=find(lables==i);
%     plot(Upspace_Mat(cluster_index,1),Upspace_Mat(cluster_index,2),'.','color','g');
%     hold on;
% end

ii=0;
for cp=[10,100,200,300,400,500,600,700,800,900,1000]
    ii=ii+1
    tic
    AproximateData_CP=[];
    AproximateData_CP_Ns=[];
    for i=1:NoofK
        i;
        indi=find(lables==i);
        currentcluster= Upspace_Mat(indi,:);
        clear indi
        tic
        [Currentsmp,~,m]= SamplingBigDatawithCP(currentcluster,cp,cp);
        MaximinTime(ii)=MaximinTime(ii)+toc
        
        AproximateData_CP=[AproximateData_CP;[currentcluster(m,:) ones(length(m),1)*i]];
        AproximateData_CP_Ns=[AproximateData_CP_Ns;[currentcluster(Currentsmp,:) ones(length(Currentsmp),1)*i]];
        
        
%         plot(currentcluster(Currentsmp1,1),currentcluster(Currentsmp1,2),'*','color','b','markers',15);hold on;
%         plot(currentcluster(Currentsmp2,1),currentcluster(Currentsmp2,2),'+','color','k','markers',15);hold on;
    end

    if size(AproximateData_CP_Ns,1)<10000
        %  DN_New(kk,ii)= indexDN(AproximateData(:,1:end-1),AproximateData(:,end),'euclidean');
        DN_New_CP_NS(ii)= DunnIndex(AproximateData_CP_Ns(:,1:end-1),AproximateData_CP_Ns(:,end),'euclidean','SL','MaxD','SmallData');
        New_DI_Time(ii)=toc;
        DN_New_CP(ii)= DunnIndex(AproximateData_CP(:,1:end-1),AproximateData_CP(:,end),'euclidean','SL','MaxD','SmallData');
    else
        %  DN_New(kk,ii)= indexDN1(AproximateData(:,1:end-1),AproximateData(:,end),'euclidean');
        DN_New_CP_NS(ii)= DunnIndex(AproximateData_CP_Ns(:,1:end-1),AproximateData_CP_Ns(:,end),'euclidean','SL','MaxD','LargeData');
        New_DI_Time(ii)=toc;
        DN_New_CP(ii)= DunnIndex(AproximateData_CP(:,1:end-1),AproximateData_CP(:,end),'euclidean','SL','MaxD','LargeData');
    end
    AppDataSize_MM_CP_NS=size(AproximateData_CP_Ns,1);
    AppDataSize_MM_CP=size(AproximateData_CP,1);
    
    clear AproximateData_CP AproximateData_CP_Ns
    
    
     %%%%%%%%%%%%%%%%%%%%%%%Approximate Maximin Sampling Method %% Didnt
     %%%%%%%%%%%%%%%%%%%%%%%work well
    AproximateData_CP_NearlyMM=[];
    for i=1:NoofK
        i;
        indi=find(lables==i);
        currentcluster= Upspace_Mat(indi,:);
        clear indi
        tic
        [Currentsmp1,~,m]= ApproximateMaximinSampling(currentcluster,cp,cp,30);
        NearlyMaximinTime(ii)=NearlyMaximinTime(ii)+toc
        
        AproximateData_CP_NearlyMM=[AproximateData_CP_NearlyMM;[currentcluster(Currentsmp1,:) ones(length(Currentsmp1),1)*i]];

%         plot(currentcluster(Currentsmp1,1),currentcluster(Currentsmp1,2),'*','color','b','markers',15);hold on;
%         plot(currentcluster(Currentsmp2,1),currentcluster(Currentsmp2,2),'+','color','k','markers',15);hold on;
    end

    if size(AproximateData_CP_NearlyMM,1)<10000
        DN_New_CP_NS_NMM(ii)= DunnIndex(AproximateData_CP_NearlyMM(:,1:end-1),AproximateData_CP_NearlyMM(:,end),'euclidean','SL','MaxD','SmallData');
        New_DI_Time_NMM(ii)=toc;
    else
        DN_New_CP_NS_NMM(ii)= DunnIndex(AproximateData_CP_NearlyMM(:,1:end-1),AproximateData_CP_NearlyMM(:,end),'euclidean','SL','MaxD','LargeData');
        New_DI_Time_NMM(ii)=toc;
    end
    AppDataSize_NMM=size(AproximateData_CP_NearlyMM,1);
     clear AproximateData_CP_NearlyMM 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
     plot(DN_Old(end)*ones(length(DN_New_CP(1:ii)),1),'-+','LineWidth',8,'color','g'); hold on;
%      plot(DN_New_CP(1:ii),'-*','LineWidth',8,'color','b'); hold on;
      plot(DN_New_CP_NS(1:ii),'-+','LineWidth',8,'color','r'); hold on;
      plot(DN_New_CP_NS_NMM(1:ii),'-*','LineWidth',8,'color','b'); hold on;

%       h_legend1=legend('OldDI','MaximinDI with C_p','MaximinDI with C_p and N_s');
       h_legend1=legend('OldDI','MaximinDI with C_p and N_s','Nearly Maximin Sampling');

      set(h_legend1,'FontSize',20,'Orientation','Horizontal');
     set(gca,'FontSize',15)
      set(gca,'Xtick',1:11,'xticklabel',{'10','100','200','300','400','500','600','700','800','900','1000'})
      xlabel('No of distingueshed points (c_{p})(n_{s} = c_{p})', 'FontSize', 30);
      ylabel({'Values'}, 'FontSize', 30)
      ylim([1  2])
%       xlim([0  100])
%     drawnow;
end













