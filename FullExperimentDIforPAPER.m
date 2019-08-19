clc
close all
clear all

data=importdata('hKernel.csv');
Data=data.data(:,[1,2]);
Labels=data.data(:,3)+1; %%Label with seperate anomalies

n=size(Data,1);
INDX=randperm(n);
samples=n;
tt=INDX(1:samples);
Upspace_Mat=Data(tt,:);
lables= Labels(tt,:);
NoofK=length(unique(lables))

figure(1);
plot(Upspace_Mat(:,1),Upspace_Mat(:,2),'.','color','r');
colors = distinguishable_colors(3)


tic %%Actual value of DI
%%following method runs iteratively to handle large size data
DN_Old= DunnIndex(Upspace_Mat,lables,'euclidean','SL','MaxD','LargeData')
%%following method is only good for small size data, N<5000-8000
%DN_Old= DunnIndex(Upspace_Mat,lables,'euclidean','SL','MaxD','SmallData')
Old_DI_Time=toc;

%%%%%%%%%%%%%%%aMMSR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%z
cp=100 %%choose some percent of data, 2-5%
ns=cp;
AproximateData_CP=[];
for i=1:NoofK
    i
    indi=find(lables==i);
    if length(indi)<cp
        cp= length(indi);
    end
    if cp<2
        cp=10;
    end
    currentcluster= Upspace_Mat(indi,:);
    [Currentsmp,~,m]= SamplingBigDatawithCP(currentcluster,cp,ns);
    AproximateData_CP=[AproximateData_CP;[currentcluster(m,:) ones(length(m),1)*i]];
end

if size(AproximateData_CP,1)<10000
    aMMSR_SL= DunnIndex(AproximateData_CP(:,1:end-1),AproximateData_CP(:,end),'euclidean','SL','MaxD','SmallData');
else
    aMMSR_SL= DunnIndex(AproximateData_CP(:,1:end-1),AproximateData_CP(:,end),'euclidean','SL','MaxD','LargeData');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%aNMMSR%%%%%%%%%%%%%%%%%%%%%
cp=100
ns=cp;
AproximateData_CP_Ns=[];
for i=1:NoofK
    i
    indi=find(lables==i);
    currentcluster= Upspace_Mat(indi,:);
    clear indi;
    [Currentsmp,~,m]= SamplingBigDatawithCP(currentcluster,cp,ns);
    Currentsmp=[Currentsmp;m];
    AproximateData_CP_Ns=[AproximateData_CP_Ns;[currentcluster(Currentsmp,:) ones(length(Currentsmp),1)*i]];
    %     plot(currentcluster(Currentsmp,1),currentcluster(Currentsmp,2),'o','color','r','markers',10);hold on;
end

if size(AproximateData_CP_Ns,1)<10000
    aNMMSR_SL= DunnIndex(AproximateData_CP_Ns(:,1:end-1),AproximateData_CP_Ns(:,end),'euclidean','SL','MaxD','SmallData');
else
    aNMMSR_SL= DunnIndex(AproximateData_CP_Ns(:,1:end-1),AproximateData_CP_Ns(:,end),'euclidean','SL','MaxD','LargeData');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=zeros(2,1);
figure;
for i=1:NoofK
    cluster_index=find(lables==i);
    plot(Upspace_Mat(cluster_index,1),Upspace_Mat(cluster_index,2),'.','color',colors(i,:));
    h(1)= plot(Upspace_Mat(cluster_index,1),Upspace_Mat(cluster_index,2),'.','color','g');
    hold on;
end
h(2)=plot(AproximateData_CP(:,1),AproximateData_CP(:,2),'o','color','k','markers',15);

legend(h,'Datapoints','i\alphaMMSR pts')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%inMMSR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk=1
R=cell(NoofK,1);
DN_New_CP_NS_iMM=0;
DN_New_CP_NS_iMM_Sorted=0;
DN=Inf;
DN1=Inf;
D=cell(NoofK,1);
AlphaPts=cell(NoofK,1);
X=cell(NoofK,1);
ApproxData=[];
for i=1:NoofK
    indi=find(lables==i);
    X{i}= Upspace_Mat(indi,:);
    initP= ceil(rand(1)*length(indi));
    D{i}=distance2(X{i}(initP,:),X{i})';
    R{i}(:,1)=D{i};
end
% clear Upspace_Mat
cp=0;
j=0;
figure;
hold on;
h_old=plot(0,0);
while cp<401 %% you can set max cp
    cp=cp+1
    ns=cp;
    j=j+1;
    smp=cell(NoofK,1);
    ApproxData2=[];
    for i=1:NoofK
        D{i}=min(D{i},R{i}(:,end));
        [max_d,m]=max(D{i});
        R{i}(:,end+1)=distance2(X{i}(m,:),X{i})';
        AlphaPts{i}(j)=m;
        ApproxData=[ApproxData; [X{i}(m,:) i]]; %%only for iMMRS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [dd,ii]=min(R{i},[],2);
        for t=1:cp
            s = find(ii==t+1);
            nt = ceil(ns*length(s)/(size(X{i},1)-length(find(ii==1))))-1; %% -1 written by Punit for Testing.
            ind = ceil(rand(nt,1)*length(s));
            smp{i}=[smp{i};s(ind)];
        end
        clear dd ii
        ApproxData2=[ApproxData2; [X{i}(smp{i},:) i*ones(length(smp{i}),1)]; [X{i}(AlphaPts{i},:) i*ones(length(AlphaPts{i}),1)]];
    end
    
    DN_New_CP_NS_iMM(j)=DunnIndex(ApproxData2(:,1:end-1),ApproxData2(:,end),'euclidean','SL','MaxD','SmallData');
    DN_New_CP_NS_iMM_Sorted=sort(DN_New_CP_NS_iMM,'descend');
    
    %  Use your criteria to terminate
    %    if (j>30) && (std(DN_New_CP_NS_iMM_Sorted(end-30:end))<=0.0001)
    %          j;
    %          break;
    %    end
    
    
    %%Ground Truth Value in DN_Old
    plot(DN_Old*ones(j,1),'-+','LineWidth',3,'color','g'); hold on;
    h= plot(DN_New_CP_NS_iMM_Sorted(1:end),'-+','LineWidth',8,'color','r'); hold on;
    
    xlabel('No. of Distinguished Points')
    ylabel('Dunn Index Value')
    legend('Original DI','iNMMSR');
    delete(h_old);
    h_old=h;
    xlim([0  500])
    drawnow;
end