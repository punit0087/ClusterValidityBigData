function DI=DunnIndex(data,labels,distance,NUMmethod,DENmethod,DataSize)
%%%Dunn's index for clustering compactness and separation measurement
% data = data matrix [n X d]
% labels = [n X 1] vector of data labels
% distance  =   string that is passed to pdist to calculate distance between
%               data (e.g. 'euclidean', 'correlation', ...). If [],
%               'euclidean' is used as default.
% NUMmethod= selection of GDI for numerator
% DENmethod= selection of GDI for numerator
% -------------------------------------------------------------------------
%%Written by Punit Rathore, PhD Student, The University of Melbourne, Australia
% -------------------------------------------------------------------------


if ~exist('distance','var') || isempty(distance)
    distance = 'euclidean';
end

i= length(unique(labels));

switch DataSize
    case 'SmallData'
        numerator=[];
        NumInd=[];
        distM = squareform(pdist(data,distance));
        switch NUMmethod  %%for DI numerator
            
            case 'SL' %% single linkage
                for i2=1:i
                    i2
                    indi=find(labels==i2);
                    indj=find(labels~=i2);
                    [temp,temp1]=min(distM(indi,indj),[],2);
                    numerator=[numerator;min(temp)];
                end
            case 'CL' %% complete linkage
                for i2=1:i
                    i2
                    indi=find(labels==i2);
                    indj=find(labels~=i2);
                    temp=max(distM(indi,indj),[],2);
                    numerator=[numerator;max(temp)];
                end
            case 'AL'
                for i2=1:i %%average linkage
                    i2
                    indi=find(labels==i2);
                    indj=find(labels~=i2);
                                
%                     if length(labels)~=length(FullDataLables)
%                     indiOrig=length(find(FullDataLables==i2));
%                     indjOrig=length(find(FullDataLables~=i2));
%                       aC=  indjOrig/indiOrig;  
%                       bC=length(indj)/length(indi);
%                       NormFact=bC/aC;
%                     else
%                        NormFact=1;
%                     end
                    
                    temp= sum(distM(indi,indj),2)./length(indj);
                    numerator=[numerator; (sum(temp)/length(indi))];
                end
        end
        num= min(numerator);
        
        switch DENmethod %% for DI denominator
            case 'MaxD' %maximum diameter
                neg_obs=zeros(size(distM,1),size(distM,2));
                
                for ix=1:i
                    ix
                    indxs=find(labels==ix);
                    neg_obs(indxs,indxs)=1;
                end
                
                dem=neg_obs.*distM;
                denominator=max(dem);
            case 'AveD' %%average diameter
                denominator=[];
                for ix=1:i
                    ix
                    indxs=find(labels==ix);
                    MeanC=mean(data(labels==ix,:)); %%Changed on 2/10/2017 from Data to data and FullDataLables to labels
                    temp=mean(pdist2(data(indxs,:),MeanC,'euclidean'));
                    denominator=[denominator;2*temp];
                end
        end
        
        dem=max(denominator);
        
    case 'LargeData'
        
        switch NUMmethod
            case 'SL'
                numerator=Inf;
                for i2=1:i
                    i2
                    indi=find(labels==i2);
                    indj=find(labels~=i2);
                    for j=1:length(indi)
                        temp=min(pdist2(data(indi(j),:),data(indj,:),'euclidean'));
                        if temp<numerator
                            numerator=temp;
                        end
                    end
                end
                
            case 'CL'
                numerator=[];
                for i2=1:i
                    i2
                    Old_Value=-Inf;
                    indi=find(labels==i2);
                    indj=find(labels~=i2);
                    for j=1:length(indi)
                        temp=max(pdist2(data(indi(j),:),data(indj,:),'euclidean'));
                        if temp>Old_Value
                            Old_Value=temp;
                        end
                    end
                    numerator=[numerator;Old_Value];
                end
                
            case 'AL'
                temp=0;
                numerator=[];
                for i2=1:i
                    i2
                    indi=find(labels==i2);
                    indj=find(labels~=i2);
                    
                    for j=1:length(indi)
                        temp1=mean(pdist2(data(indi(j),:),data(indj,:),'euclidean'));
                        temp= ((j-1)*temp+temp1)/j;
                    end
                    numerator=[numerator;temp];
                end
        end
        num= min(numerator);
        
        switch DENmethod
            case 'MaxD'
                denominator=-Inf;
                for ix=1:i
                    ix
                    indxs=find(labels==ix);
                    for j=1:length(indxs)
                        temp=max(pdist2(data(indxs(j),:),data(indxs,:),'euclidean'));
                        if temp>denominator
                            denominator=temp;
                        end
                    end
                end
                
            case 'AveD'
                denominator=[];
                for ix=1:i
                    ix
                    indxs=find(labels==ix);
                    
                    temp=0;
                    for j=1:length(indxs)
                        temp1=pdist2(data(indxs(j),:),mean(data(indxs,:)),'euclidean');
                        temp= ((j-1)*temp+temp1)/j;
                    end
                    
                    denominator=[denominator;2*temp];
                end
        end
        dem=max(denominator);
        
end
DI=num/dem;
end
