function [m,Rp,max_distance_Rp] = ApproximatesamplePlus(X,cp,subsamplepercent)
subsamplesize=ceil(subsamplepercent*size(X,1)/100);
[n,p]=size(X);
m=ones(cp,1);
max_distance_Rp=zeros(cp,1);
m(1)=ceil(rand(1)*size(X,1)); %%radonmly choose first point
Subsamples= ceil(rand(subsamplesize,1)*size(X,1)) %% select random p1% (subsample) of the original data
d=distance2(X(m(1),:),X(Subsamples,:))'; %% calculate distance of subsample from first point
Rp(:,1)=d;
for t=2:cp,
    d=min(d,Rp(:,t-1));
    [max_distance_Rp(t),m(t)]=max(d);
    Subsamples= ceil(rand(subsamplesize,1)*size(X,1)); %% select random p1% (subsample) of the original data
    Rp(:,t)=distance2(X(m(t),:),X(Subsamples,:))';
end;