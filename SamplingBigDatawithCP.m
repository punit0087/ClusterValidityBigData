function [smp,rp,m] = SamplingBigDatawithCP(x, cp, ns )
[n,p]=size(x);

[m,rp]=samplePlus(x,cp);
m=m(2:end);
rp=rp(:,2:end);

[d,i]=min(rp,[],2);
smp=[];
if size(x,1)<cp
    smp=m;
else
    for t=1:length(m),
        t;
        s = find(i==t);
        nt = ceil(ns*length(s)/n) ; %% -1 written by Punit for Testing.
        ind = ceil(rand(nt,1)*length(s));
        smp=[smp; s(ind)];
    end;
end
smp=unique(smp); %% Written by Punit.
end