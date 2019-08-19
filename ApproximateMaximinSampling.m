function [smp,rp,m] = ApproximateMaximinSampling(x, cp, ns,subsamplepercent)
[n,p]=size(x);

[m,rp]=ApproximatesamplePlus(x,cp,subsamplepercent);
% m=m(2:end);
% rp=rp(:,2:end);

[d,i]=min(rp,[],2);
smp=[];
if size(x,1)<cp
    smp=m;
else
    for t=1:length(m),
        t;
        s = find(i==t);
        nt = ceil(ns*length(s)/n) ; %% -1 written by Punit for Testing.
        [d_t,ind] = sort(d(s));
        smp=[smp; s(ind(1:nt))];
        %     ind = ceil(rand(nt,1)*length(s));
        %     smp=[smp; s(ind)];
    end;
end
smp=unique(smp); %% Written by Punit.
m=unique(m);
end