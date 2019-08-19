function [smp,rp,m,temp_cell] = SamplingBigData(x, cp, ns )
[n,p]=size(x);

temp_cell=cell(1,cp);

[m,rp]=samplePlus(x,cp);

[~,i]=min(rp,[],2);
smp=[];
for t=1:cp,
    s = find(i==t); 
    nt = ceil(ns*length(s)/n); %% -1 written by Punit for Testing. 
    ind = ceil(rand(nt,1)*length(s));
    smp=[smp; s(ind)];
    temp_cell{t}=s(ind);
end;
smp=unique(smp); %% Written by Punit.
end