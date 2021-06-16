function [neta,thresh]=otusthresh(a)
    gry=1:256;
b=imhist(a,256);
prob=b/length(a(:));
meantot=gry*prob;
for k=1:256
w0=sum(prob(1:k));
w1=1-w0;
 meancls0=(gry(1:k)*prob(1:k))/w0;
meancls1=(gry(k+1:256)*prob(k+1:256))/w1;
totvar=((gry-(meantot*ones(1,256))).^2)*prob;
varcls0=((gry(1:k)-(meancls0*ones(1,k))).^2)*prob(1:k);
varcls0(k)=varcls0/sum(prob(1:k));
varcls1=((gry(k+1:256)-(meancls1*ones(1,256-k))).^2)*prob(k+1:256);
varcls1(k)=varcls1/sum(prob(k+1:256));
varbtw=sum(prob(1:k))*sum(prob(k+1:256))*((meancls1-meancls0)^2);
neta(k)=varbtw/totvar;
end
l=max(neta);
fresh=find(neta==l);
thresh=mean(fresh);