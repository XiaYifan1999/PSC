function [X, Y] = normr(x,y)
% x, y are normalized separately 

x = double(x);
y = double(y);

n=size(x,1);
m=size(y,1);

xm=mean(x);
ym=mean(y);

x=x-repmat(xm,n,1);
y=y-repmat(ym,m,1);

xscale=1./max(abs(x));
yscale=1./max(abs(y));

X=x.*repmat(xscale,n,1)*0.4;
Y=y.*repmat(yscale,m,1)*0.4;

X=X+[0.5,0.5];
Y=Y+[0.5,0.5];