function [x]=RisolviBidiagInf(d1,d2,y)
n=length(y);
x=zeros(n,1);
x(1)=y(1)/d1(1);
 
for k=2:n
    x(k)=(y(k)-d2(k-1)*x(k-1))/d1(k);
end
