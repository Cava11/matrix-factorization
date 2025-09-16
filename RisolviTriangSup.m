%FUNZIONA

function [x]=RisolviTriangSup(U,y)
n=length(y);
x=zeros(n,1);
x(n,1)=y(n)/U(n,n);
for k=n-1:-1:1
    x(k,1)=(1/U(k,k))*(y(k)-U(k,k+1:n)*x(k+1:n));
end
end
    
    