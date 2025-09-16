function [u]=RisolviThomas(A,f)
n=length(f);
L=zeros(n,n);
U=L;
L=L+diag(ones(n,1));
a=diag(A);
b=diag(A,-1);
c=diag(A,1);
alpha=zeros(n,1);
beta=zeros(n,1);

alpha(1)=a(1);

for i=2:n
    beta(i)=b(i-1)/alpha(i-1);
    alpha(i)=a(i)-beta(i)*c(i-1);
end

y=zeros(n,1);
y(1)=f(1);

for i=2:n
    y(i)=f(i)-beta(i)*y(i-1);
end

u=zeros(n,1);
u(n)=y(n)/alpha(n);

for i=n-1:-1:1
    u(i)=(y(i)-c(i)*u(i+1))/alpha(i);
end

    
    
    
    
    
    
    



