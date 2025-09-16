%QR per tutte, anche rettangolari n>=m

%FUNZIONA

function [Q1,R1]=fattorizzo_QR_GramShmidt(A)
[n,m]=size(A);
Q1=zeros(n,m);
R1=zeros(m,m);
R1(1,1)=norm(A(1:n,1));
Q1(1:n,1)=A(1:n,1)/R1(1,1);

for k=2:m
    R1(1:(k-1),k)=Q1(1:n,1:k-1)'*A(1:n,k);
    q=A(1:n,k)-Q1(1:n,1:k-1)*R1(1:(k-1),k);
    R1(k,k)=norm(q);
    if R1(k,k)<1e-14
        fprintf('linear dependence\n'); break
    else
        Q1(1:n,k)=q/R1(k,k);
    end
end
end

    