%Sia A ? R n×n di tipo Hessenberg superiore che ammette fattorizzazione LU

%FATTORIZZAZIONE LU

%FUNZIONA

function [L,U]=fattorizzazione_HessSup(A)
%determinare la fatt LU di una matrice hess. sup. , L sarÃ  bidiag
%implemento gauss
n=length(A);
U=A;
mol=zeros(n,n);
for k=1:n-1
    m=U(k+1,k)/U(k,k);
    mol(k+1,k)=m;
    U(k+1,k+1:n)=U(k+1,k+1:n)-m*U(k,k+1:n);
end
U=triu(U);
L=eye(n)+mol;
end
    
