%Fattorizzazione LU con Gauss
%Non ci sono HP su A

%FUNZIONA


function [L,U]=fattattorizzazione_LU_Gauss(A)
n=size(A,1);
U=A;
mol=zeros(n,n);
for k=1:n-1
    m=U(k+1:n,k)/U(k,k);
    mol(k+1:n,k)=m;
    U(k+1:n,k+1:n)=U(k+1:n,k+1:n)-m*U(k,k+1:n);
    
end
U=triu(U);
L=eye(n)+mol;
end