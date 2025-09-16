%FUNZIONA

function [L] = fattorizzazione_LU_Cholesky(A)
%A deve essere simmetrica e definita positiva
%fattorizzo A=L*L' con L triang inf, non singolare e L(i,i)>0
%U=L'

L=zeros(size(A));
n=length(A);
L(1,1) = sqrt(A(1,1));
for j = 2:n
    for i = 1:j-1
        L(i,j) = (1/L(i,i))*(A(i,j)-sum(L(1:i-1,i).*L(1:i-1,j)));
    end
    L(j,j) = sqrt(A(j,j) - sum(L(1:j-1,j).^2));
end
L = L';
end