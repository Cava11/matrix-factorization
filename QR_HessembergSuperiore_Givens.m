%Scomposizione QR_HessembergSuperiore

%Con Givens

%Algoritmo di base QR_Givens

%FUNZIONA

%Restituisce la fattorizzazione QR con Q ortogonale e R triangolare
%superiore
function [Q,R]=QR_HessembergSuperiore_Givens(A)
n=size(A,1);

R=A;
U=eye(n);
for k=1:n-1
    x=R(k:k+1,k);
    normx=norm(x);
    if normx==0
        G=eye(2);
    else
        c=x(1)/normx;
        s=x(2)/normx;
        G=[c s;-s c];
    end
    R(k:k+1,k:n)=G*R(k:k+1,k:n);
    R(k:k+1,k)=normx;
    U(k:k+1,1:n)=G*U(k:k+1,1:n);

end
Q=U';
R=triu(R);

end
