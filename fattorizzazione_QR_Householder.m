%Fattorizzazione QR con Householder

%Per matrici rettangolari nxm con n>=m
%Si avrà P_m*...P_1*A=R; P_i sono le matrici di trasformazioni di
%Householder

%Le P sono simmetriche,ortogonali e possono annullare componenti di un
%vettore

%FUNZIONA

function [Q,R]=fattorizzazione_QR_Householder(A)
[n,m]=size(A);

R=A;
U=eye(n,n);

for k=1:m
    x=R(k:n,k);
    alfa=-sign(x(1))*norm(x);
    e1=eye(n-k+1,1);
    v=x-alfa*e1;
    
    beta=2/(v'*v);
    R(k:n,k:m)=R(k:n,k:m)-beta*v*(v'*R(k:n,k:m));
    
    %Ora aggiorniamo il prodotto delle amtrici ortogonali
    U(k:n,1:n)=U(k:n,1:n)-beta*v*(v'*U(k:n,1:n));   
end

%U è ortogonale per costruzione
%La U mi dà UA=R quindi per avere A=QR impongo Q=U'
Q=U';
end
    
    
    
