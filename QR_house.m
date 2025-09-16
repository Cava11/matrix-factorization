%QR per matrici quadrate o rettangolari con n>=m; dove abbiamo:
%Q=P_m*___*P_1
%P:=I-beta v v' con beta:=2/norm(v)^2
%P è una modifica di rango una della matrice Identità, P è simmetrica e
%ortogonale, è una riflessione, può annullare componenti di un vettore
%ad ogni iterazione k si costruisce una P che annulla la colonna k fino
%alla riga n-k+1

%FUNZIONA

function [Q,R]=QR_house(A)
[n,m]=size(A);
R=A;
U=eye(n);
for k=1:m
    x=R(k:n,k);
    alpha = - sign(x(1))*norm(x);
    e1=eye(n-k+1,1);
    v=x-alpha*e1;
    beta = 2/(v'*v);
    R(k:n,k:m)=R(k:n,k:m) - beta*v*((v'*R(k:n,k:m)));
    U(k:n,1:n)=U(k:n,1:n) - beta*v*((v'*U(k:n,1:n)));
end
Q=U';
end