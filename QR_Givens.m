%G è una matrice 2x2, ortogonale 

%A matrice n x m qualsiasi

%FUNZIONA

function [Q,R]=QR_Givens(A)
[n,m]=size(A);
R=A;
U=eye(n);
for k=1:m
    %if [R(n-1,k),R(n,k)]==[0,0]
    %    rigadipartenza=n-1;
    %else
    %rigadipartenza=n;
    for i=n:-1:k+1
        x=R(i-1:i,k);
        normx=norm(x);
        if normx==0
            G=eye(2);
        else
            c=x(1)/normx;
            s=x(2)/normx;
            G=[c s;-s c];
        end
        R(i-1:i,k+1:m)=G*R(i-1:i,k+1:m);
        R(i-1,k)=normx;
        R(i,k)=0;
        U(i-1:i,1:n)=G*U(i-1:i,1:n);
    end
end
Q=U';
R=triu(R);

end
