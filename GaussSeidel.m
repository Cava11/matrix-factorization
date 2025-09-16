%Approssimazione iterativa A=P-N
%Splitting A=-E+D-F
%E:=triangolare strettamente inferiore
%D:=diagonale di A (non singolare)
%F:=triangolare strettamente superiore

%P=D-E
%N:=F

%B_GS=P^-1 * N

%A deve essere a diagonale strettamente dominante

%FUNZIONA

function [x,its,normr] = GaussSeidel(A,f,x0,maxit,tol)
r = f-A*x0;
x = x0;
normr0=norm(r);
normr(1)=normr0;
P=tril(A);
its=1;

while its<maxit
    x = x + P\r; % Usare propria funzione di eliminaz.Gauss
    r = f - A*x;
    its = its+1;
    normr(its,1)=norm(r)/normr0;
    %disp([its,normr(its)])
    if norm(r)/normr0 < tol
        break,end
    
end