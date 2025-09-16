%Metodo iterativo di Jacobi

%Siamo interessati a risolvere Ax=f, A in R^(nxn)

%rk=residuo all'iterazione k

%Di solito consideriamo: 
%x0=zeros(n,1)
%maxit=5000
%tol=1e-8

function [x,its,resnormvec]=Jacobi(A,f,x0,maxit, tol)
% x=vettore soluzione
% its=numero di iterazioni effettuate
% resnormvec=contiene norm(residuo_k) per ogni k=0,_,its

%Consideriamo lo splitting A=P-N, scegliamo x0 e calcoliamo r0

%A=-E+D-F, E stret_tri_inf, F stret_tri_sup, D diag. --> P=D , N=E+F (tutte n*n)

%metodo iterativo stazionario, splitting
%A=P-N -> Ax=Px-Nx=b -> Px=Nx+b, P non singolare

x=x0;

r=f-A*x;

its=0;
resnormvec=[norm(r),];

P=diag(A);

for k=1:maxit
    x=x+r./P;
    r=f-A*x;
    resnormvec(k)=norm(r);
    if resnormvec(k)<tol
        break, end
    its=its+1;
end
end



