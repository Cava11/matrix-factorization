%Fattorizzazione QR di Givens per le hessemberg inferiori

%FUNZIONA

function [Q,L]=QR_HessembergInferiore_Givens(A)
A=A';
[Q,R]=QR_HessembergSuperiore_Givens(A);
Q=Q';
L=R';
end
