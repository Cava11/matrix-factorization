%Sherman Morrison
%Per A ed f del precedente esercizio, e dato il sistema
%Bu = f, dove B `e uguale ad A, tranne la prima colonna,
%che `e data da B(i, 1) = (n?1)2/(4?2), i=1,...,n.
%i) Determina esplicitamente w, v 
%in modo che valga B = A + wvT.
%ii) Per n = 100, scrivi uno script per una procedura che risolva il sistema
%Bu = f, sfruttando in modo opportuno la seguente formula di Sherman-Morrison:
%(A+wv')^(?1 )=A(?1)?A(?1)w(1+v'A(?1)w)(?1)v'A(?1),
%in modo da dover solo  risolvere sistemi lineari con A, mediante
%l’algoritmo di Thomas.

%utile per la risoluzione di (A+uv')x=b se so risolvere Ax=b

%FUNZIONA 

function [x,normSM]=Sherman_Morrison(L,u,v,b) %caso L qualsiasi
%risolvere Ax=b con A=L+uv' variazione di rango uno di una non singolare
% Risolvi i due sistemi
w2=L\b; %L^-1 * b, 'implementa con il codice che preferichi per migliorare
w1=L\u; %L^-1 * u,   il costo dell' algoritmo'
theta1=1+v'*w1;
theta2=v'*w2;
% applica la formula di SM
x = w2-w1*(theta2/theta1); %x=L^-1*b
% controlla l'accuratezza senza creare B esplicitamente
normSM=norm( b - L*x - u*(v'*x));


