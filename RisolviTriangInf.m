
function [x] = RisolviTriangInf(L,y)
%dato il sistema Ux=y ritrovare x
n=length(y);
x(1,1)=y(1,1)/L(1,1);
for k = 2:n
    x(k,1)=(y(k,1)-L(k,1:k-1)*x(1:k-1,1))/L(k,k);
end
end