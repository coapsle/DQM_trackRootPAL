function x = utri(U,y)
n=size(y,1);
x=zeros(n,1);
for j = n:-1:2
    x(j) = y(j)/U(j,j);
    y(1:j-1) = y(1:j-1) - x(j) * U(1:j-1,j);
end
x(1) = y(1)/U(1,1);