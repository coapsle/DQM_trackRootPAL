%Xi;
%matrixC1;
C2=zeros(N);
for i=1:N
    for j=1:N
        for k=1:N
            C2(i,j)=C2(i,j)+C1(i,k)*C1(k,j);
        end
    end
end
