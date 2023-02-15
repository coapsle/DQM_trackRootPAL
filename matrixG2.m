%Xi;
%matrixC1;
G2=zeros(N);
for i=1:N
    for j=1:N
        for k=1:N
            G2(i,j)=G2(i,j)+G1(i,k)*G1(k,j);
        end
    end
end
