%Xi;
%N=15;
R1=zeros(N^2);
for i=1:N
    for j=1:N
        for k=1:N
            R1((i-1)*N+k,(j-1)*N+k)=C1(i,j);
        end
    end
end