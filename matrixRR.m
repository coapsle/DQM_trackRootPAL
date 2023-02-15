%Xi;
%N=15;
RR1=zeros(N^2);
for i=1:N
    for j=1:N
        for k=1:N
            RR1((i-1)*N+k,(j-1)*N+k)=C2(i,j);
        end
    end
end