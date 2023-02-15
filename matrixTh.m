%Xi;
%N=15;
Th1=zeros(N^2);
for i=1:N
    for k=1:N
        for l=1:N
            Th1((i-1)*N+k,(i-1)*N+l)=G1(k,l);
        end
    end
end