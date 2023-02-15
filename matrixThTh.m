%Xi;
%N=15;
ThTh1=zeros(N^2);
for i=1:N
    for k=1:N
        for l=1:N
            ThTh1((i-1)*N+k,(i-1)*N+l)=G2(k,l);
        end
    end
end