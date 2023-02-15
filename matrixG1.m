%Xi;
%N=15;
G1=ones(N);
for i=1:N
    for j=1:N
        if j~=i
            G1(i,j)=M1t(1,i)/((th(1,i)-th(1,j))*M1t(1,j));
        else
            G1(i,j)=0;
        end
    end
end
for i=1:N
    for j=1:N
        if j~=i
            G1(i,i)=G1(i,i)-G1(i,j);
        end
    end
end