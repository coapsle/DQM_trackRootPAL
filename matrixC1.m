%Xi;
%N=15;
C1=ones(N);
for i=1:N
    for j=1:N
        if j~=i
            C1(i,j)=M1r(1,i)/((r(1,i)-r(1,j))*M1r(1,j));
        else
            C1(i,j)=0;
        end
    end
end
for i=1:N
    for j=1:N
        if j~=i
            C1(i,i)=C1(i,i)-C1(i,j);
        end
    end
end