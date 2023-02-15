%N=15;
%ri;
M1r=ones(1,N);
for i=1:N
    for k=1:N
        if k~=i
            M1r(1,i)=M1r(1,i)*(r(1,i)-r(1,k));
        else
            M1r(1,i)=M1r(1,i);
        end
    end
end