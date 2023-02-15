%N=15;
%ri;
M1t=ones(1,N);
for i=1:N
    for k=1:N
        if k~=i
            M1t(1,i)=M1t(1,i)*(th(1,i)-th(1,k));
        else
            M1t(1,i)=M1t(1,i);
        end
    end
end