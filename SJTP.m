function f = SJTP(inputArg1,inputArg2)%------ ˵����inputArg1 Ϊm*1�������� �� inputArg2 Ϊm*n�ľ���-----
lineN=size(inputArg1,1);
f=zeros(lineN);
for i=1:1:lineN
    for j=1:1:lineN
        f(i,j)=inputArg1(i,1).*inputArg2(i,j);
    end
end 
end
