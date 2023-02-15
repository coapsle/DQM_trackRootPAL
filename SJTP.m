function f = SJTP(inputArg1,inputArg2)%------ 说明：inputArg1 为m*1的列向量 和 inputArg2 为m*n的矩阵-----
lineN=size(inputArg1,1);
f=zeros(lineN);
for i=1:1:lineN
    for j=1:1:lineN
        f(i,j)=inputArg1(i,1).*inputArg2(i,j);
    end
end 
end
