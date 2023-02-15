err=1e-16;
% x0=dispc;
tol=1;
nn=1
while tol>err
%     r = x0 - myfun(x0)/dmyfun(x0);
    JacobiYuan;
    jisuanYuanEqsJg;
    rr=dispc-JacobiY\(YuanEqsJg-matrixbc3);
    tol=norm(rr-dispc);
    nn=nn+1
    dispc=rr;
    disp6=zeros(N^2,5);
for i=1:5
    for j=1:N^2
        disp6(j,i)=dispc((i-1)*N^2+j,1);
    end
end
dispuc=zeros(N^2,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果分列存储
for i=1:N^2
    dispuc(i,1)=disp6(i,1);
end
dispvc=zeros(N^2,1);
for i=1:N^2
    dispvc(i,1)=disp6(i,2);
end
dispwc=zeros(N^2,1);
for i=1:N^2
    dispwc(i,1)=disp6(i,3);
end
dispphirc=zeros(N^2,1);
for i=1:N^2
    dispphirc(i,1)=disp6(i,4);
end
dispphithc=zeros(N^2,1);
for i=1:N^2
    dispphithc(i,1)=disp6(i,5);
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果分列存储
    if(nn>500)
        disp('迭代步数太多，可能不收敛！');
        return;
    end
end

ddispcha=disp5-disp6;
ddisperfanshu=norm(ddispcha);
ddispzuidawucha=max(ddispcha);
ddispcha3=disp5(:,3)-disp6(:,3);
ddisp3erfanshu=norm(ddispcha3);