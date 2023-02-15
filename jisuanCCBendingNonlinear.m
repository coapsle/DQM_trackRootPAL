err=1e-9;
% x0=dispc;
tol=1;
nn=1
while tol>err
%     r = x0 - myfun(x0)/dmyfun(x0);
    JacobiYuan;
    jisuanYuanEqsJg;
    rr=dispc-JacobiY\(YuanEqsJg-Om1);
    tol=norm(rr-dispc);
    nn=nn+1;
    if mod(nn,20)==0
        nn
    end
    dispc=rr;
    disp7=zeros(N^2,5);
    [valueBI,rowBI]=max(dispwc);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%求w中最大值
for i=1:5
    for j=1:N^2
        disp7(j,i)=dispc((i-1)*N^2+j,1);
    end
end
dispuc=zeros(N^2,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果分列存储
for i=1:N^2
    dispuc(i,1)=disp7(i,1);
end
dispvc=zeros(N^2,1);
for i=1:N^2
    dispvc(i,1)=disp7(i,2);
end
dispwc=zeros(N^2,1);
for i=1:N^2
    dispwc(i,1)=disp7(i,3);
end
dispphirc=zeros(N^2,1);
for i=1:N^2
    dispphirc(i,1)=disp7(i,4);
end
dispphithc=zeros(N^2,1);
for i=1:N^2
    dispphithc(i,1)=disp7(i,5);
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果分列存储
    if((max(dispwc)/h)>50)
        disp('最大位移过大！');
        return;
    end
    if (rcond(JacobiY)<1e-30)
        disp('矩阵接近奇异，可能不收敛！');
        return
    end
    if(nn>50)
        disp('迭代步数太多，可能不收敛！');
        return;
    end
end

for i=1:N^2
    dispvc(i,1)=0;
end

figure
plot3(x,y,dispuc)
% figure
% plot3(x,y,dispvc)
% figure
% plot3(x,y,dispwc)
% figure
% plot3(x,y,dispphirc)
% figure
% plot3(x,y,dispphithc)

% figure
% plot3(x,y,dispuc)
% 
% figure
% plot3(x,y,dispphirc)
% figure%绘制w位移图
% plot3(x,y,dispwc)%绘制w位移图
% ddispcha=disp5-disp7;
% ddisperfanshu=norm(ddispcha);
% ddispzuidawucha=max(ddispcha);
% ddispcha3=disp5(:,3)-disp7(:,3);
% ddisp3erfanshu=norm(ddispcha3);