dispZTc(1:5*N^2,1)=dispc;%��λ�ƴ�5*N^2����չ��5*N^2+1��
dispZTc(5*N^2+1,1)=DeltaT;%��λ�ƴ�5*N^2����չ��5*N^2+1��
Om2(1:5*N^2,1)=Om1;%���Ⱥ��ұ߽����5*N^2����չ��5*N^2+1��
Om2(5*N^2+1,1)=0;%���Ⱥ��ұ߽����5*N^2����չ��5*N^2+1��
DeltaS=0.1;%���廡�������صĲ���DeltaS
err=1e-12;
% x0=dispc;
tol=1;
nn=1
while tol>err
%     r = x0 - myfun(x0)/dmyfun(x0);
    JacobiZhengti;
    jisuanZhengtiEqsJg;
    rrzt=dispZTc-JacobiZT\(ZhengtiEqsJg-Om2);
    tol=norm(rrzt-dispZTc);
    nn=nn+1
    dispZTc=rrzt;
    disp7=zeros(N^2,5);
    DeltaT=dispZTc(5*N^2+1,1);
for i=1:5
    for j=1:N^2
        disp7(j,i)=dispZTc((i-1)*N^2+j,1);
    end
end
dispuc=zeros(N^2,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������д洢
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
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������д洢
    if(nn>500)
        disp('��������̫�࣬���ܲ�������');
        return;
    end
end

ddispcha=disp5-disp7;
ddisperfanshu=norm(ddispcha)
ddispzuidawucha=max(ddispcha)
ddispcha3=disp5(:,3)-disp7(:,3);
ddisp3erfanshu=norm(ddispcha3)