%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BCCC,h/Ra=0.01
% clear;
% clc;
s1=zeros(70,1);
% for N=5:2:23
% for ii=0:4:20
%syms omega %��Ƶ��omega����Ϊ���ű�����������Ƶ��
%syms DeltaT %��������P����Ϊ���ű���������ٽ������غ�
% for ii=1:2:20
% syms DeltaT
% n=0;
% N=15;
xsdw=i;
% II=eye(N^2);
% Omega=100;
% q=0;
% Ra=0.2;
% parameter;
% h=0.01*Ra;
% Rb=0.4*Ra;
% Rb=0.1*Ra;
ii=1;
% Rb=0.05*ii*Ra
% % ri;
% % thi;
% % M1R;
% % M1Theta;
%DeltaT=0;
% kw=0*1e8;
% kg=0.9*1e5;
% kw=5e15;
% kg=500;
% Astiffness;
% Bstiffness;
% % Cstiffness;
% % matrixC1;
% % matrixC2;
% % matrixG1;
% % matrixG2;
% % matrixR;
% % matrixRR;
% % matrixTh;
% % matrixThTh;
% % matrixThR;
% % matrixLrLth;
% matrixK;
matrixK2dceil;%K2d�а����ṹ����ն������ת����ĸն���
% matrixKf;
% matrixKfceil;
% matrixKT;
% % matrixKTceil;
% vectorOm;
% % vectorOmceil;
% vectorQn;
matrixC;
matrixM;
matrixCOmega=C*Omega;
matrixkuai1=K2d+DeltaT*KT1;
matrixkuai2=matrixCOmega;
matrixkuai3=M;
for i=1:N %���߽������еķ���ϵ������
    for j=1:5*N^2
        matrixkuai1(i,j)=0;
        matrixkuai1((N-1)*N+i,j)=0;
        matrixkuai1(N^2+i,j)=0;
        matrixkuai1(N^2+(N-1)*N+i,j)=0;
        matrixkuai1(2*N^2+i,j)=0;
        matrixkuai1(2*N^2+(N-1)*N+i,j)=0;
        matrixkuai1(3*N^2+i,j)=0;
        matrixkuai1(3*N^2+(N-1)*N+i,j)=0;
        matrixkuai1(4*N^2+i,j)=0;
        matrixkuai1(4*N^2+(N-1)*N+i,j)=0;
    end
end
for j=1:N %���߽������еķ���ϵ�������������ղ�ֵ��������
    matrixkuai1(j,j)=1;
    matrixkuai1((N-1)*N+j,(N-1)*N+j)=1;
    matrixkuai1(N^2+j,N^2+j)=1;
    matrixkuai1(N^2+(N-1)*N+j,N^2+(N-1)*N+j)=1;
    matrixkuai1(2*N^2+j,2*N^2+j)=1;
    matrixkuai1(2*N^2+(N-1)*N+j,2*N^2+(N-1)*N+j)=1;
    matrixkuai1(3*N^2+j,3*N^2+j)=1;
    matrixkuai1(3*N^2+(N-1)*N+j,3*N^2+(N-1)*N+j)=1;
    matrixkuai1(4*N^2+j,4*N^2+j)=1;
    matrixkuai1(4*N^2+(N-1)*N+j,4*N^2+(N-1)*N+j)=1;
end
for i=1:N %���߽������еĿ���������صķ���ϵ������
    for j=1:5*N^2
        matrixkuai2(i,j)=0;
        matrixkuai2((N-1)*N+i,j)=0;
        matrixkuai2(N^2+i,j)=0;
        matrixkuai2(N^2+(N-1)*N+i,j)=0;
        matrixkuai2(2*N^2+i,j)=0;
        matrixkuai2(2*N^2+(N-1)*N+i,j)=0;
        matrixkuai2(3*N^2+i,j)=0;
        matrixkuai2(3*N^2+(N-1)*N+i,j)=0;
        matrixkuai2(4*N^2+i,j)=0;
        matrixkuai2(4*N^2+(N-1)*N+i,j)=0;
    end
end
for i=1:N %���߽������е���������صķ���ϵ������
    for j=1:5*N^2
        matrixkuai3(i,j)=0;
        matrixkuai3((N-1)*N+i,j)=0;
        matrixkuai3(N^2+i,j)=0;
        matrixkuai3(N^2+(N-1)*N+i,j)=0;
        matrixkuai3(2*N^2+i,j)=0;
        matrixkuai3(2*N^2+(N-1)*N+i,j)=0;
        matrixkuai3(3*N^2+i,j)=0;
        matrixkuai3(3*N^2+(N-1)*N+i,j)=0;
        matrixkuai3(4*N^2+i,j)=0;
        matrixkuai3(4*N^2+(N-1)*N+i,j)=0;
    end
end
%��״̬�����еľ������װ
matrixZTK1=zeros(10*N^2);
matrixZTK1(1:5*N^2,5*N^2+1:10*N^2)=-matrixkuai1;
matrixZTK1(5*N^2+1:10*N^2,1:5*N^2)=-matrixkuai1;
matrixZTK1(5*N^2+1:10*N^2,5*N^2+1:10*N^2)=-xsdw*matrixkuai2;
matrixZTK2=zeros(10*N^2);
matrixZTK2(1:5*N^2,1:5*N^2)=-matrixkuai1;
matrixZTK2(5*N^2+1:10*N^2,5*N^2+1:10*N^2)=-matrixkuai3;
[svector,s]=eig(-inv(matrixZTK1)*matrixZTK2);%��������ֵ������Ƶ�ʣ�������������ģ̬��
sdiag=diag(s);
jjj=1;
for jj=1:size(sdiag)
%     if ((real(sdiag(jj))>1e-10)||(real(sdiag(jj))<-1e-10))% &&(abs(imag(sdiag(jj))/real(sdiag(jj)))<2e-2))
%     if (real(sdiag(jj))>1e-10&&(abs(imag(sdiag(jj))/real(sdiag(jj)))<2e-2))
%     if real(sdiag(jj))>1e-10
    if real(sdiag(jj))>1e-10&&(abs(imag(sdiag(jj))/real(sdiag(jj)))<1e-1)
        s1(jjj,ii)=(1/sdiag(jj))/(2*pi);
        s1(jjj,ii+1)=jj;
%         s11(jjj,ii)=1/real(sdiag(jj));
%         s11(jjj,ii+1)=jj;
%         s2(jjj,(N+1)/2-2)=s1(jjj,(N+1)/2-2)*12*(1+0.34)*(8.2e-5)*(Ra/h).^2;
        jjj=jjj+1;
    end
end
% s=vpa(solve(det(matrixbc1+matrixbc2*DeltaT)))
% ss=double(s)
% s1(1:size(ss),ii)=ss;
% s2=vpa(s)*12*(1+0.34)*(8.2e-5)*(Ra/(2*h)).^2;
% xlswrite('critialBucklingTemperaturePCNchange.xlsx',s1);
% xlswrite('critialBucklingTemperaturePCNchangelambda.xlsx',s2);
% end

% [valueb,rowb]=max(dispwc);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��w�����ֵ
x=zeros(N^2,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������������ת��Ϊֱ������ϵ������
y=zeros(N^2,1);
for i=1:N
    for j=1:N
        x((i-1)*N+j,1)=r(1,i)*cos(th(1,j));
        y((i-1)*N+j,1)=r(1,i)*sin(th(1,j));
    end
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������������ת��Ϊֱ������ϵ������

% lnum=1140;
% disp6=zeros(N^2,5);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������д洢
% for i=1:N^2
%     disp6(i,1)=svector(i,lnum);
%     disp6(i,2)=svector(i+N^2,lnum);
%     disp6(i,3)=svector(i+2*N^2,lnum);
%     disp6(i,4)=svector(i+3*N^2,lnum);
%     disp6(i,5)=svector(i+4*N^2,lnum);    
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������д洢
% 
% lnum=1140;
% disp61=zeros(N^2,5);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������д洢
% for i=1:N^2
%     disp61(i,1)=sign(real(svector(i,lnum)))*norm(svector(i,lnum));
%     disp61(i,2)=sign(real(svector(i+N^2,lnum)))*norm(svector(i+N^2,lnum));
%     disp61(i,3)=sign(real(svector(i+2*N^2,lnum)))*norm(svector(i+2*N^2,lnum));
%     disp61(i,4)=sign(real(svector(i+3*N^2,lnum)))*norm(svector(i+3*N^2,lnum));
%     disp61(i,5)=sign(real(svector(i+4*N^2,lnum)))*norm(svector(i+4*N^2,lnum));    
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������д洢

% x=x';
% disp=zeros(N^2,5);
% for i=1:5
%     for j=1:N^2
%         disp(j,i)=ss((i-1)*N^2+j,1);
%     end
% end
% dispu=zeros(N); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������д洢
% for i=1:N
%     for j=1:N
%         dispu(i,j)=disp((i-1)*N+j,1);
%     end
% end
% dispv=zeros(N);
% for i=1:N
%     for j=1:N
%         dispv(i,j)=disp((i-1)*N+j,2);
%     end
% end
% dispw=zeros(N);
% for i=1:N
%     for j=1:N
%         dispw(i,j)=disp((i-1)*N+j,3);
%     end
% end
% dispphir=zeros(N);
% for i=1:N
%     for j=1:N
%         dispphir(i,j)=disp((i-1)*N+j,4);
%     end
% end
% dispphith=zeros(N);
% for i=1:N
%     for j=1:N
%         dispphith(i,j)=disp((i-1)*N+j,5);
%     end
% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������д洢
% x=zeros(N); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������������ת��Ϊֱ������ϵ������
% y=zeros(N);
% for i=1:N
%     for j=1:N
%         x(i,j)=r(1,i)*cos(th(1,j));
%         y(i,j)=r(1,i)*sin(th(1,j));
%     end
% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������������ת��Ϊֱ������ϵ������
% dispuc=zeros(N^2,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������д洢
% for i=1:N^2
%     dispuc(i,1)=disp(i,1);
% end
% dispvc=zeros(N^2,1);
% for i=1:N^2
%     dispvc(i,1)=disp(i,2);
% end
% dispwc=zeros(N^2,1);
% for i=1:N^2
%     dispwc(i,1)=disp(i,3);
% end
% dispphirc=zeros(N^2,1);
% for i=1:N^2
%     dispphirc(i,1)=disp(i,4);
% end
% dispphithc=zeros(N^2,1);
% for i=1:N^2
%     dispphithc(i,1)=disp(i,5);
% end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������д洢
% [value,row]=max(dispwc);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��w�����ֵ
% x=zeros(N^2,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������������ת��Ϊֱ������ϵ������
% y=zeros(N^2,1);
% for i=1:N
%     for j=1:N
%         x((i-1)*N+j,1)=r(1,i)*cos(th(1,j));
%         y((i-1)*N+j,1)=r(1,i)*sin(th(1,j));
%     end
% end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������������ת��Ϊֱ������ϵ������
% s=solve(det(matrixbc1+matrixbc2.*DeltaT));%��������P����Ϊ���ű���������ٽ������غ�
% %vpa(s);
% s2=vpa(s)*12*(1+0.34)*(8.2e-5)*(Ra/(2*h)).^2;%*Ra.^2/(1765e15*h.^3/(12*(1-0.3.^2)));
% %s1=(min(double(s2)))
% s1(1,ii)=(min(double(s2)));
% s1(1,ii)
% xlswrite('CCRaRbhRa0.01toCT.xlsx',s1);
% s1;
% end