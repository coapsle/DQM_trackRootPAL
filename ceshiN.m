%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BCCC,h/Ra=0.01
clear;
clc;
s1=zeros(1,14);
s2=zeros(1,14);
for N=31:2:51
%for ii=1:1:14
%syms omega %将频率omega假设为符号变量，求解固有频率
%syms DeltaT %将面内力P假设为符号变量，求解临界屈曲载荷
DeltaT=20;
n=3;
% N=19;
% h=0.005;
% ii=10;
Omega=10;
q=100;
parameter;
% laGPL=0.01;
% Rb=0.05*ii*Ra
ri;
thi;
M1R;
M1Theta;
%DeltaT=0;
kw=0;
kg=0;
%kw=5e15;
%kg=500;
Astiffness;
% Bstiffness;
% Cstiffness;
A22=A11;
B22=B11;
D22=D11;
matrixC1;
matrixC2;
matrixG1;
matrixG2;
matrixR;
matrixRR;
matrixTh;
matrixThTh;
matrixThR;
matrixLrLth;
% matrixK;
matrixKceil;
% matrixKf;
matrixKfceil;
% matrixKT;
matrixKTceil;
% vectorOm;
vectorOmceil;
vectorQn;
%matrixL1;
%matrixTemp;
%matrixKg;
%matrixPp;%计算临界屈曲温升不需要引入面内载荷
%matrixHada1;
%matrixHada2;%计算临界屈曲温升不需要引入惯性项
matrixbc1=K1+Kf1;
matrixbc2=KT1;
matrixbc3=(Om1+qn1);
for i=1:N %将边界条件中的方程系数归零
    for j=1:5*N^2
        matrixbc1(i,j)=0;
        matrixbc1((N-1)*N+i,j)=0;
        matrixbc1(N^2+i,j)=0;
        matrixbc1(N^2+(N-1)*N+i,j)=0;
        matrixbc1(2*N^2+i,j)=0;
        matrixbc1(2*N^2+(N-1)*N+i,j)=0;
        matrixbc1(3*N^2+i,j)=0;
        matrixbc1(3*N^2+(N-1)*N+i,j)=0;
        matrixbc1(4*N^2+i,j)=0;
        matrixbc1(4*N^2+(N-1)*N+i,j)=0;
    end
end
for j=1:N %将边界条件中的方程系数采用拉格朗日插值函数插入
    matrixbc1(j,j)=1;
    matrixbc1((N-1)*N+j,(N-1)*N+j)=1;
    matrixbc1(N^2+j,N^2+j)=1;
    matrixbc1(N^2+(N-1)*N+j,N^2+(N-1)*N+j)=1;
    matrixbc1(2*N^2+j,2*N^2+j)=1;
    matrixbc1(2*N^2+(N-1)*N+j,2*N^2+(N-1)*N+j)=1;
    matrixbc1(3*N^2+j,3*N^2+j)=1;
    matrixbc1(3*N^2+(N-1)*N+j,3*N^2+(N-1)*N+j)=1;
    matrixbc1(4*N^2+j,4*N^2+j)=1;
    matrixbc1(4*N^2+(N-1)*N+j,4*N^2+(N-1)*N+j)=1;
end
for i=1:N %将边界条件中的温度相关的方程系数归零
    for j=1:5*N^2
        matrixbc2(i,j)=0;
        matrixbc2((N-1)*N+i,j)=0;
        matrixbc2(N^2+i,j)=0;
        matrixbc2(N^2+(N-1)*N+i,j)=0;
        matrixbc2(2*N^2+i,j)=0;
        matrixbc2(2*N^2+(N-1)*N+i,j)=0;
        matrixbc2(3*N^2+i,j)=0;
        matrixbc2(3*N^2+(N-1)*N+i,j)=0;
        matrixbc2(4*N^2+i,j)=0;
        matrixbc2(4*N^2+(N-1)*N+i,j)=0;
    end
end
for i=1:N
    matrixbc3(i,1)=0;
    matrixbc3((N-1)*N+i,1)=0;
    matrixbc3(N^2+i,1)=0;
    matrixbc3(N^2+(N-1)*N+i,1)=0;
    matrixbc3(2*N^2+i,1)=0;
    matrixbc3(2*N^2+(N-1)*N+i,1)=0;
    matrixbc3(3*N^2+i,1)=0;
    matrixbc3(3*N^2+(N-1)*N+i,1)=0;
    matrixbc3(4*N^2+i,1)=0;
    matrixbc3(4*N^2+(N-1)*N+i,1)=0;
end
% x=inv(matrixbc1+matrixbc2*DeltaT)*(Om1+qn1);
ss=(matrixbc1+matrixbc2*DeltaT)\(matrixbc3);
% ss=(sym(matrixbc1+matrixbc2*DeltaT))\(sym(matrixbc3));
% x=x';
disp=zeros(N^2,5);
for i=1:5
    for j=1:N^2
        disp(j,i)=ss((i-1)*N^2+j,1);
    end
end
% dispu=zeros(N); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果分行存储
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
% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果分行存储
% x=zeros(N); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将极坐标坐标轴转换为直角坐标系下坐标
% y=zeros(N);
% for i=1:N
%     for j=1:N
%         x(i,j)=r(1,i)*cos(th(1,j));
%         y(i,j)=r(1,i)*sin(th(1,j));
%     end
% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将极坐标坐标轴转换为直角坐标系下坐标
dispuc=zeros(N^2,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果分列存储
for i=1:N^2
    dispuc(i,1)=disp(i,1);
end
dispvc=zeros(N^2,1);
for i=1:N^2
    dispvc(i,1)=disp(i,2);
end
dispwc=zeros(N^2,1);
for i=1:N^2
    dispwc(i,1)=disp(i,3);
end
dispphirc=zeros(N^2,1);
for i=1:N^2
    dispphirc(i,1)=disp(i,4);
end
dispphithc=zeros(N^2,1);
for i=1:N^2
    dispphithc(i,1)=disp(i,5);
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果分列存储
[s1(1,(N-29)/2),s2(1,(N-29)/2)]=max(dispwc);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%求w中最大值
x=zeros(N^2,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将极坐标坐标轴转换为直角坐标系下坐标
y=zeros(N^2,1);
for i=1:N
    for j=1:N
        x((i-1)*N+j,1)=r(1,i)*cos(th(1,j));
        y((i-1)*N+j,1)=r(1,i)*sin(th(1,j));
    end
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将极坐标坐标轴转换为直角坐标系下坐标
% s=solve(det(matrixbc1+matrixbc2.*DeltaT));%将面内力P假设为符号变量，求解临界屈曲载荷
% %vpa(s);
% s2=vpa(s)*12*(1+0.34)*(8.2e-5)*(Ra/(2*h)).^2;%*Ra.^2/(1765e15*h.^3/(12*(1-0.3.^2)));
% %s1=(min(double(s2)))
% s1(1,ii)=(min(double(s2)));
% s1(1,ii)
% xlswrite('CCRaRbhRa0.01toCT.xlsx',s1);
% s1;
 end