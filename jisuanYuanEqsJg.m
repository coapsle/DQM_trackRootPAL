%计算原非线性方程组的结果F(x^k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %计算原非线性方程组中的等号左边的非线性项结果
    %赋值非线性方程中非线性项结果
matrixNLI=zeros(5*N^2,1);
% 建立ri矩阵
% II=eye(N^2);
% ri1=zeros(N^2,1);
% for i=1:1:N^2
%     ri1(i,1)=1/r(1,ceil(i/N));
% end
% ri2=zeros(N^2,1);
% for i=1:1:N^2
%     ri2(i,1)=1/(r(1,ceil(i/N))^2);
% end
% ri3=zeros(N^2,1);
% for i=1:1:N^2
%     ri3(i,1)=1/(r(1,ceil(i/N))^3);
% end
% ri4=zeros(N^2,1);
% for i=1:1:N^2
%     ri4(i,1)=1/(r(1,ceil(i/N))^4);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 赋值（1,1）*(N^2*N^2)非线性方程中非线性项结果
matrixNLI(1:N^2,1)=A11*((R1*dispwc).*(RR1*dispwc)+1/2*ri1.*(R1*dispwc).*(R1*dispwc))...
    +A12*(-1/2*ri3.*(Th1*dispwc).*(Th1*dispwc)+ri2.*(Th1*dispwc).*(ThR*dispwc)-1/2*ri1.*(R1*dispwc).*(R1*dispwc))...
    -A22*1/2*ri3.*(Th1*dispwc).*(Th1*dispwc)...
    +A44*(ri2.*(Th1*dispwc).*(ThR*dispwc)+ri2.*(R1*dispwc).*(ThTh1*dispwc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 赋值（2,1）*(N^2,1)非线性方程中非线性项结果
matrixNLI(N^2+1:2*N^2,1)=A12*ri1.*(R1*dispwc).*(ThR*dispwc)...
    +A12*ri3.*(Th1*dispwc).*(ThTh1*dispwc)...
    +A44*(ri2.*(R1*dispwc).*(Th1*dispwc)+ri1.*(Th1*dispwc).*(RR1*dispwc)+ri1.*(R1*dispwc).*(ThR*dispwc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 赋值（3,1）*(N^2,1)非线性方程中非线性项结果
matrixNLI(2*N^2+1:3*N^2,1)=A11*(ri1.*(R1*dispuc).*(R1*dispwc)+1/2*ri1.*(R1*dispwc).*(R1*dispwc).*(R1*dispwc)...
    +(RR1*dispuc).*(R1*dispwc)+(R1*dispuc).*(RR1*dispwc)+3/2*(R1*dispwc).*(R1*dispwc).*(RR1*dispwc))...
    +A12*(ri1.*(ThR*dispvc).*(R1*dispwc)+ri1.*(R1*dispuc).*(R1*dispwc)...
    -1/2*ri3.*(R1*dispwc).*(Th1*dispwc).*(Th1*dispwc)+2*ri2.*(R1*dispwc).*(Th1*dispwc).*(ThR*dispwc)...
    +ri1.*(Th1*dispvc).*(RR1*dispwc)+ri1.*(II*dispuc).*(RR1*dispwc)...
    +1/2*ri2.*(RR1*dispwc).*(Th1*dispwc).*(Th1*dispwc)...
    +ri2.*(ThR*dispuc).*(Th1*dispwc)+ri2.*(R1*dispuc).*(ThTh1*dispwc)...
    +1/2*ri2.*(R1*dispwc).*(R1*dispwc).*(ThTh1*dispwc))...
    +A22*(ri3.*(ThTh1*dispvc).*(Th1*dispwc)+ri3.*(Th1*dispuc).*(Th1*dispwc)...
    +ri3.*(Th1*dispvc).*(ThTh1*dispwc)+ri3.*(II*dispuc).*(ThTh1*dispwc)...
    +3/2*ri4.*(Th1*dispwc).*(Th1*dispwc).*(ThTh1*dispwc))...
    +A44*(-ri3.*(Th1*dispuc).*(Th1*dispwc)+ri2.*(ThR*dispuc).*(Th1*dispwc)+ri1.*(RR1*dispvc).*(Th1*dispwc)...
    +ri3.*(II*dispvc).*(Th1*dispwc)-ri2.*(R1*dispvc).*(Th1*dispwc)...
    -ri3.*(R1*dispwc).*(Th1*dispwc).*(Th1*dispwc)+ri2.*(RR1*dispwc).*(Th1*dispwc).*(Th1*dispwc)...
    +ri2.*(ThTh1*dispuc).*(R1*dispwc)+ri1.*(ThR*dispvc).*(R1*dispwc)-ri2.*(Th1*dispvc).*(R1*dispwc)...
    +ri2.*(R1*dispwc).*(R1*dispwc).*(ThTh1*dispwc)+2*ri2.*(Th1*dispuc).*(ThR*dispwc)...
    +2*ri1.*(R1*dispvc).*(ThR*dispwc)-2*ri2.*(II*dispvc).*(ThR*dispwc)...
    +4*ri2.*(R1*dispwc).*(Th1*dispwc).*(ThR*dispwc))...
    +B11*(ri1.*(R1*dispphirc).*(R1*dispwc)+(RR1*dispphirc).*(R1*dispwc)+(R1*dispphirc).*(RR1*dispwc))...
    +B12*(ri1.*(R1*dispphirc).*(R1*dispwc)+ri1.*(ThR*dispphithc).*(R1*dispwc)...
    +ri1.*(II*dispphirc).*(RR1*dispwc)+ri1.*(Th1*dispphithc).*(RR1*dispwc)...
    +ri2.*(ThR*dispphirc).*(Th1*dispwc)+ri2.*(R1*dispphirc).*(ThTh1*dispwc))...
    +B22*(ri3.*(Th1*dispphirc).*(Th1*dispwc)+ri3.*(ThTh1*dispphithc).*(Th1*dispwc)...
    +ri3.*(II*dispphirc).*(ThTh1*dispwc)+ri3.*(Th1*dispphithc).*(ThTh1*dispwc))...
    +B44*(-ri3.*(Th1*dispphirc).*(Th1*dispwc)+ri2.*(ThR*dispphirc).*(Th1*dispwc)...
    +ri1.*(RR1*dispphithc).*(Th1*dispwc)+ri3.*(II*dispphithc).*(Th1*dispwc)...
    -ri2.*(R1*dispphithc).*(Th1*dispwc)+ri2.*(ThTh1*dispphirc).*(R1*dispwc)...
    +ri1.*(ThR*dispphithc).*(R1*dispwc)-ri2.*(Th1*dispphithc).*(R1*dispwc)...
    +2*ri2.*(Th1*dispphirc).*(ThR*dispwc)+2*ri1.*(R1*dispphithc).*(ThR*dispwc)...
    -2*ri2.*(II*dispphithc).*(ThR*dispwc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 赋值（4,1）*(N^2,1)非线性方程中非线性项结果
matrixNLI(3*N^2+1:4*N^2,1)=B11*((R1*dispwc).*(RR1*dispwc)+1/2*ri1.*(R1*dispwc).*(R1*dispwc))...
    +B12*(-1/2*ri3.*(Th1*dispwc).*(Th1*dispwc)+ri2.*(Th1*dispwc).*(ThR*dispwc)-1/2*ri1.*(R1*dispwc).*(R1*dispwc))...
    -B22*1/2*ri3.*(Th1*dispwc).*(Th1*dispwc)...
    +B44*(ri2.*(Th1*dispwc).*(ThR*dispwc)+ri2.*(R1*dispwc).*(ThTh1*dispwc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 赋值（5,1）*(N^2,1)非线性方程中非线性项结果
matrixNLI(4*N^2+1:5*N^2,1)=B12*ri1.*(R1*dispwc).*(ThR*dispwc)...
    +B22*ri3.*(Th1*dispwc).*(ThTh1*dispwc)...
    +B44*(ri2.*(Th1*dispwc).*(R1*dispwc)+ri1.*(Th1*dispwc).*(RR1*dispwc)+ri1.*(R1*dispwc).*(ThR*dispwc));
    %赋值非线性方程中非线性项结果
    %计算原非线性方程组中的等号左边的非线性项结果
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %计算原非线性方程组中的等号左边的线性项结果
    %赋值非线性方程中线性项结果
matrixKK=matrixbc1+matrixbc2*DeltaT;
matrixLiAns=matrixKK*dispc;
    %赋值非线性方程中线性项结果
    %计算原非线性方程组中的等号左边的线性项结果
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %将原非线性方程组中的等号左边的非线性项结果与线性项结果相加
YuanEqsJg=matrixLiAns+matrixNLI;
    %将原非线性方程组中的等号左边的非线性项结果与线性项结果相加
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %将边界条件赋值到原非线性方程组中的等号左边
YuanEqsJg(1:N,1)=dispc(1:N,1);
YuanEqsJg((N-1)*N+1:N^2,1)=dispc((N-1)*N+1:N^2,1);
YuanEqsJg(N^2+1:N^2+N,1)=dispc(N^2+1:N^2+N,1);
YuanEqsJg(N^2+(N-1)*N+1:2*N^2,1)=dispc(N^2+(N-1)*N+1:2*N^2,1);
YuanEqsJg(2*N^2+1:2*N^2+N,1)=dispc(2*N^2+1:2*N^2+N,1);
YuanEqsJg(2*N^2+(N-1)*N+1:3*N^2,1)=dispc(2*N^2+(N-1)*N+1:3*N^2,1);
YuanEqsJg(3*N^2+1:3*N^2+N,1)=dispc(3*N^2+1:3*N^2+N,1);
YuanEqsJg(3*N^2+(N-1)*N+1:4*N^2,1)=dispc(3*N^2+(N-1)*N+1:4*N^2,1);
YuanEqsJg(4*N^2+1:4*N^2+N,1)=dispc(4*N^2+1:4*N^2+N,1);
YuanEqsJg(4*N^2+(N-1)*N+1:5*N^2,1)=dispc(4*N^2+(N-1)*N+1:5*N^2,1);
    %将边界条件赋值到原非线性方程组中的等号左边
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

