%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%赋值x*DeltaT*
dispcxing=dispc;
DeltaTxing=DeltaT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%赋值x*pian
    %赋值原方程组Jacobian矩阵中x相关项
JacobiYuan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %赋值原方程组Jacobian矩阵中DeltaT相关项
JacobiDeltaT=zeros(5*N^2,1);
    %赋值原方程组Jacobian矩阵中DeltaT相关项中（3,1）*(N^2*1)矩阵块
JacobiDeltaT(2*N^2+1:3*N^2,1)=-A11alpha*(ri1.*(R1*dispwc)+RR1*dispwc)...
    -A12alpha*(RR1*dispwc+ri1.*(R1*dispwc)+ri2.*(ThTh1*dispwc))...
    -A22alpha*(ri2.*(ThTh1*dispwc));
    %将边界条件赋值到原方程组Jacobian矩阵中DeltaT相关项中（3,1）*(N^2*1)矩阵块中
JacobiDeltaT(2*N^2+1:2*N^2+N,1)=0;
JacobiDeltaT(2*N^2+(N-1)*N+1:3*N^2,1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matrixZ=-JacobiY\JacobiDeltaT;%计算列向量Z
% if dispwc(rowBI)>=0
%     JacobiAlphaPian=-(1+matrixZ.'*matrixZ)^(1/2);%计算附加方程对DeltaT的偏导alpha'
% else
%     JacobiAlphaPian=(1+matrixZ.'*matrixZ)^(1/2);%计算附加方程对DeltaT的偏导alpha'
% end
JacobiAlphaPian=-(1+matrixZ.'*matrixZ)^(1/2);%计算附加方程对DeltaT的偏导alpha'
% if matrixZlao'*matrixZ<0
%     JacobiAlphaPian=(1+matrixZ.'*matrixZ)^(1/2)
% end
JacobiXpian=matrixZ*JacobiAlphaPian;%计算附加方程对x的偏导x'