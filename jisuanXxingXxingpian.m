%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ֵx*DeltaT*
dispcxing=dispc;
DeltaTxing=DeltaT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ֵx*pian
    %��ֵԭ������Jacobian������x�����
JacobiYuan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %��ֵԭ������Jacobian������DeltaT�����
JacobiDeltaT=zeros(5*N^2,1);
    %��ֵԭ������Jacobian������DeltaT������У�3,1��*(N^2*1)�����
JacobiDeltaT(2*N^2+1:3*N^2,1)=-A11alpha*(ri1.*(R1*dispwc)+RR1*dispwc)...
    -A12alpha*(RR1*dispwc+ri1.*(R1*dispwc)+ri2.*(ThTh1*dispwc))...
    -A22alpha*(ri2.*(ThTh1*dispwc));
    %���߽�������ֵ��ԭ������Jacobian������DeltaT������У�3,1��*(N^2*1)�������
JacobiDeltaT(2*N^2+1:2*N^2+N,1)=0;
JacobiDeltaT(2*N^2+(N-1)*N+1:3*N^2,1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matrixZ=-JacobiY\JacobiDeltaT;%����������Z
% if dispwc(rowBI)>=0
%     JacobiAlphaPian=-(1+matrixZ.'*matrixZ)^(1/2);%���㸽�ӷ��̶�DeltaT��ƫ��alpha'
% else
%     JacobiAlphaPian=(1+matrixZ.'*matrixZ)^(1/2);%���㸽�ӷ��̶�DeltaT��ƫ��alpha'
% end
JacobiAlphaPian=-(1+matrixZ.'*matrixZ)^(1/2);%���㸽�ӷ��̶�DeltaT��ƫ��alpha'
% if matrixZlao'*matrixZ<0
%     JacobiAlphaPian=(1+matrixZ.'*matrixZ)^(1/2)
% end
JacobiXpian=matrixZ*JacobiAlphaPian;%���㸽�ӷ��̶�x��ƫ��x'