%��ֵ�Ⱥ��Ҳೣ������
Om1=zeros(5*N^2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N^2  %��ֵ������صĵȺ��Ҳೣ������
    Om1(i,1)=-2*I0*r(1,ceil(i/N))*Omega^2;
    Om1(3*N^2+i,1)=-2*I1*r(1,ceil(i/N))*Omega^2;
end