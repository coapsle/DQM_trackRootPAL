%��ֵJacobian����
JacobiMNL=zeros(5*N^2);
%����ri����
II=eye(N^2);
ri1=zeros(N^2,1);
for i=1:1:N^2
    ri1(i,1)=1/r(1,ceil(i/N));
end
ri2=zeros(N^2,1);
for i=1:1:N^2
    ri2(i,1)=1/(r(1,ceil(i/N))^2);
end
ri3=zeros(N^2,1);
for i=1:1:N^2
    ri3(i,1)=1/(r(1,ceil(i/N))^3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ֵ��1,3��*(N^2*N^2)�����
JacobiMNL(1:N^2,2*N^2+1:3*N^2)=A11*(SJTP(RR1*dispwc,R1)+SJTP(R1*dispwc,RR1)...
    +1/2*SJTP(ri1,(2*SJTP(R1*dispwc,R1))))...
    +A12*(-1/2*SJTP(ri3,2*SJTP(Th1*dispwc,Th1))...
    +SJTP(ri2,(SJTP(ThR*dispwc,Th1)+SJTP(Th1*dispwc,ThR)))...
    -1/2*SJTP(ri1,2*SJTP(R1*dispwc,R1)))...
    -A22*1/2*SJTP(ri3,2*SJTP(Th1*dispwc,Th1))...
    +A44*(SJTP(ri2,SJTP(ThR*dispwc,Th1)+SJTP(Th1*dispwc,ThR))...
    +SJTP(ri2,SJTP(ThTh1*dispwc,R1)+SJTP(R1*dispwc,ThTh1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ֵ��2,3��*(N^2*N^2)�����
JacobiMNL(N^2+1:2*N^2,2*N^2+1:3*N^2)=A12*SJTP(ri1,SJTP(ThR*dispwc,R1)+SJTP(R1*dispwc,ThR))...
    +A22*SJTP(ri3,SJTP(ThTh1*dispwc,Th1)+SJTP(Th1*dispwc,ThTh1))...
    +A44*(SJTP(ri2,SJTP(Th1*dispwc,R1)+SJTP(R1*dispwc,Th1))...
    +SJTP(ri1,SJTP(RR1*dispwc,Th1)+SJTP(Th1*dispwc,RR1))...
    +SJTP(ri1,SJTP(ThR*dispwc,R1)+SJTP(R1*dispwc,ThR)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ֵ��3,1��*(N^2*N^2)�����
JacobiMNL(2*N^2+1:3*N^2,1:N^2)=A11*(SJTP(ri1,SJTP(R1*dispwc,R1))...
    +SJTP(R1*dispwc,RR1)+SJTP(RR1*dispwc,R1))...
    +A12*(SJTP(ri1,SJTP(R1*dispwc,R1))+SJTP(ri1,SJTP(RR1*dispwc,II))...
    +SJTP(ri2,SJTP(R1*dispwc,ThR))+SJTP(ri2,SJTP(ThTh1*dispwc,R1)))...
    +A22*(SJTP(ri3,SJTP(Th1*dispwc,Th1))+SJTP(ri3,SJTP(ThTh1*dispwc,II)))...
    +A44*(-SJTP(ri3,SJTP(Th1*dispwc,Th1))+SJTP(ri2,SJTP(Th1*dispwc,ThR))...
    +SJTP(ri2,SJTP(R1*dispwc,ThTh1))+2*SJTP(ri2,SJTP(ThR*dispwc,Th1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ֵ��3,2��*(N^2*N^2)�����
JacobiMNL(2*N^2+1:3*N^2,N^2+1:2*N^2)=A12*(SJTP(ri1,SJTP(R1*dispwc,ThR))...
    +SJTP(ri1,SJTP(RR1*dispwc,Th1)))...
    +A22*(SJTP(ri3,SJTP(Th1*dispwc,ThTh1))+SJTP(ri3,SJTP(ThTh1*dispwc,Th1)))...
    +A44*(SJTP(ri1,SJTP(Th1*dispwc,RR1))+SJTP(ri3,SJTP(Th1*dispwc,II))...
    -SJTP(ri2,SJTP(Th1*dispwc,R1))+SJTP(ri1,SJTP(R1*dispwc,ThR))...
    -SJTP(ri2,SJTP(R1*dispwc,Th1))+2*SJTP(ri1,SJTP(ThR*dispwc,R1))...
    -2*SJTP(ri2,SJTP(ThR*dispwc,II)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ֵ��3,3��*(N^2*N^2)�����
JacobiMNL(2*N^2+1:3*N^2,2*N^2+1:3*N^2)=A11*(SJTP(ri1,SJTP(R1*dispuc,R1))...
    +SJTP(RR1*dispuc,R1)+SJTP(R1*dispuc,RR1))...
    +A12*(SJTP(ri1,SJTP(R1*dispuc,R1))+SJTP(ri1,SJTP(II*dispuc,RR1))...
    +SJTP(ri2,SJTP(ThR*dispuc,Th1))+SJTP(R1*dispuc,ThTh1))...
    +A22*(SJTP(ri3,SJTP(Th1*dispuc,Th1))+SJTP(ri3,SJTP(II*dispuc,ThTh1)))...
    +A44*(-SJTP(ri3,SJTP(Th1*dispuc,Th1))+SJTP(ri2,SJTP(ThR*dispuc,Th1))...
    +SJTP(ri2,SJTP(ThTh1*dispuc,R1))+2*SJTP(ri2,SJTP(Th1*dispuc,ThR)))...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    +A12*(SJTP(ri1,SJTP(ThR*dispvc,R1))+SJTP(ri1,SJTP(Th1*dispvc,RR1)))...
    +A12*(SJTP(ri3,SJTP(ThTh1*dispvc,Th1))+SJTP(ri3,SJTP(Th1*dispvc,ThTh1)))...
    +A44*(SJTP(ri1,SJTP(RR1*dispvc,Th1))+SJTP(ri3,SJTP(II*dispvc,Th1))...
    -SJTP(ri2,SJTP(R1*dispvc,Th1))+SJTP(ri1,SJTP(ThR*dispvc,R1))...
    -SJTP(ri2,SJTP(Th1*dispvc,R1))+2*SJTP(ri1,SJTP(R1*dispvc,ThR))...
    -2*SJTP(ri2,SJTP(II*dispvc,ThR)))...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    +B11*(SJTP(ri1,SJTP(R1*dispphirc,R1))+SJTP(RR1*dispphirc,R1)...
    +SJTP(R1*dispphirc,RR1))...
    +B12*(SJTP(ri1,SJTP(R1*dispphirc,R1))+SJTP(ri1,SJTP(II*dispphirc,RR1))...
    +SJTP(ri2,SJTP(ThR*dispphirc,Th1))+SJTP(ri2,SJTP(R1*dispphirc,ThTh1)))...
    +B22*(SJTP(ri3,SJTP(Th1*dispphirc,Th1))+SJTP(ri3,SJTP(II*dispphirc,ThTh1)))...
    +B44*(-SJTP(ri3,SJTP(Th1*dispphirc,Th1))+SJTP(ri2,SJTP(ThR*dispphirc,Th1))...
    +SJTP(ri2,SJTP(ThTh1*dispphirc,R1))+2*SJTP(ri2,SJTP(Th1*dispphirc,ThR)))...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    +B12*(SJTP(ri1,SJTP(ThR*dispphithc,R1))+SJTP(ri1,SJTP(Th1*dispphithc,RR1)))...
    +B22*(SJTP(ri3,SJTP(ThTh1*dispphithc,Th1))+SJTP(ri3,SJTP(Th1*dispphithc,ThTh1)))...
    +B44*(SJTP(ri1,SJTP(RR1*dispphithc,Th1))+SJTP(ri3,SJTP(II*dispphithc,Th1))...
    -SJTP(ri1,SJTP(R1*dispphithc,Th1))+SJTP(ri1,SJTP(ThR*dispphithc,R1))...
    -SJTP(ri2,SJTP(Th1*dispphithc,R1))+2*SJTP(ri1,SJTP(Th1*dispphithc,ThR))...
    -2*SJTP(ri2,SJTP(II*dispphithc,ThR)))...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    +A11*(3/2*SJTP(ri1,SJTP((R1*dispwc).*(R1*dispwc),R1))...
    +3/2*(2*SJTP((R1*dispwc).*(RR1*dispwc),R1)+SJTP((R1*dispwc).*(R1*dispwc),RR1)))...
    +A12*(-1/2*SJTP(ri3,2*SJTP((Th1*dispwc).*(R1*dispwc),Th1)+SJTP((Th1*dispwc).*(Th1*dispwc),R1))...
    +2*SJTP(ri2,SJTP((Th1*dispwc).*(ThR*dispwc),R1)+SJTP((R1*dispwc).*(ThR*dispwc),Th1)...
    +SJTP((R1*dispwc).*(Th1*dispwc),ThR))...
    +1/2*SJTP(ri2,2*SJTP((Th1*dispwc).*(RR1*dispwc),Th1)+SJTP((Th1*dispwc).*(Th1*dispwc),RR1))...
    +1/2*SJTP(ri2,2*SJTP((R1*dispwc).*(ThTh1*dispwc),R1)+SJTP((R1*dispwc).*(R1*dispwc),ThTh1)))...
    +A22*(3/2*SJTP(ri4,2*SJTP((Th1*dispwc).*(ThTh1*dispwc),Th1)+SJTP((Th1*dispwc).*(Th1*dispwc),ThTh1)))...
    +A44*(-SJTP(ri3,2*SJTP((R1*dispwc).*(Th1*dispwc),Th1)+SJTP((Th1*dispwc).*(Th1*dispwc),R1))...
    +SJTP(ri2,2*SJTP((RR1*dispwc).*(Th1*dispwc),Th1)+SJTP((Th1*dispwc).*(Th1*dispwc),RR1))...
    +SJTP(ri2,2*SJTP((R1*dispwc).*(ThTh1*dispwc),R1)+SJTP((R1*dispwc).*(R1*dispwc),ThTh1))...
    +4*SJTP(ri4,SJTP((Th1*dispwc)*(ThR*dispwc),R1)+SJTP((R1*dispwc).*(ThR*dispwc),Th1)...
    +SJTP((R1*dispwc).*(Th1*dispwc),ThR)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ֵ��3,4��*(N^2*N^2)�����
JacobiMNL(2*N^2+1:3*N^2,3*N^2+1:4*N^2)=B11*(SJTP(ri1,SJTP(R1*dispwc,R1))+SJTP(R1*dispwc,RR1)...
    +SJTP(RR1*dispwc,R1))...
    +B12*(SJTP(ri1,SJTP(R1*dispwc,R1))+SJTP(ri1,SJTP(RR1*dispwc,II))...
    +SJTP(ri2,SJTP(Th1*dispwc,ThR))+SJTP(ri2,SJTP(ThTh1*dispwc,R1)))...
    +B22*(SJTP(ri3,SJTP(Th1*dispwc,Th1))+SJTP(ri3,SJTP(ThTh1*dispwc,II)))...
    +B44*(-SJTP(ri3,SJTP(Th1*dispwc,Th1))+SJTP(ri2,SJTP(Th1*dispwc,ThR))...
    -SJTP(ri2,SJTP(Th1*dispwc,R1))+SJTP(ri2,SJTP(R1*dispwc,ThTh1))...
    +2*SJTP(ri2,SJTP(ThR*dispwc,Th1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ֵ��3,5��*(N^2*N^2)�����
JacobiMNL(2*N^2+1:3*N^2,4*N^2+1:5*N^2)=B12*(SJTP(ri1,SJTP(R1*dispwc,ThR))...
    +SJTP(RR1*dispwc,Th1))...
    +B22*(SJTP(ri3,SJTP(Th1*dispwc,ThTh1))+SJTP(ri3,SJTP(ThTh1*dispwc,Th1)))...
    +B44*(SJTP(ri1,SJTP(Th1*dispwc,RR1))+SJTP(ri3,SJTP(Th1*dispwc,II))...
    -SJTP(ri2,SJTP(Th1*dispwc,R1))+SJTP(ri1,SJTP(R1*dispwc,ThR))...
    -SJTP(ri2,SJTP(R1*dispwc,Th1))+2*SJTP(ri1,SJTP(ThR*dispwc,R1))...
    -2*SJTP(ri2,SJTP(ThR*dispwc,II)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ֵ��4,3��*(N^2*N^2)�����
JacobiMNL(3*N^2+1:4*N^2,2*N^2+1:3*N^2)=B11*(SJTP(RR1*dispwc,R1)+SJTP(R1*dispwc,RR1)...
    +1/2*SJTP(ri1,2*SJTP(R1*dispwc,R1)))...
    +B12*(-1/2*SJTP(ri3,2*SJTP(Th1*dispwc,Th1))...
    +SJTP(ri2,SJTP(ThR*dispwc,Th1)+SJTP(Th1*dispwc,ThR))...
    -1/2*SJTP(ri1,2*SJTP(R1*dispwc,R1)))...
    -B22*(1/2*SJTP(ri3,2*SJTP(Th1*dispwc,Th1)))...
    +B44*(SJTP(ri2,SJTP(ThR*dispwc,Th1)+SJTP(Th1*dispwc,ThR))...
    +SJTP(ri2,SJTP(ThTh1*dispwc,R1)+SJTP(R1*dispwc,ThTh1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ֵ��5,3��*(N^2*N^2)�����
JacobiMNL(4*N^2+1:5*N^2,2*N^2+1:3*N^2)=B12*SJTP(ri1,2*SJTP((R1*dispwc).*(ThR*dispwc),R1)...
    +SJTP((R1*dispwc).*(R1*dispwc),ThR))...
    +B22*SJTP(ri3,SJTP(ThTh1*dispwc,Th1)+SJTP(Th1*dispwc,ThTh1))...
    +B44*(SJTP(ri2,SJTP(R1*dispwc,Th1)+SJTP(Th1*dispwc,R1))...
    +SJTP(ri1,SJTP(RR1*dispwc,Th1)+SJTP(Th1*dispwc,RR1))...
    +SJTP(ri1,SJTP(ThR*dispwc,R1)+SJTP(R1*dispwc,ThR)));