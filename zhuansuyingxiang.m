% jisuanCCBending;
% JacobiYuan;
% NewtonIterNLBending;
% JacobiYuan;
% Jacobizhengti;
% NewtonIterPostBuckling;
clear;
clc;
alphapie=zeros(1,1);
DeltaT=70;
n=0
N=13
q=-1030000;
Omega=1000;
kw=0;
kg=0;
parameter;
% Astiffness;
% Bstiffness;
Cstiffness;
jisuanCCBendingforPostbuckling;
NewtonIterBucklingInitial;
nnnn=1
matrixZlao=zeros(5*N^2,1);
while DeltaT<100
    jisuanXxingXxingpian;
    alphapie(nnnn,1)=JacobiAlphaPian;
    matrixZlao=matrixZ;
    NewtonIterPostBuckling;
%     if(nn>200)
%         disp('迭代步数太多，可能不收敛�?);
%         break;
%     end
    sss1(nnnn,1)=DeltaT;
    sss1(nnnn,2)=dispwc(rowBI);
%     sss1(nnnn,2)=max(dispwc);
    if mod(nnnn,1)==0
        sss2(1,nnnn/1)=DeltaT;
%         sss2(2,nnnn/10)=max(dispwc);
        sss2(2,nnnn/1)=dispwc(rowBI);
        sss2(3:size(dispZTc)+2,nnnn/1)=dispZTc;
    end
    if mod(nnnn,10)==0
        sss3(1,nnnn/10)=DeltaT;
%         sss3(2,nnnn/10)=max(dispwc);
        sss3(2,nnnn/10)=dispwc(rowBI);
        sss3(3:2+N^2,nnnn/10)=dispZTc(2*N^2+1:3*N^2);
    end
    nnnn=nnnn+1
    DeltaT
end
xlswrite('postbucklingT70PathPA2.xlsx',sss1);
xlswrite('postbucklingT70PathdispPA2.xlsx',sss2);
xlswrite('postbucklingT70PathdispwPA2.xlsx',sss3);