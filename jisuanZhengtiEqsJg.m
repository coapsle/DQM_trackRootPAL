%计算扩展后整体非线性方程组的结果G(x^k)
%将原方程组的结果复制到整体方程组结果的前5*N^2行
ZhengtiEqsJg(1:5*N^2,1)=YuanEqsJg;
%将附加方程的结果赋值到整体方程组的5*N^2+1行
ZhengtiEqsJg(5*N^2+1,1)=(dispc-dispc)'*JacobiXpian+(DeltaT-DeltaT)*JacobiAlphaPian-DeltaS;
