% Ra=2;
% Rb=1;
% N=15;
r=ones(1,N);
for k=1:1:N
    r(1,k)=Rb+(Ra-Rb)/2*(1-cos(pi*(k-1)/(N-1)));
end