% Ra=2;
% Rb=1;
% N=15;
th=ones(1,N);
for k=1:1:N
%     th(1,k)=pi*(1-cos(0.99*pi*(k-1)/(N-1))); %非对称显示，不适用
    th(1,k)=pi*(1-cos(pi*(k-1)/(N-1)));
%     th(1,k)=2*pi*(k-1)/(N); %非对称显示，不适用
%     th(1,k)=2*pi*(k-1)/(N-1);
end