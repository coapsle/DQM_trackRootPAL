x1=th;
y1=csd;
ymax=max(y1);
ymin=min(y1);
if abs(ymax)>abs(ymin)
    y1=y1/abs(ymax);
else
    y1=y1/abs(ymin);
end
% y(1)=0;
% y(N)=0;
cftool
%%%%�Ͽ�
%%%%
param=coeffvalues(fittedmodel);%������ʽ�еĲ���
vpa(param)