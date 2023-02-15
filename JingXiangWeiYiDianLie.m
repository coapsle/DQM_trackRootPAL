jnum=3;
for i=1:N
    csdw(i,1)=disp6((i-1)*N+jnum,3);
    csdphir(i,1)=disp6((i-1)*N+jnum,4);
    csdphitheta(i,1)=disp6((i-1)*N+jnum,5);
    xxx(i,1)=x((i-1)*N+jnum);
    yyy(i,1)=y((i-1)*N+jnum);
end
csdwmax=max(csdw);
csdwmin=min(csdw);
if abs(csdwmax)>abs(csdwmin)
    csdw=csdw/abs(csdwmax);
else
    csdw=csdw/abs(csdwmin);
end
csdphirmax=max(csdphir);
csdphirmin=min(csdphir);
if abs(csdphirmax)>abs(csdphirmin)
    csdphir=csdphir/abs(csdphirmax);
else
    csdphir=csdphir/abs(csdphirmin);
end
csdphithetamax=max(csdphitheta);
csdphithetamin=min(csdphitheta);
if abs(csdphithetamax)>abs(csdphithetamin)
    csdphitheta=csdphitheta/abs(csdphithetamax);
else
    csdphitheta=csdphitheta/abs(csdphithetamin);
end
rp=r';
figure
plot3(xxx,yyy,csdw)
figure
plot(r,csdw)
figure
plot3(xxx,yyy,csdphir)
figure
plot(r,csdphir)
figure
plot3(xxx,yyy,csdphitheta)
figure
plot(r,csdphitheta)


% for i=1:N
%     us(i,1)=dispuc((i-1)*N+jnum);
% end
% figure
% plot(r,us)
