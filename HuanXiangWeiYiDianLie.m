csd=disp6(ceil(N/4)*N+1:ceil(N/4)*N+N,4);
xx=x(1:N);
yy=y(1:N);
figure
plot3(xx,yy,csd)
figure
plot(th,csd)
