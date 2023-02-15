for nii=1:30
lnum=s1(nii,2);
% DeltaT=real(1/s(lnum,lnum))
disp6=zeros(N^2,5);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果分列存储
for i=1:N^2
    disp6(i,1)=real(svector(i,lnum));
    disp6(i,2)=real(svector(i+N^2,lnum));
    disp6(i,3)=real(svector(i+2*N^2,lnum));
    disp6(i,4)=real(svector(i+3*N^2,lnum));
    disp6(i,5)=real(svector(i+4*N^2,lnum));    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果分列存储
disp61=zeros(N^2,5);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果取范数分列存储
for i=1:N^2
    disp61(i,1)=sign(real(svector(i,lnum)))*norm(svector(i,lnum));
    disp61(i,2)=sign(real(svector(i+N^2,lnum)))*norm(svector(i+N^2,lnum));
    disp61(i,3)=sign(real(svector(i+2*N^2,lnum)))*norm(svector(i+2*N^2,lnum));
    disp61(i,4)=sign(real(svector(i+3*N^2,lnum)))*norm(svector(i+3*N^2,lnum));
    disp61(i,5)=sign(real(svector(i+4*N^2,lnum)))*norm(svector(i+4*N^2,lnum));    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果取范数分列存储

% scatter(x,y,5,disp6(:,3))%散点图

% figure
% [X,Y,Z]=griddata(x,y,disp6(:,3),linspace(min(x),max(x))',linspace(min(y),max(y)),'v4');%插值
% pcolor(X,Y,Z);shading interp%伪彩色图

figure
plot3(x,y,disp6(:,3))
% figure,contourf(X,Y,Z) %等高线图
% figure,mesh(X,Y,Z)%三维曲面
% axis([-0.5 0.5 -0.5 0.5 -0.01 0.01]);

% figure
% [X,Y,Z]=griddata(x,y,disp61(:,3),linspace(min(x),max(x))',linspace(min(y),max(y)),'v4');%插值
% pcolor(X,Y,Z);shading interp%伪彩色图

figure
plot3(x,y,disp61(:,3))
% figure,contourf(X,Y,Z) %等高线图
% figure,mesh(X,Y,Z)%三维曲面
% axis([-0.5 0.5 -0.5 0.5 -0.01 0.01]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%为后屈曲路径上初始值做计算准备
% dispc=0.01*svector(:,lnum);
% disp5=zeros(N^2,5);
% for i=1:5
%     for j=1:N^2
%         disp5(j,i)=svector((i-1)*N^2+j,lnum);
%     end
% end
% qnw=zeros(N^2,1);
% for i=1:N^2
%     qnw(i,1)=qn1(2*N^2+i,1);
% end
% 
% dispuc=zeros(N^2,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果分列存储
% for i=1:N^2
%     dispuc(i,1)=disp5(i,1);
% end
% dispvc=zeros(N^2,1);
% for i=1:N^2
%     dispvc(i,1)=disp5(i,2);
% end
% dispwc=zeros(N^2,1);
% for i=1:N^2
%     dispwc(i,1)=disp5(i,3);
% end
% dispphirc=zeros(N^2,1);
% for i=1:N^2
%     dispphirc(i,1)=disp5(i,4);
% end
% dispphithc=zeros(N^2,1);
% for i=1:N^2
%     dispphithc(i,1)=disp5(i,5);
% end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将结果分列存储
% %建立ri矩阵
% II=eye(N^2);
% ri1=zeros(N^2,1);
% for i=1:1:N^2
%     ri1(i,1)=1/r(1,ceil(i/N));
% end
% ri2=zeros(N^2,1);
% for i=1:1:N^2
%     ri2(i,1)=1/(r(1,ceil(i/N))^2);
% end
% ri3=zeros(N^2,1);
% for i=1:1:N^2
%     ri3(i,1)=1/(r(1,ceil(i/N))^3);
% end
% ri4=zeros(N^2,1);
% for i=1:1:N^2
%     ri4(i,1)=1/(r(1,ceil(i/N))^4);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%