function zhi=CfuI2(e0,laGPL,h,rhoGpL,rhom);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FX=@(z)1-1.121*(1-(1-e0*cos(pi*z/h)).^(1/2.3));
% ccbi1=quadl(FX,-h/2,h/2);
% FY=@(z)(1-cos(pi*z/h)).*(1-1.121*(1-(1-e0*cos(pi*z/h)).^(1/2.3)));
% ccbi2=quadl(FY,-h/2,h/2);
sA=laGPL*rhom/(laGPL*rhom+rhoGpL-laGPL*rhoGpL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F=@(z)z.^2.*(rhoGpL.*sA+rhom.*(1-sA)).*(1-1.121.*(1-(1-e0.*cos(pi.*z/h)).^(1/2.3)));
zhi=quad(F,-h/2,h/2);
end