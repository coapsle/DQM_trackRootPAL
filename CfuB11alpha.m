function zhi=CfuB11alpha(e0,laGPL,h,bbi,tbi,EGPL,nuGpL,rhoGpL,aaGpL,em,num,rhom,aam,LGPL);
%syms z
bGPL=LGPL/bbi;tGPL=LGPL/tbi;
xiL=2*LGPL/tGPL;xiB=2*bGPL/tGPL;
etaL=(EGPL-em)/(EGPL+xiL*em);
etaB=(EGPL-em)/(EGPL+xiB*em);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sA=laGPL*rhom/(laGPL*rhom+rhoGpL-laGPL*rhoGpL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aa=aaGpL*sA+aam*(1-sA);
% Q11=@(z)((em*((3+3*xiL*etaL*sA)/(8-8*etaL*sA)+(5+5*xiB*etaB*sA)/(8-8*etaB*sA)))*(1-e0*cos(pi*z/h)))./(1-(0.342*(nuGpL*sA+num*(1-sA))*(1-1.121*(1-(1-e0*cos(pi*z/h)).^(1/2.3))).^2+(0.526*(nuGpL*sA+num*(1-sA))-0.221)*(1-1.121*(1-(1-e0*cos(pi*z/h)).^(1/2.3)))+0.132*(nuGpL*sA+num*(1-sA))+0.221).^2);
% Q12=(em.*((5+5.*etaB.*sA.*xiB)./(8-8.*etaB.*sA)+(3+3.*etaL.*sA.*xiL)./(8-8.*etaL.*sA)))*(nuGpL*sA+num*(1-sA))/(1-(nuGpL*sA+num*(1-sA))^2);
F=@(z)z.*(z.*0+((em*((3+3*xiL*etaL*sA)/(8-8*etaL*sA)+(5+5*xiB*etaB*sA)/(8-8*etaB*sA)))*(1-e0*cos(pi*z/h)))./(1-(0.342*(nuGpL*sA+num*(1-sA))*(1-1.121*(1-(1-e0*cos(pi*z/h)).^(1/2.3))).^2+(0.526*(nuGpL*sA+num*(1-sA))-0.221)*(1-1.121*(1-(1-e0*cos(pi*z/h)).^(1/2.3)))+0.132*(nuGpL*sA+num*(1-sA))+0.221).^2).*aa);
% zhi=double(int(F,-h/2,h/2));
zhi=quad(F,-h/2,h/2);
end