function zhi=CfuB11(e0,laGPL,h,bbi,tbi,EGPL,nuGpL,rhoGpL,em,num,rhom,LGPL);
bGPL=LGPL/bbi;tGPL=LGPL/tbi;
xiL=2*LGPL/tGPL;xiB=2*bGPL/tGPL;
etaL=(EGPL-em)/(EGPL+xiL*em);
etaB=(EGPL-em)/(EGPL+xiB*em);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sA=laGPL*rhom/(laGPL*rhom+rhoGpL-laGPL*rhoGpL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FA11=@(z)z.*(((em*((3+3*xiL*etaL*sA)/(8-8*etaL*sA)+(5+5*xiB*etaB*sA)/(8-8*etaB*sA)))*(1-e0*cos(pi*z/h)))./(1-(0.342*(nuGpL*sA+num*(1-sA))*(1-1.121*(1-(1-e0*cos(pi*z/h)).^(1/2.3))).^2+(0.526*(nuGpL*sA+num*(1-sA))-0.221)*(1-1.121*(1-(1-e0*cos(pi*z/h)).^(1/2.3)))+0.132*(nuGpL*sA+num*(1-sA))+0.221).^2));
A11=quad(FA11,-h/2,h/2);
zhi=A11;
end