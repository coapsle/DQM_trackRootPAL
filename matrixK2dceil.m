%¸³ÖµÏßÐÔÏµÊý¾ØÕóÖÐµÄ¸Õ¶ÈÕó
K2d=zeros(5*N^2);
% wspr=R1*dispwc;wsprr=RR1*dispwc;
% wspth=Th1*dispwc;wspthth=ThTh1*dispwc;
% wspthr=ThR*dispwc;
uspr=R1*dispuc;usprr=RR1*dispuc;
uspth=Th1*dispuc;uspthth=ThTh1*dispuc;
uspthr=ThR*dispuc;
% vspr=R1*dispvc;vsprr=RR1*dispvc;
% vspth=Th1*dispvc;vspthth=ThTh1*dispvc;
% vspthr=ThR*dispvc;
% phirspr=R1*dispphirc;phirsprr=RR1*dispphirc;
% phirspth=Th1*dispphirc;phirspthth=ThTh1*dispphirc;
% phirspthr=ThR*dispphirc;
% phithspr=R1*dispphithc;phithsprr=RR1*dispphithc;
% phithspth=Th1*dispphithc;phithspthth=ThTh1*dispphithc;
% phithspthr=ThR*dispphithc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N^2  %¸³Öµ£¨1,1£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i,j)=A11*(RR1(i,j)+1/r(1,ceil(i/N))*R1(i,j))-A22*1/r(1,ceil(i/N))^2*LL(i,j)+A44*1/r(1,ceil(i/N))^2*ThTh1(i,j)+I0*Omega^2*LL(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨1,2£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i,j+N^2)=A12*1/r(1,ceil(i/N))*ThR(i,j)-A22*1/r(1,ceil(i/N))^2*Th1(i,j)+A44*(1/r(1,ceil(i/N))*ThR(i,j)-1/r(1,ceil(i/N))^2*Th1(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨1,4£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i,j+3*N^2)=B11*(RR1(i,j)+1/r(1,ceil(i/N))*R1(i,j))-B22*1/r(1,ceil(i/N))^2*LL(i,j)+B44*1/r(1,ceil(i/N))^2*ThTh1(i,j)+I1*Omega^2*LL(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨1,5£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i,j+4*N^2)=B12*1/r(1,ceil(i/N))*ThR(i,j)-B22*1/r(1,ceil(i/N))^2*Th1(i,j)+B44*(1/r(1,ceil(i/N))*ThR(i,j)-1/r(1,ceil(i/N))^2*Th1(i,j));        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N^2  %¸³Öµ£¨2,1£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+N^2,j)=A12*1/r(1,ceil(i/N))*ThR(i,j)+A22*1/r(1,ceil(i/N))^2*Th1(i,j)+A44*(1/r(1,ceil(i/N))^2*Th1(i,j)+1/r(1,ceil(i/N))*ThR(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨2,2£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+N^2,j+N^2)=A22*1/r(1,ceil(i/N))^2*ThTh1(i,j)+A44*(RR1(i,j)-1/r(1,ceil(i/N))^2*LL(i,j)+1/r(1,ceil(i/N))*R1(i,j))+I0*Omega^2*LL(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨2,4£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+N^2,j+3*N^2)=B12*1/r(1,ceil(i/N))*ThR(i,j)+B22*1/r(1,ceil(i/N))^2*Th1(i,j)+B44*(1/r(1,ceil(i/N))^2*Th1(i,j)+1/r(1,ceil(i/N))*ThR(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨2,5£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+N^2,j+4*N^2)=B22*1/r(1,ceil(i/N))^2*ThTh1(i,j)+B44*(RR1(i,j)-1/r(1,ceil(i/N))^2*LL(i,j)+1/r(1,ceil(i/N))*R1(i,j))+I1*Omega^2*LL(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N^2  %¸³Öµ£¨3,3£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+2*N^2,j+2*N^2)=Ks*A55*(RR1(i,j)+1/r(1,ceil(i/N))*R1(i,j))+Ks*A66*1/r(1,ceil(i/N))^2*ThTh1(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨3,4£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+2*N^2,j+3*N^2)=Ks*A55*(R1(i,j)+1/r(1,ceil(i/N))*LL(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨3,5£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+2*N^2,j+4*N^2)=Ks*A66*1/r(1,ceil(i/N))*Th1(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N^2  %¸³Öµ£¨4,1£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+3*N^2,j)=B11*(RR1(i,j)+1/r(1,ceil(i/N))*R1(i,j))-B22*1/r(1,ceil(i/N))^2*LL(i,j)+B44*1/r(1,ceil(i/N))^2*ThTh1(i,j)+I1*Omega^2*LL(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨4,2£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+3*N^2,j+N^2)=B12*1/r(1,ceil(i/N))*ThR(i,j)-B22*1/r(1,ceil(i/N))^2*Th1(i,j)+B44*(1/r(1,ceil(i/N))*ThR(i,j)-1/r(1,ceil(i/N))^2*Th1(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨4,3£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+3*N^2,j+2*N^2)=-Ks*A55*R1(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨4,4£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+3*N^2,j+3*N^2)=-Ks*A55*LL(i,j)+D11*(RR1(i,j)+1/r(1,ceil(i/N))*R1(i,j))-D22*1/r(1,ceil(i/N))^2*LL(i,j)+D44*1/r(1,ceil(i/N))^2*ThTh1(i,j)+I2*Omega^2*LL(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨4,5£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+3*N^2,j+4*N^2)=D12*1/r(1,ceil(i/N))*ThR(i,j)-D22*1/r(1,ceil(i/N))^2*Th1(i,j)+D44*(1/r(1,ceil(i/N))*ThR(i,j)-1/r(1,ceil(i/N))^2*Th1(i,j));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N^2  %¸³Öµ£¨5,1£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+4*N^2,j)=B12*1/r(1,ceil(i/N))*ThR(i,j)+B22*1/r(1,ceil(i/N))^2*Th1(i,j)+B44*(1/r(1,ceil(i/N))^2*Th1(i,j)+1/r(1,ceil(i/N))*ThR(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨5,2£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+4*N^2,j+N^2)=B22*1/r(1,ceil(i/N))^2*ThTh1(i,j)+B44*(RR1(i,j)-1/r(1,ceil(i/N))^2*LL(i,j)+1/r(1,ceil(i/N))*R1(i,j))+I1*Omega^2*LL(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨5,3£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+4*N^2,j+2*N^2)=-Ks*A66*1/r(1,ceil(i/N))*Th1(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨5,4£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+4*N^2,j+3*N^2)=D12*1/r(1,ceil(i/N))*ThR(i,j)+D22*1/r(1,ceil(i/N))^2*Th1(i,j)+D44*(1/r(1,ceil(i/N))^2*Th1(i,j)+1/r(1,ceil(i/N))*ThR(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨5,5£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+4*N^2,j+4*N^2)=-Ks*A66*LL(i,j)+D22*1/r(1,ceil(i/N))^2*ThTh1(i,j)+D44*(RR1(i,j)-1/r(1,ceil(i/N))^2*LL(i,j)+1/r(1,ceil(i/N))*R1(i,j))+I2*Omega^2*LL(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%¸³Öµº¬³õÊ¼¹¹ÐÍus,vs,ws,phirs,phithetasÓ°ÏìµÄ¾ØÕó¿é
% for i=1:N^2  %¸³Öµ£¨1,3£©*(N^2*N^2)¾ØÕó¿é
%     for j=1:N^2
%         K2d(i,j+2*N^2)=A11*(wspr(i)*RR1(i,j)+wsprr(i)*R1(i,j)+ri1(i)*wspr(i)*R1(i,j))...
%             +A12*(-ri3(i)*wspth(i)*Th1(i,j)+ri2(i)*(wspth(i)*ThR(i,j)+wspthr(i)*Th1(i,j))...
%             -ri1(i)*wspr(i)*R1(i,j))...
%             -A22*ri3(i)*wspth(i)*Th1(i,j)...
%             +A44*(ri2(i)*(wspthr(i)*Th1(i,j)+wspth(i)*ThR(i,j))...
%             +ri2(i)*(wspr(i)*ThTh1(i,j)+wspthth(i)*R1(i,j)));
%     end
% end
% for i=1:N^2  %¸³Öµ£¨2,3£©*(N^2*N^2)¾ØÕó¿é
%     for j=1:N^2
%         K2d(i+N^2,j+2*N^2)=A12*ri1(i)*(wspr(i)*ThR(i,j)+wspthr(i)*R1(i,j))...
%             +A22*ri3(i)*(wspth(i)*ThTh1(i,j)+wspthth(i)*Th1(i,j))...
%             +A44*(ri2(i)*(wspr(i)*Th1(i,j)+wspth(i)*R1(i,j))...
%             +ri1(i)*(wsprr(i)*Th1(i,j)+wspth(i)*RR1(i,j))...
%             +ri1(i)*(wspr(i)*ThR(i,j)+wspthr(i)*R1(i,j)));
%     end
% end
% for i=1:N^2  %¸³Öµ£¨4,3£©*(N^2*N^2)¾ØÕó¿é
%     for j=1:N^2
%         K2d(i+3*N^2,j+2*N^2)=K2d(i+3*N^2,j+2*N^2)+B11*(ri1(i)*wspr(i)*R1(i,j)...
%             +(wspr(i)*RR1(i,j)+wsprr(i)*R1(i,j)))...
%             +B12*(-ri3(i)*wspth(i)*Th1(i,j)-ri1(i)*wspr(i)*R1(i,j)...
%             +ri2(i)*(wspth(i)*ThR(i,j)+wspthr(i)*Th1(i,j)))...
%             +B22*ri3(i)*wspth(i)*Th1(i,j)...
%             +B44*(ri2(i)*(wspthr(i)*Th1(i,j)+wspth(i)*ThR(i,j))...
%             +ri2(i)*(wspr(i)*ThTh1(i,j)+wspthth(i)*R1(i,j)));
%     end
% end
% for i=1:N^2  %¸³Öµ£¨5,3£©*(N^2*N^2)¾ØÕó¿é
%     for j=1:N^2
%         K2d(i+4*N^2,j+2*N^2)=K2d(i+4*N^2,j+2*N^2)+B12*ri1(i)*(wspr(i)*ThR(i,j)+wspthr(i)*R1(i,j))...
%             +B22*ri3(i)*(wspth(i)*ThTh1(i,j)+wspthth(i)*Th1(i,j))...
%             +B44*(ri2(i)*(wspr(i)*Th1(i,j)+wspth(i)*R1(i,j))...
%             +ri1(i)*(wsprr(i)*Th1(i,j)+wspth(i)*RR1(i,j))...
%             +ri1(i)*(wspr(i)*ThR(i,j)+wspthr(i)*R1(i,j)));
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:N^2  %¸³Öµ£¨3,1£©*(N^2*N^2)¾ØÕó¿é
%     for j=1:N^2
%         K2d(i+2*N^2,j)=A11*(ri1(i)*wspr(i)*R1(i,j)+wspr(i)*RR1(i,j)+wsprr(i)*R1(i,j))...
%             +A12*(ri1(i)*wspr(i)*R1(i,j)+ri1(i)*wsprr(i)*LL(i,j)...
%             +ri2(i)*wspth(i)*ThR(i,j)+ri2(i)*wspthth(i)*R1(i,j))...
%             +A22*(ri3(i)*wspth(i)*Th1(i,j)+ri3(i)*wspthth(i)*LL(i,j))...
%             +A44*(-ri3(i)*wspth(i)*Th1(i,j)+ri2(i)*wspth(i)*ThR(i,j)...
%             +ri2(i)*wspr(i)*ThTh1(i,j)+2*ri2(i)*wspthr(i)*Th1(i,j));
%     end
% end
% for i=1:N^2  %¸³Öµ£¨3,2£©*(N^2*N^2)¾ØÕó¿é
%     for j=1:N^2
%         K2d(i+2*N^2,j+N^2)=A12*(ri1(i)*wspr(i)*ThR(i,j)+ri1(i)*wsprr(i)*Th1(i,j))...
%             +A22*(ri3(i)*wspth(i)*ThTh1(i,j)+ri3(i)*wspthth(i)*Th1(i,j))...
%             +A44*(ri1(i)*wspth(i)*RR1(i,j)+ri3(i)*wspth(i)*LL(i,j)...
%             -ri2(i)*wspth(i)*R1(i,j)+ri1(i)*wspr(i)*ThR(i,j)...
%             -ri2(i)*wspr(i)*Th1(i,j)+2*ri1(i)*wspthr(i)*R1(i,j)-2*ri2(i)*wspthr(i)*LL(i,j));
%     end
% end
% for i=1:N^2  %¸³Öµ£¨3,4£©*(N^2*N^2)¾ØÕó¿é
%     for j=1:N^2
%         K2d(i+2*N^2,j+3*N^2)=K2d(i+2*N^2,j+3*N^2)+B11*(ri1(i)*wspr(i)*R1(i,j)...
%             +wspr(i)*RR1(i,j)+wsprr(i)*R1(i,j))...
%             +B12*(ri1(i)*wspr(i)*R1(i,j)+ri1(i)*wsprr(i)*LL(i,j)...
%             +ri2(i)*wspth(i)*ThR(i,j)+ri2(i)*wspthth(i)*R1(i,j))...
%             +B22*(ri3(i)*wspth(i)*Th1(i,j)+ri3(i)*wspthth(i)*LL(i,j))...
%             +B44*(-ri3(i)*wspth(i)*Th1(i,j)+ri2(i)*wspth(i)*ThR(i,j)...
%             +ri2(i)*wspr(i)*ThTh1(i,j)+2*ri2(i)*wspthr(i)*Th1(i,j));
%     end
% end
% for i=1:N^2  %¸³Öµ£¨3,5£©*(N^2*N^2)¾ØÕó¿é
%     for j=1:N^2
%         K2d(i+2*N^2,j+4*N^2)=K2d(i+2*N^2,j+4*N^2)+B12*(ri1(i)*wspr(i)*ThR(i,j)+ri1(i)*wsprr(i)*Th1(i,j))...
%             +B22*(ri3(i)*wspth(i)*ThTh1(i,j)+ri3(i)*wspthth(i)*Th1(i,j))...
%             +B44*(ri1(i)*wspth(i)*RR1(i,j)+ri3(i)*wspth(i)*LL(i,j)...
%             -ri2(i)*wspth(i)*R1(i,j)+ri1(i)*wspr(i)*ThR(i,j)...
%             -ri2(i)*wspr(i)*Th1(i,j)+2*ri1(i)*wspthr(i)*R1(i,j)-2*ri2(i)*wspthr(i)*LL(i,j));
%     end
% end
for i=1:N^2  %¸³Öµ£¨3,3£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K2d(i+2*N^2,j+2*N^2)=K2d(i+2*N^2,j+2*N^2)...
            +A11*(ri1(i)*uspr(i)*R1(i,j)+usprr(i)*R1(i,j)+uspr(i)*RR1(i,j))...
            +A12*(ri1(i)*uspr(i)*R1(i,j)+ri2(i)*uspthr(i)*Th1(i,j)+ri2(i)*uspr(i)*ThTh1(i,j)...
            +ri1(i)*dispuc(i)*RR1(i,j))...
            +A22*(ri3(i)*uspth(i)*Th1(i,j)+ri3(i)*dispuc(i)*ThTh1(i,j))...
            +A44*(-ri3(i)*uspth(i)*Th1(i,j)+ri2(i)*uspthr(i)*Th1(i,j)+ri2(i)*uspthth(i)*R1(i,j)...
            +2*ri2(i)*uspth(i)*ThR(i,j));
%             +A11*(...
%             +3/2*ri1(i)*wspr(i)*wspr(i)*R1(i,j)...
%             +3/2*(wspr(i)*wspr(i)*RR1(i,j)+2*wspr(i)*wsprr(i)*R1(i,j)))...
%             +A12*(ri1(i)*vspthr(i)*R1(i,j)...
%             -1/2*ri3(i)*(wspth(i)*wspth(i)*R1(i,j)+2*wspth(i)*wspr(i)*Th1(i,j))...
%             +2*ri2(i)*(wspth(i)*wspr(i)*ThR(i,j)+wspr(i)*wspthr(i)*Th1(i,j)+wspth(i)*wspthr(i)*R1(i,j))...
%             +ri1(i)*vspth(i)*RR1(i,j)...
%             +1/2*ri2(i)*(wspth(i)*wspth(i)*RR1(i,j)+2*wspth(i)*wsprr(i)*Th1(i,j))...
%             ...
%             +1/2*ri2(i)*(wspr(i)*wspr(i)*ThTh1(i,j)+2*wspr(i)*wspthth(i)*R1(i,j)))...
%             +A22*(ri3(i)*vspthth(i)*Th1(i,j)+ri3(i)*vspth(i)*ThTh1(i,j)...
%             +3/2*ri4(i)*(wspth(i)*wspth(i)*ThTh1(i,j)+2*wspth(i)*wspthth(i)*Th1(i,j)))...
%             +A44*(ri1(i)*vsprr(i)*Th1(i,j)...
%             +ri3(i)*dispvc(i)*Th1(i,j)-ri2(i)*vspr(i)*Th1(i,j)...
%             -ri3(i)*(wspth(i)*wspth(i)*R1(i,j)+2*wspth(i)*wspr(i)*Th1(i,j))...
%             +ri2(i)*(wspth(i)*wspth(i)*RR1(i,j)+2*wspth(i)*wsprr(i)*Th1(i,j))...
%             +ri1(i)*vspthr(i)*R1(i,j)-ri2(i)*vspth(i)*R1(i,j)...
%             +ri2(i)*(wspr(i)*wspr(i)*ThTh1(i,j)+2*wspr(i)*wspthth(i)*R1(i,j))...
%             +2*ri1(i)*vspr(i)*ThR(i,j)-2*ri2(i)*dispvc(i)*ThR(i,j)...
%             +4*ri2(i)*(wspth(i)*wspr(i)*ThR(i,j)+wspth(i)*wspthr(i)*R1(i,j)+wspr(i)*wspthr(i)*Th1(i,j)))...
%             +B11*(ri1(i)*phirspr(i)*R1(i,j)+phirsprr(i)*R1(i,j)+phirspr(i)*RR1(i,j))...
%             +B12*(ri1(i)*phirspr(i)*R1(i,j)+ri1(i)*phithspthr(i)*R1(i,j)+ri1(i)*dispphirc(i)*RR1(i,j)...
%             +ri1(i)*phithspth(i)*RR1(i,j)+ri2(i)*phirspthr(i)*Th1(i,j)+ri2(i)*phirspr(i)*ThTh1(i,j))...
%             +B22*(ri3(i)*phirspth(i)*Th1(i,j)+ri3(i)*phithspthth(i)*Th1(i,j)...
%             +ri3(i)*dispphirc(i)*ThTh1(i,j)+ri3(i)*phithspth(i)*ThTh1(i,j))...
%             +B44*(-ri3(i)*phirspth(i)*Th1(i,j)+ri2(i)*phirspthr(i)*Th1(i,j)...
%             +ri1(i)*phithsprr(i)*Th1(i,j)+ri3(i)*dispphithc(i)*Th1(i,j)...
%             -ri2(i)*phithspr(i)*Th1(i,j)+ri2(i)*phirspthth(i)*R1(i,j)+ri1(i)*phithspthr(i)*R1(i,j)...
%             -ri2(i)*phithspth(i)*R1(i,j)+2*ri2(i)*phirspth(i)*ThR(i,j)...
%             +2*ri1(i)*phithspr(i)*ThR(i,j)-2*ri2(i)*dispphithc(i)*ThR(i,j));
    end
end
%%%%¸³Öµº¬³õÊ¼¹¹ÐÍus,vs,ws,phirs,phithetasÓ°ÏìµÄ¾ØÕó¿é