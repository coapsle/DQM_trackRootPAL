%¸³ÖµÏßÐÔÏµÊý¾ØÕóÖÐµÄ¸Õ¶ÈÕó
K1=zeros(5*N^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N^2  %¸³Öµ£¨1,1£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i,j)=A11*(RR1(i,j)+1/r(1,ceil(i/N))*R1(i,j))-A22*1/r(1,ceil(i/N))^2*LL(i,j)+A44*1/r(1,ceil(i/N))^2*ThTh1(i,j)+I0*Omega^2*LL(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨1,2£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i,j+N^2)=A12*1/r(1,ceil(i/N))*ThR(i,j)-A22*1/r(1,ceil(i/N))^2*Th1(i,j)+A44*(1/r(1,ceil(i/N))*ThR(i,j)-1/r(1,ceil(i/N))^2*Th1(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨1,4£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i,j+3*N^2)=B11*(RR1(i,j)+1/r(1,ceil(i/N))*R1(i,j))-B22*1/r(1,ceil(i/N))^2*LL(i,j)+B44*1/r(1,ceil(i/N))^2*ThTh1(i,j)+I1*Omega^2*LL(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨1,5£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i,j+4*N^2)=B12*1/r(1,ceil(i/N))*ThR(i,j)-B22*1/r(1,ceil(i/N))^2*Th1(i,j)+B44*(1/r(1,ceil(i/N))*ThR(i,j)-1/r(1,ceil(i/N))^2*Th1(i,j));        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N^2  %¸³Öµ£¨2,1£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+N^2,j)=A12*1/r(1,ceil(i/N))*ThR(i,j)+A22*1/r(1,ceil(i/N))^2*Th1(i,j)+A44*(1/r(1,ceil(i/N))^2*Th1(i,j)+1/r(1,ceil(i/N))*ThR(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨2,2£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+N^2,j+N^2)=A22*1/r(1,ceil(i/N))^2*ThTh1(i,j)+A44*(RR1(i,j)-1/r(1,ceil(i/N))^2*LL(i,j)+1/r(1,ceil(i/N))*R1(i,j))+I0*Omega^2*LL(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨2,4£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+N^2,j+3*N^2)=B12*1/r(1,ceil(i/N))*ThR(i,j)+B22*1/r(1,ceil(i/N))^2*Th1(i,j)+B44*(1/r(1,ceil(i/N))^2*Th1(i,j)+1/r(1,ceil(i/N))*ThR(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨2,5£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+N^2,j+4*N^2)=B22*1/r(1,ceil(i/N))^2*ThTh1(i,j)+B44*(RR1(i,j)-1/r(1,ceil(i/N))^2*LL(i,j)+1/r(1,ceil(i/N))*R1(i,j))+I1*Omega^2*LL(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N^2  %¸³Öµ£¨3,3£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+2*N^2,j+2*N^2)=Ks*A55*(RR1(i,j)+1/r(1,ceil(i/N))*R1(i,j))+Ks*A66*1/r(1,ceil(i/N))^2*ThTh1(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨3,4£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+2*N^2,j+3*N^2)=Ks*A55*(R1(i,j)+1/r(1,ceil(i/N))*LL(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨3,5£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+2*N^2,j+4*N^2)=Ks*A66*1/r(1,ceil(i/N))*Th1(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N^2  %¸³Öµ£¨4,1£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+3*N^2,j)=B11*(RR1(i,j)+1/r(1,ceil(i/N))*R1(i,j))-B22*1/r(1,ceil(i/N))^2*LL(i,j)+B44*1/r(1,ceil(i/N))^2*ThTh1(i,j)+I1*Omega^2*LL(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨4,2£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+3*N^2,j+N^2)=B12*1/r(1,ceil(i/N))*ThR(i,j)-B22*1/r(1,ceil(i/N))^2*Th1(i,j)+B44*(1/r(1,ceil(i/N))*ThR(i,j)-1/r(1,ceil(i/N))^2*Th1(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨4,3£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+3*N^2,j+2*N^2)=-Ks*A55*R1(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨4,4£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+3*N^2,j+3*N^2)=-Ks*A55*LL(i,j)+D11*(RR1(i,j)+1/r(1,ceil(i/N))*R1(i,j))-D22*1/r(1,ceil(i/N))^2*LL(i,j)+D44*1/r(1,ceil(i/N))^2*ThTh1(i,j)+I2*Omega^2*LL(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨4,5£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+3*N^2,j+4*N^2)=D12*1/r(1,ceil(i/N))*ThR(i,j)-D22*1/r(1,ceil(i/N))^2*Th1(i,j)+D44*(1/r(1,ceil(i/N))*ThR(i,j)-1/r(1,ceil(i/N))^2*Th1(i,j));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N^2  %¸³Öµ£¨5,1£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+4*N^2,j)=B12*1/r(1,ceil(i/N))*ThR(i,j)+B22*1/r(1,ceil(i/N))^2*Th1(i,j)+B44*(1/r(1,ceil(i/N))^2*Th1(i,j)+1/r(1,ceil(i/N))*ThR(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨5,2£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+4*N^2,j+N^2)=B22*1/r(1,ceil(i/N))^2*ThTh1(i,j)+B44*(RR1(i,j)-1/r(1,ceil(i/N))^2*LL(i,j)+1/r(1,ceil(i/N))*R1(i,j))+I1*Omega^2*LL(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨5,3£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+4*N^2,j+2*N^2)=-Ks*A66*1/r(1,ceil(i/N))*Th1(i,j);
    end
end
for i=1:N^2  %¸³Öµ£¨5,4£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+4*N^2,j+3*N^2)=D12*1/r(1,ceil(i/N))*ThR(i,j)+D22*1/r(1,ceil(i/N))^2*Th1(i,j)+D44*(1/r(1,ceil(i/N))^2*Th1(i,j)+1/r(1,ceil(i/N))*ThR(i,j));
    end
end
for i=1:N^2  %¸³Öµ£¨5,5£©*(N^2*N^2)¾ØÕó¿é
    for j=1:N^2
        K1(i+4*N^2,j+4*N^2)=-Ks*A66*LL(i,j)+D22*1/r(1,ceil(i/N))^2*ThTh1(i,j)+D44*(RR1(i,j)-1/r(1,ceil(i/N))^2*LL(i,j)+1/r(1,ceil(i/N))*R1(i,j))+I2*Omega^2*LL(i,j);
    end
end