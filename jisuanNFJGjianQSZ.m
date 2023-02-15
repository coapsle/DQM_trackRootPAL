qsz=24;
s2=zeros(70,2);
for i=1:(size(s1)-qsz)
    s2(i,1)=s1(i+qsz,1);
    s2(i,2)=s1(i+qsz,2);
end