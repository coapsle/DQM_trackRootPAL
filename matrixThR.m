%Xi;
%N=15;
ThR=zeros(N^2);
% for i=1:N
%     for j=1:N
%         if j~=i
%             ThR(i,j)=M1r(1,i)/((r(1,i)-r(1,j))*M1r(1,j));
%         else
%             ThR(i,j)=0;
%         end
%     end
% end
for i=1:N
    for j=1:N
        for k=1:N
            for l=1:N
                ThR((i-1)*N+j,(k-1)*N+l)=C1(i,k)*G1(j,l);
            end
        end
    end
end