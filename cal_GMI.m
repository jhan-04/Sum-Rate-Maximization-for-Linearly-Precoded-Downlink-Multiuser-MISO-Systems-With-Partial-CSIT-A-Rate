function [GMI]=cal_GMI(M,A,B,C,D,p)


for k=1:M
    GMI_p(k)=log2((p'*C(:,:,k)*p)/(p'*D(:,:,k)*p));
    GMI_c(k)=log2((p'*A(:,:,k)*p)/(p'*B(:,:,k)*p));
end

GMI=sum(GMI_p)+min(GMI_c);

end