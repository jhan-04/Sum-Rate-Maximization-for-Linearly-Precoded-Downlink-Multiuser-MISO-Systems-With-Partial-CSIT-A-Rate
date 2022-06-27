function  [A,B,C,D]=cal_ABCD(Nt,M,Pt,N0,h,sigma_e)




%%%%A_k
for k=1:M
    H(:,:,k)=h(:,k)*h(:,k)';
    for j=1:M+1
        a(Nt*(j-1)+1:j*Nt,Nt*(j-1)+1:j*Nt,k)=H(:,:,k);
    end
    %A(:,:,k)=[]
    A(:,:,k)=a(:,:,k)+sigma_e(k)^2*eye(Nt*(M+1))+(N0/Pt)*eye(Nt*(M+1));
end

%%%%B_k
for k=1:M
    %H(:,:,k)=h(:,k)*h(:,k)';
    for j=1:M
        b(Nt*(j-1)+1:j*Nt,Nt*(j-1)+1:j*Nt,k)=H(:,:,k);
    end
    b(Nt*(M+1),Nt*(M+1),k)=0;
    %A(:,:,k)=[]
    B(:,:,k)=b(:,:,k)+sigma_e(k)^2*eye(Nt*(M+1))+(N0/Pt)*eye(Nt*(M+1));
end

%%%%C_k
for k=1:M
    %H(:,:,k)=h(:,k)*h(:,k)';
    for j=1:M
        c(Nt*(j-1)+1:j*Nt,Nt*(j-1)+1:j*Nt,k)=H(:,:,k)+sigma_e(k)^2*eye(Nt);
    end
    c(Nt*(M+1),Nt*(M+1),k)=0;
    %A(:,:,k)=[]
    C(:,:,k)=c(:,:,k)+(N0/Pt)*eye(Nt*(M+1));
end


%%%%D_k
for k=1:M
    %H(:,:,k)=h(:,k)*h(:,k)';
    for j=1:M
        
        d(Nt*(j-1)+1:j*Nt,Nt*(j-1)+1:j*Nt,k)=H(:,:,k)+sigma_e(k)^2*eye(Nt);
        
        if j==k
            d(Nt*(j-1)+1:j*Nt,Nt*(j-1)+1:j*Nt,k)=sigma_e(k)^2*eye(Nt);
        end
        
    end
    d(Nt*(M+1),Nt*(M+1),k)=0;
    %A(:,:,k)=[]
    D(:,:,k)=d(:,:,k)+(N0/Pt)*eye(Nt*(M+1));
end


for k=1:M
    A(:,:,k)=(A(:,:,k)+A(:,:,k)')/2;
    B(:,:,k)=(B(:,:,k)+B(:,:,k)')/2;
    C(:,:,k)=(C(:,:,k)+C(:,:,k)')/2;
    D(:,:,k)=(D(:,:,k)+D(:,:,k)')/2;
end


end
