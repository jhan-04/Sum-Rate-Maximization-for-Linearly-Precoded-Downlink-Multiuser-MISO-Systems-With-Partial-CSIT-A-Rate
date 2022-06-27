

function [Rs_set]=cal_ach_rate(Nt,M,h,N0,p)


pp=reshape(p,[Nt,M+1]);%p=m=Nt X (M+1)vector


for k=1:M
    inter(k)=0;
    for i=1:M
        if i~=k
            inter(k)=inter(k)+abs(h(:,k)'*pp(:,i))^2;
        end
    end
    R_ach(k)= log2(1+abs(h(:,k)'*pp(:,k))^2/(inter(k)+N0));
    Rc_ach(k)=log2(1+abs(h(:,k)'*pp(:,M+1))^2/(inter(k)+abs(h(:,k)'*pp(:,k))^2+N0));
    
end

R_ach(M+1)=min(Rc_ach);
Rs_set=R_ach.';
end





