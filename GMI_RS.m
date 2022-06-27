function  [pc,pk,GMI,SR]=GMI_RS(Nt,K,H_h,H_m,Pt,M,N0,sigma_e)


[A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
[X,result]=cal_X(Nt,K,Pt,A,B,C,D);
%result_set(1:length(result))=result;
[p_max]=find_p(K,Pt,A,B,C,D,X);%norm(p_max(:,i,j))^2;
pc=p_max(Nt*K+1:Nt*(K+1));
pk=reshape(p_max(1:Nt*K),Nt,K);
[GMI]=cal_GMI(K,A,B,C,D,p_max);
for m=1:M
    [Rs_set]=cal_ach_rate(Nt,K,H_m(:,:,m),N0,p_max);
    ach_rate_set(m)=sum(Rs_set);
end

SR=mean(ach_rate_set);
end