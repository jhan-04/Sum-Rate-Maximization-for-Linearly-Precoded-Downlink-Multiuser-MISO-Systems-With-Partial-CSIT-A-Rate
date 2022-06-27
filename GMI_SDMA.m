function  [pc,pk,GMI_nc,SR_SDMA]=GMI_SDMA(Nt,K,H_h,H_m,Pt,M,N0,sigma_e)


[A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
[X_nc,result_nc]=cal_X_no_common(Nt,K,Pt,A,B,C,D);
%result_nc_set(1:length(result_nc))=result_nc;
%[GMI_X_nc]=cal_GMI_withX(K,A,B,C,D,X_nc);
[p_nc]=find_p(K,Pt,A,B,C,D,X_nc);norm(p_nc)^2;
[GMI_nc]=cal_GMI(K,A,B,C,D,p_nc);

pc=p_nc(Nt*K+1:Nt*(K+1));
pk=reshape(p_nc(1:Nt*K),Nt,K);
for m=1:M
    [Rs_set]=cal_ach_rate(Nt,K,H_m(:,:,m),N0,p_nc);
    ach_rate_set(m)=sum(Rs_set);
end

SR_SDMA=mean(ach_rate_set);
end