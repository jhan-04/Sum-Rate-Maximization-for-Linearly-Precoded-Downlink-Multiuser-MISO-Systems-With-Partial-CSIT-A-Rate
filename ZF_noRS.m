function  [pc,pk,SR_ZF]=ZF_noRS(Nt,K,H_h,H_m,Pt,M,N0)
%this function only work in two-user system 

rho=1-abs(H_h(:,1)'/norm(H_h(:,1))*H_h(:,2)/norm(H_h(:,2)))^2;
[P_ZF]=SDMA_ZF(Pt,H_h(:,1),H_h(:,2),rho);
pk_zf=H_h*inv(H_h'*H_h);
for k=1:K
    mag(k)=norm(pk_zf(:,k));
end

ZF_p=[reshape(pk_zf.*(sqrt(P_ZF)./mag')',[],1);zeros(Nt,1)];
%[GMI_ZF]=cal_GMI(K,A,B,C,D,ZF_p);
pc=zeros(1,Nt);
pk=pk_zf.*(sqrt(P_ZF)./mag')';
ZF_p;
for m=1:M
    [Rs_set_ZF]=cal_ach_rate(Nt,K,H_m(:,:,m),N0,ZF_p);
    ach_rate_mrt_ZF(m)=sum(Rs_set_ZF);
end
SR_ZF=mean(ach_rate_mrt_ZF);
end