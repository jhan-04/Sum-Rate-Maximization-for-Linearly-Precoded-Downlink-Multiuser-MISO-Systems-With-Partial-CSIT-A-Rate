function  [ESR]=Alg_1_cal_ESR(pc,pk,H_h,e_m,M,K,N0)

        for m=1:M
            H=H_h+e_m(:,m);
            for k=1:K
                T_ck_mm(k,m)=abs(H(:,k)'*pc)^2+sum(abs(H(:,k)'*pk).^2)+N0;
                T_k_mm(k,m)=sum(abs(H(:,k)'*pk).^2)+N0;
                S_ck_mm(k,m)=abs(H(:,k)'*pc)^2;
                S_k_mm(k,m)=abs(H(:,k)'*pk(:,k)).^2;
                I_ck_mm(k,m)=T_k_mm(k,m);
                I_k_mm(k,m)=I_ck_mm(k,m)-S_k_mm(k,m);
                Rs_ck_mm(k,m)=log2(1+S_ck_mm(k,m)/I_ck_mm(k,m));
                Rs_k_mm(k,m)=log2(1+S_k_mm(k,m)/I_k_mm(k,m));
            end
            Rs_ck(m)=min(Rs_ck_mm(:,m));
            Rs_k(m)= sum(Rs_k_mm(:,m));
        end
        %%%%
        ESR=mean(Rs_ck+Rs_k);

end