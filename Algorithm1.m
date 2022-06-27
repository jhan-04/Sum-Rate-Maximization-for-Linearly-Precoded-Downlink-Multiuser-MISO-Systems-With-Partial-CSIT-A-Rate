clear all

M=100;%sample
K=2;%users
Nt=3;%Transmitter antenna
sigma_2=0.063;

snr=10;

for hn=1:1
    hn
    H_h1=(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);
    %H=(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);
    %     e=sqrt(sigma_2)*(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);
    %     H=H_h+e;
    %H_h=H-e;
    % H_m=H_h+e_m;
    
    
    
    for s=1:1:length(snr)
                SNR=snr(s);
        Pt=10^(SNR/10);
        sigma_2=Pt^(-0.1);
        H_h=sqrt(1-sigma_2)*H_h1;

        N0=1;
        %%%
        pk=Pt^(-0.1)*(H_h/norm(H_h));
        [X,Y,Z]=svd(H_h);
        pc=(Pt-Pt^(-0.1))*X(:,1)/norm(X(:,1));
        %%%
        n=2;
        A=[];
        A(2)=0;
        A(1)=10;
        e=0.01;
        while abs(A(n)-A(n-1))>=e
            
            for m=1:M
                e_m(:,:,m)=sqrt(sigma_2)*(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);
                H_m(:,:,m)=H_h+e_m(:,m);
                
                
                for k=1:K
                    T_ck(k,m)=abs(H_m(:,k,m)'*pc)^2+sum(abs(H_m(:,k,m)'*pk).^2)+N0;
                    T_k(k,m)=sum(abs(H_m(:,k,m)'*pk).^2)+N0;
                    S_ck(k,m)=abs(H_m(:,k,m)'*pc)^2;
                    S_k(k,m)=abs(H_m(:,k,m)'*pk(:,k)).^2;
                    I_ck(k,m)=T_k(k,m);
                    I_k(k,m)=I_ck(k,m)-S_k(k,m);
                    %%%
                    g_ck_m(k,m)=pc'*H_m(:,k,m)/T_ck(k,m);
                    g_k_m(k,m)=pk(:,k)'*H_m(:,k,m)/T_k(k,m);
                    u_ck_m(k,m)=T_ck(k,m)/I_ck(k,m);
                    u_k_m(k,m)=T_k(k,m)/I_k(k,m);
                    %%%
                    t_ck_m(k,m)=u_ck_m(k,m)*abs(g_ck_m(k,m))^2;
                    t_k_m(k,m)=u_k_m(k,m)*abs(g_k_m(k,m))^2;
                    psi_ck_m(:,:,k,m)=t_ck_m(k,m)*H_m(:,k,m)*H_m(:,k,m)';
                    psi_k_m(:,:,k,m)=t_k_m(k,m)*H_m(:,k,m)*H_m(:,k,m)';
                    f_ck_m(:,k,m)=u_ck_m(k,m)*H_m(:,k,m)*g_ck_m(k,m)';
                    f_k_m(:,k,m)=u_k_m(k,m)*H_m(:,k,m)*g_k_m(k,m)';
                    v_ck_m(k,m)=log2(u_ck_m(k,m));
                    v_k_m(k,m)=log2(u_k_m(k,m));
                end
            end
            %averaging=bar{variable in paper}
            for k=1:K
                t_ck(k)=mean(t_ck_m(k,:));
                t_k(k)=mean(t_k_m(k,:));
                psi_ck(:,:,k)=round(mean(psi_ck_m(:,:,k,:),4),13);
                psi_k(:,:,k)=(mean(psi_k_m(:,:,k,:),4)+(mean(psi_k_m(:,:,k,:),4))')/2;
                f_ck(:,k)=mean(f_ck_m(:,k,:),3);
                f_k(:,k)=mean(f_k_m(:,k,:),3);
                v_ck(k)=mean(v_ck_m(k,:));
                v_k(k)=mean(v_k_m(k,:));
                %%%
                u_ck(k)=mean(u_ck_m(k,:));
                u_k(k)=mean(u_k_m(k,:));
            end
            
            clearvars a c ac
            
            cvx_begin quiet
            %cvx_begin
            variable xi_c(1,1) 
            variable P(Nt*(K+1),1) complex
            for k=1:K
                for i=1:K
                    a(k,i)=P((Nt*i+1):Nt*(i+1),1)'*psi_k(:,:,k)*P((Nt*i+1):Nt*(i+1),1);
                end
                c(k)=real(f_k(:,k)'*P((Nt*k+1):Nt*(k+1),1));
            end
            
            xx=xi_c+sum(sum(a))+sum(N0*t_k+u_k-v_k)+sum(-2.*c);
            minimize(xx)
            subject to
            for k=1:K
                for i=1:k
                    ac(k,i)=P(Nt*i+1:Nt*(i+1),1)'*psi_ck(:,:,k)*P((Nt*i+1):Nt*(i+1),1);
                end
                P(1:Nt,1)'*psi_ck(:,:,k)*P(1:Nt,1)+sum(ac(k,:))+N0*t_ck-2*real(f_ck(:,k)'*P(1:Nt,1))+u_ck-v_ck<=xi_c;
            end
            P'*P<=Pt;
            cvx_end
            
            %%%
            PP=reshape(P,Nt,K+1);
            pc=PP(:,1);
            pk=PP(:,2:K+1);
            
            n=n+1;
            A(n)=xx;
%             P'*P
%             Pt
        end
        A
        %%%
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
        ESR(hn,s)=mean(Rs_ck+Rs_k);
        
    end
    
end