function  [pc,pk,Bound]=Alg_1_con(Nt,K,pc,pk,H_h,Pt,N0,sigma_2)
%%%
n=2;
A=[];
A(2)=0;
A(1)=10;
e=0.01;
while abs(A(n)-A(n-1))>=e
    
    
    for k=1:K
        for i=1:K
            ss(i)=pk(:,i)'*(H_h(:,k)*H_h(:,k)'+sigma_2*eye(Nt))*pk(:,i);
        end
        T_k(k)=real(sum(ss)+N0);
        T_ck(k)=real(pc'*(H_h(:,k)*H_h(:,k)'+sigma_2*eye(Nt))*pc+T_k(k));
        
        %             S_ck(k)=abs(H_h(:,k)'*pc)^2;
        %             S_k(k)=abs(H_h(:,k)'*pk(:,k)).^2;
        %             I_ck(k)=T_k(k);
        %             I_k(k)=I_ck(k)-S_k(k);
        %%%
        g_ck(k)=pc'*H_h(:,k)/T_ck(k);
        g_k(k)=pk(:,k)'*H_h(:,k)/T_k(k);
        u_ck(k)=(1-abs(H_h(:,k)'*pc)^2/T_ck(k))^(-1);
        u_k(k)=(1-abs(H_h(:,k)'*pk(:,k))^2/T_k(k))^(-1);
        %%%
        t_ck(k)=u_ck(k)*abs(g_ck(k))^2;
        t_k(k)=u_k(k)*abs(g_k(k))^2;
        psi_ck(:,:,k)=((t_ck(k)*H_h(:,k)*H_h(:,k)')+(t_ck(k)*H_h(:,k)*H_h(:,k)')')/2;
        psi_k(:,:,k)=((t_k(k)*H_h(:,k)*H_h(:,k)')+(t_k(k)*H_h(:,k)*H_h(:,k)')')/2;
        f_ck(:,k)=u_ck(k)*H_h(:,k)*g_ck(k)';
        f_k(:,k)=u_k(k)*H_h(:,k)*g_k(k)';
        v_ck(k)=log2(u_ck(k));
        v_k(k)=log2(u_k(k));
    end
    
    
    
    clearvars a c ac
    
    cvx_begin quiet
    %cvx_begin
    variable xi_c(1,1) 
    variable P(Nt*(K+1),1) complex
    % variable a(K,K) complex
    % variable c(K)
    % variable ac(K,K) complex
    % variable cc(K)
    % PP=reshape(P,Nt,K+1);
    % PK=PP(:,2:K+1);
    % Pc=PP(:,1);
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
    
end
A
for k=1:K
    for i=1:K
        ss(i)=pk(:,i)'*(H_h(:,k)*H_h(:,k)'+sigma_2*eye(Nt))*pk(:,i);
    end
    T_k(k)=real(sum(ss)+N0);
    T_ck(k)=real(pc'*(H_h(:,k)*H_h(:,k)'+sigma_2*eye(Nt))*pc+T_k(k));
    R_ck(k)=-log2(1-abs(H_h(:,k)'*pc)^2/T_ck(k));
    R_k(k)=-log2(1-abs(H_h(:,k)'*pk(:,k))^2/T_k(k));
end
Bound=min(R_ck)+sum(R_k);
end