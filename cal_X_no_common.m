function  [X_nc,result_set]=cal_X_no_common(Nt,M,Pt,A,B,C,D)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%initialization
%Xs(:,:,1)=eye(Nt*(M+1),Nt*(M+1));
%a0(1)=0;
%b0(:,1)=abs(randn(M,1))*3;
%c0(1)=0;
d0(:,1)=rand(M,1);

s=1;
result(1)=0;

stop=0;
ee=0.01;
for k=1:M
    CC(:,:,k)=C(1:Nt*M,1:Nt*M,k);
    DD(:,:,k)=D(1:Nt*M,1:Nt*M,k);
end

%repeat
while ~stop
    
    s=s+1;
    if d0(1,s-1)~=d0(1,s-1)
        d0(:,s-1)=(rand(M,1));
    end
    
    cvx_begin quiet
    %cvx_begin
    %variable X(Nt*(M+1),Nt*(M+1)) hermitian
    variable X_nc(Nt*(M),Nt*(M)) hermitian
    variable Lc(1) nonnegative
    %     variable a(M,1) nonnegative
    %     variable b(M,1) nonnegative
    variable c(M,1) nonnegative
    variable d(M,1) nonnegative
    xx=sum(c)-sum(d);
    maximize(xx)
    
    subject to
    for i=1:M
        %a(i)-b(i)>=Lc
        %trace(A(:,:,i)*X)>=exp(a(i))
        trace(CC(:,:,i)*X_nc)>=exp(c(i))
        %linearize
        %trace(B(:,:,i)*X)<=exp(b0(i,s-1))*(b(i)-b0(i,s-1)+1)
        trace(DD(:,:,i)*X_nc)<=exp(d0(i,s-1))*(d(i)-d0(i,s-1)+1)
    end
    %trace(X_nc)<=Pt;
    X_nc==hermitian_semidefinite(Nt*(M));
    
    cvx_end
    
    %     b0(:,s)=b;
    d0(:,s)=d;
    
    %     result(s)=Lc+sum(c)-sum(d);
    %     result_Lc(s)=Lc;
    %     result_c(s)=sum(c);
    %     result_d(s)=sum(d);
    
    % until
    if norm(d0(:,s)-d0(:,s-1))^2<ee %norm(result(s)-result(s-1))^2<0.0001
        stop=1;
    end
    
    
    result(s)=sum(c)-sum(d);
    if result(s-1)>result(s)&&s>2
        fprintf('NOT monotonically increasing %d in No- Rs \n',result(s-1)-result(s))
        result
    end
    X_nc1=X_nc;
    X_nc(Nt*(M+1),Nt*(M+1))=0;
    
    
    
    for k=1:M
        GMI_p(k)=log2(trace(C(:,:,k)*X_nc))-log2((trace(D(:,:,k)*X_nc)));
        %GMI_c(k)=log2(trace(A(:,:,k)*X_nc)/(trace(B(:,:,k)*X_nc)));
    end
    
    GMI(s)=sum(GMI_p);
    
    
    
end

GMI;
result/log(2);
result_final=result(s);
result_set=result.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end