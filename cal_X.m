function  [X,result_set]=cal_X(Nt,M,Pt,A,B,C,D)



%initialization
%Xs(:,:,1)=eye(Nt*(M+1),Nt*(M+1));
%a0(1)=0;
 b0(:,1)=(rand(M,1));%randn(M,1)*0+5;%
% %c0(1)=0;
 d0(:,1)=(rand(M,1));%randn(M,1)*0+5;%
%bb=(rand(Nt*(M+1),1));%randn(M,1)*0+5;%
%c0(1)=0;
%dd=(rand(Nt*(M+1),1));%randn(M,1)*0+5;%


s=1;
result(1)=0;

stop=0;
ee=0.01;

%repeat
while ~stop
    
    s=s+1;
    if b0(1,s-1)~=b0(1,s-1)
        b0(:,s-1)=(rand(M,1));
        d0(:,s-1)=(rand(M,1));
    end
    
    cvx_begin quiet
    %cvx_begin
    variable X(Nt*(M+1),Nt*(M+1)) complex semidefinite
    variable Lc(1) nonnegative
    variable a(M,1) nonnegative
    variable b(M,1) nonnegative
    variable c(M,1) nonnegative
    variable d(M,1) nonnegative
    xx=Lc+sum(c)-sum(d);
    maximize(xx)
    subject to
    for i=1:M
        a(i)-b(i)>=Lc;
        trace(A(:,:,i)*X)>=exp(a(i));
        trace(C(:,:,i)*X)>=exp(c(i));
        %linearize
        %         i
        %         exp(b0(i,s-1))
        %         -b0(i,s-1)+1
        trace(B(:,:,i)*X)<=exp(b0(i,s-1))*(b(i)-b0(i,s-1)+1);
        trace(D(:,:,i)*X)<=exp(d0(i,s-1))*(d(i)-d0(i,s-1)+1);
    end
    %trace(X)<=Pt;
    X== hermitian_semidefinite(Nt*(M+1));
    
    cvx_end
    
    
    %zzz=[a b c d]
    
    b0(:,s)=b;
    d0(:,s)=d;
    result(s)=Lc+sum(c)-sum(d);
    result_Lc(s)=Lc;
    result_c(s)=sum(c);
    result_d(s)=sum(d);
    
    % until
    if norm(b0(:,s)-b0(:,s-1))^2+norm(d0(:,s)-d0(:,s-1))^2<ee%norm(result(s)-result(s-1))^2<0.0001
        stop=1;
    end
    
    if result(s-1)>result(s)&&s>2
        fprintf('NOT monotonically increasing %d in Rs \n',result(s-1)-result(s))
        result
    end
    %
    %[GMI(s)]=cal_GMI_withX(M,A,B,C,D,X);
    
end

% %GMI
% %result/log(2)
result_final=result(s);
result_set=result.';

end