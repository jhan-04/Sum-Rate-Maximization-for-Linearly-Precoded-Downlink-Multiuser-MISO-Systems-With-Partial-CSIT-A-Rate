clear all

M=100;%sample
K=2;%users
Nt=3;%Transmitter antenna
sigma_2=0.063;
sigma_e1=sqrt(sigma_2);
sigma_e=zeros(K)+sigma_e1;%all cahannel has sam channel errror covariance
SNR=5:6:35;

alpha=0.1;
sample=2;
for i=1:sample
    i
    H_h1=(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);
    %H=(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);
    %e=sqrt(sigma_2)*(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);
    %H_h=H-e;
    
    
    
    for j=1:1:length(SNR)
        snr=SNR(j)
        Pt=10^(snr/10);
        sigma_2=Pt^(-alpha);
        H_h=sqrt(1-sigma_2)*H_h1;
        N0=1;
        %%%
        for k=1:K
            pk(:,k)=sqrt(Pt^(-alpha))/K*(H_h(:,k)/norm(H_h(:,k)));
        end
       % pk=sqrt(Pt^(-alpha))*(H_h/norm(H_h));
        [X,Y,Z]=svd(H_h);
        pc=sqrt(Pt-Pt^(-alpha))*X(:,1)/norm(X(:,1));
        
        e_m=sqrt(sigma_2)*(randn(Nt,K,M)+1i*randn(Nt,K,M))/sqrt(2);
        H_m=H_h+e_m;
         
        %%%MMSE
        [pc_1,pk_1]=Alg_1(Nt,K,pc,pk,H_m,Pt,M,N0);
        ESR(i,j)=Alg_1_cal_ESR(pc_1,pk_1,H_h,e_m,M,K,N0);%paper
%         [pc_noRS,pk_noRS]=Alg_1_noRS(Nt,K,pc,pk,H_m,Pt,M,N0);
%         ESR_noRS(i,j)=Alg_1_cal_ESR(pc_noRS,pk_noRS,H_h,e_m,M,K,N0);%paper
        
        [pc_con,pk_con,Bound(i,j)]=Alg_1_con(Nt,K,pc,pk,H_h,Pt,N0,sigma_2);
        ESR_con(i,j)=Alg_1_cal_ESR(pc_con,pk_con,H_h,e_m,M,K,N0);%paper
%         
%         
%         %%%GMI
%         %RS
%         [pc_G,pk_G,GMI(i,j),SR1(i,j)]=GMI_RS(Nt,K,H_h,H_m,Pt,M,N0,sigma_e);%proposed
%         SR2(i,j)=Alg_1_cal_ESR(pc_G,pk_G,H_h,e_m,M,K,N0);
%         %SDMA
%         [pc_SDMA,pk_SDMA,GMI_nc(i,j),SR_SDMA(i,j)]=GMI_SDMA(Nt,K,H_h,H_m,Pt,M,N0,sigma_e);
%         
%         
%         %%%no-RS-ZF
%         [pc_ZF,pk_ZF,SR_ZF(i,j)]=ZF_noRS(Nt,K,H_h,H_m,Pt,M,N0);
        
        
        
    end
    
end
hold on
plot(SNR,mean(Bound),'k-','LineWidth',2,...
    'MarkerSize',8)

plot(SNR,mean(ESR),'g--o','LineWidth',2,...
    'MarkerSize',7)
plot(SNR,mean(ESR_con),'k--o','LineWidth',2,...
    'MarkerSize',7)
grid


legend('Lower Bound-MMSE','RS-MMSE','RS-MMSE-using-bound')