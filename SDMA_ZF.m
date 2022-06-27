function [P]=SDMA_ZF(P,h1,h2,rho)


Nt=length(h1);
% mu=P/2+(1/(2*rho))*(1/norm(h2)^2+1/norm(h1)^2);
%  P1_p=mu-(1/rho)*(1/norm(h1)^2);
%  P2_p=mu-(1/rho)*(1/norm(h2)^2);
 a_max=100000;
 a_min=-100000;
 while (a_max-a_min)>0.0001
     a=(a_max+a_min)/2;
      P1_x=max(a-(1/rho)*(1/norm(h1)^2),0);
     P2_x=max(a-(1/rho)*(1/norm(h2)^2),0);
     
     if (P1_x+P2_x)>P
         a_max=a;
     else
         a_min=a;
     end
     
 end
 
 P=[P1_x;P2_x];
%Rs_x=log2(1+norm(h1)^2*rho*P1_x)+log2(1+norm(h2)^2*rho*P2_x);

    
end