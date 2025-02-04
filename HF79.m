function HF=HF79(x,sparse,exact,fin_dif_2,h)
% Function that computes the Hessian of function 79
% sparse= bool. True= computes the sparse version
% exact= bool. True= computes the exact version, False= computes the approximated version with finite differences
% fin_dif_2= bool. True if exact=false and finite differences are computed using as increment h*abs(x_j) for the derivative with respect to j
% h= increment for the approximated version (if exact=true put h=0)

x_next=[x(2:end);0];
x_prev=[0;x(1:end-1)];
n= length(x);

f=(3-x/10).*x+1-x_prev-2*x_next;
diag0=5-f/5+(3-x/5).^2;
diag_1= -2*(3-x(1:end-1)/5)-(3-x_next(1:end-1)/5);
if sparse % sparse
    diag_up=[0;diag_1] ;
    diag_down=[diag_1 ; 0];
    HF=spdiags([diag_down, diag0,diag_up],[-1,0,1],n,n); %exact version
    if ~exact  %approximation with finite difference (not exact)
        if fin_dif_2 % version of finite differences with abs(xj)
            HF= HF + (h^2)*spdiags(x.^2,0,n,n)/100 ;
            vec_diag= abs(x(1:end-1))/5 + abs(x(2:end))/10;
            HF= HF + spdiags(h*[0;vec_diag],1,n,n) +spdiags(h*[vec_diag;0],-1,n,n);
            % SE NON SERVONO LE RIGHE SOTTO ELIMINA
%             d_0=diag(HF,0);
%             mult_diag=abs(x(1:end-1)).*abs(x(2:end));
%             d_1=mult_diag.*diag(HF,1);
%             HF= spdiags(d_0,0,n,n)+spdiags([0;d_1],1,n,n)+spdiags([d_1;0],-1,n,n);
            
        else % classic version of finite differences
            HF= HF + spdiags((h^2)*ones(n,1)/100,0,n,n);
            new_coef=3*h/10;
            HF= HF + spdiags(new_coef*ones(n,1),1,n,n) +spdiags(new_coef*ones(n,1),-1,n,n); 
        end
        
    end
    
else % NOT sparse
    HF=diag(diag0)+diag(diag_1,1)+diag(diag_1,-1); %exact version
    if ~exact  %approximation with finite difference (not exact)

       if fin_dif_2 % version of finite differences with abs(xj)
           HF= HF + (h^2)*diag(x.^2)/100 ;
           vec_diag= abs(x(1:end-1))/5 + abs(x(2:end))/10;
           HF= HF + diag(h*vec_diag,1) +diag(h*vec_diag,-1);
           % SE NON SERVONO LE RIGHE SOTTO ELIMINA
%            d_0=diag(HF);
%            mult_diag=abs(x(1:end-1)).*abs(x(2:end));
%            d_1=mult_diag.*diag(HF,1);
%            HF= diag(d_0)+diag(d_1,1)+diag(d_1,-1);
           
       else % classic version of finite differences
           HF= HF + (h^2)*eye(n)/100 ;
           new_coef=3*h/10;
           HF= HF + diag(new_coef*ones(n-1,1),1) +diag(new_coef*ones(n-1,1),-1);
    end
end

end