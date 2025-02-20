function HF=HF16New(x,sparse,exact,fin_dif_2,h,JF)
% Function that computes the Hessian of function 79
% sparse= bool. True= computes the sparse version
% exact= bool. True= computes the exact version, False= computes the approximated version with finite differences
% fin_dif_2= bool. True if exact=false and finite differences are computed using as increment h*abs(x_j) for the derivative with respect to j
% h= increment for the approximated version (if exact=true put h=0)
% JF= function handle of the gradient

n=length(x);
if exact %exact version
    indices=(1:n)';
    D=indices.*cos(x)-2*sin(x);
    if sparse %sparse version
        HF=spdiags(D,0,n,n);
    else % NOT sparse version
        HF=diag(D); 
    end
    HF(n,n)=(n-1)*sin(x(n))+n*cos(x(n));
else  %approximation with finite difference (not exact)
    if fin_dif_2 % version of finite differences with abs(xj)
          d=(JF(x+h.*abs(x).*ones(n,1))-JF(x-h.*abs(x).*ones(n,1)))./(2*h*abs(x));
          if sparse %sparse version
              HF =spdiags(d,0,n,n); 
          else % NOT sparse
              HF=diag(d);
          end
    else % classic version of finite differences
          d=(JF(x+h*ones(n,1))-JF(x-h*ones(n,1)))/(2*h);
          if sparse %sparse version
              HF =spdiags(d,0,n,n); 
          else % NOT sparse
              HF=diag(d);
          end
    end
end

end