function HF=HF16(x,sparse,exact,fin_dif_2,h)
% Function that computes the Hessian of function 79
% sparse= bool. True= computes the sparse version
% exact= bool. True= computes the exact version, False= computes the approximated version with finite differences
% fin_dif_2= bool. True if exact=false and finite differences are computed using as increment h*abs(x_j) for the derivative with respect to j
% h= increment for the approximated version (if exact=true put h=0)

n=length(x);
indices=(1:n)';
D=indices.*cos(x)-2*sin(x);
if sparse %sparse version
     HF=spdiags(D,0,n,n);%exact version
    if ~exact %approximation with finite difference (not exact)
        if fin_dif_2 % version of finite differences with abs(xj)
           %vec_diag=((1-cos(h*abs(x)))./(x.^2))/(h^2); % no taylor expansion
           vec_diag= 1-(h^2)*(x.^2)/12 ; % with taylor expansion
           new_diag=vec_diag.*diag(HF,0);
           HF=spdiags(new_diag,0,n,n); 
        else % classic version of finite differences
           HF=(1-(h^2)/12)*HF ;  % as h is small i use cos(h) taylor expansion to avoid numerical cancellation
        end
        % elementi non diagonali nulli anche con differenze finite
    end
else % NOT sparse version
    HF=diag(D); %exact version
    if ~exact %approximation with finite difference (not exact)
       if fin_dif_2 % version of finite differences with abs(xj)
           % vec_diag=((1-cos(h*abs(x)))./(x.^2))/(h^2); % no taylor expansion
           vec_diag= 1-(h^2)*(x.^2)/12 ; % with taylor expansion
           HF=diag(vec_diag.*diag(HF)); 
       else % classic version of finite differences
           HF=(1-(h^2)/12)*HF; % as h is small i use cos(h) taylor expansion to avoid numerical cancellation
       end
       % elementi non diagonali nulli anche con differenze finite
    end

end

HF(n,n)=(n-1)*sin(x(n))+n*cos(x(n));
if ~exact  %approximation with finite difference (not exact)
    HF(n,n)=(1-(h^2)/12)*HF(n,n);
end


end