function HF=HF16(x,sparse,exact,h)
% Function that computes the Hessian of function 79
% sparse= bool. True= computes the sparse version
% exact= bool. True= computes the exact version, False= computes the approximated version with finite differences
% h= increment for the approximated version (if exact=true put h=0)

n=length(x);
indices=(1:n)';
D=indices.*cos(x)-2*sin(x);
if sparse %sparse version
     HF=spdiags(D,0,n,n);%exact version
    if ~exact %approximation with finite difference (not exact)
        HF=(1-(h^2)/12)*HF ; 
        % elementi non diagonali nulli anche con differenze finite
    end
else % NOT sparse version
    HF=diag(D); %exact version
    if ~exact %approximation with finite difference (not exact)
        HF=(1-(h^2)/12)*HF;  
        % elementi non diagonali nulli anche con differenze finite
    end

end

HF(n,n)=(n-1)*sin(x(n))+n*cos(x(n));
if ~exact  %approximation with finite difference (not exact)
    HF(n,n)=(1-(h^2)/12)*HF(n,n);
end


end