function HF=HF79New(x,sparse,exact,fin_dif_2,h)
% Function that computes the Hessian of function 79
% sparse= bool. True= computes the sparse version
% exact= bool. True= computes the exact version, False= computes the approximated version with finite differences
% fin_dif_2= bool. True if exact=false and finite differences are computed using as increment h*abs(x_j) for the derivative with respect to j
% h= increment for the approximated version (if exact=true put h=0)
% JF= function handle of the gradient

n= length(x);
if exact %exact version
    x_next=[x(2:end);0];
    x_prev=[0;x(1:end-1)];
    f=(3-x/10).*x+1-x_prev-2*x_next;
    diag0=5-f/5+(3-x/5).^2;
    diag_1= -2*(3-x(1:end-1)/5)-(3-x_next(1:end-1)/5);
    if sparse % sparse
        diag_up=[0;diag_1] ;
        diag_down=[diag_1 ; 0];
        HF=spdiags([diag_down, diag0,diag_up],[-1,0,1],n,n); %exact version
    else % NOT sparse
        HF=diag(diag0)+diag(diag_1,1)+diag(diag_1,-1); %exact version
    end
else %approximation with finite difference (not exact)
    if fin_dif_2 % version of finite differences with abs(xj)
        z1=zeros(n,1);
        z2=zeros(n,1);
        z3=zeros(n,1);
        z1(1:3:end)=ones(length(1:3:end),1);
        z2(2:3:end)=ones(length(2:3:end),1);
        z3(3:3:end)=ones(length(3:3:end),1);
        z1=h*z1.*abs(x);
        z2=h*z2.*abs(x);
        z3=h*z3.*abs(x);
        F1=(JF(x+z1)-JF(x-z1));
        F2=(JF(x+z2)-JF(x-z2));
        F3=(JF(x+z3)-JF(x-z3));
        F1=[0;F1];
        rem_n=rem(n,3);
        if rem_n==0
            F3=[F3;0];
        elseif rem_n==1
            F1=[F1;0];
        else
            F2=[F2;0];
        end
        m1=reshape(F1,3,[])';
        m2=reshape(F2,3,[])';
        m3=reshape(F3,3,[])';
        M=zeros(n,3);
        M(1:3:end,:)=m1;
        M(2:3:end,:)=m2;
        M(3:3:end,:)=m3;
        for col=1:3
            M(:,col)=M(:,col).*(1./abs(x))/h;
        end
        if sparse % sparse
            HF=spdiags(M,[1,0,-1],n,n);
        else % NOT sparse
            HF=diag(M(:,2))+diag(M(2:end,1),1)+diag(M(1:end-1,1),-1);
        end
    else % classic version of finite differences
        z1=zeros(n,1);
        z2=zeros(n,1);
        z3=zeros(n,1);
        z1(1:3:end)=ones(length(1:3:end),1);
        z2(2:3:end)=ones(length(2:3:end),1);
        z3(3:3:end)=ones(length(3:3:end),1);
        F1=(JF(x+h*z1)-JF(x-h*z1))/(2*h);
        F2=(JF(x+h*z2)-JF(x-h*z2))/(2*h);
        F3=(JF(x+h*z3)-JF(x-h*z3))/(2*h);
        F1=[0;F1];
        rem_n=rem(n,3);
        if rem_n==0
            F3=[F3;0];
        elseif rem_n==1
            F1=[F1;0];
        else
            F2=[F2;0];
        end
        m1=reshape(F1,3,[])';
        m2=reshape(F2,3,[])';
        m3=reshape(F3,3,[])';
        M=zeros(n,3);
        M(1:3:end,:)=m1;
        M(2:3:end,:)=m2;
        M(3:3:end,:)=m3;
        if sparse % sparse
            HF=spdiags(M,[1,0,-1],n,n);
        else % NOT sparse
            HF=diag(M(:,2))+diag(M(2:end,1),1)+diag(M(1:end-1,1),-1);
        end
end


end