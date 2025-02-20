
    % Initialization of zj and j
    zj = z0; 
    j= 0; 
    
    if ~exact %approximation with finite difference (not exact)
        if fin_dif_2
            Azj= (gradf(xk+h.*abs(xk).*zj)-gradk)./(h*abs(xk));
        else
            Azj= (gradf(xk+h*zj)-gradk)/h;
        end
    end

       if ~exact %approximation with finite difference (not exact)
           
           if fin_dif_2
               z= (gradf(xk+h.*abs(xk).*p)-gradk)./(h*abs(xk));
           else
                z=(gradf(xk+h*p)-gradk)/h;
           end
       end
      
       if ~exact %approximation with finite difference (not exact)
            if fin_dif_2  
                z_new= ((gradf(xk+h.*abs(xk).*p)-gradk)./(h*abs(xk)))';
            else
                z_new=((gradf(xk+h*p)-gradk)/h)';
            end
       end