
function [ord] = give_convergence_ms(f,dfy,y_real,y0,t0,T,a,b,b_1,tol,itmax)
    for j = 1:5
        dt = (0.1)/(2^(j-1));
        time = t0:dt:T;
        u_real = y_real(time);
        u_method = multistep_geral(a,b,b_1,t0,T,dt,f,y0,5000,1e-4,dfy);
        for i=1:length(time)-1
            error_method(i) = norm(u_real(i)-u_method(i),'inf');
        end
        error(j) = max(error_method);
        if j>1
            ord(j-1)=log2(error(j-1)/error(j));
        end
    end
end
