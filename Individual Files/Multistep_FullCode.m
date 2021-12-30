
function [u] = multistep_geral(a,b,b_1,t0,T,dt,f,y0,max_it,tol,dfy)
    q = length(a);
    time = t0:dt:T;
    u(:,1) = y0;
    p=q-1;
    for i = 1:q-1 
        u(:,i+1) = heun(time,i,dt,f,u(:,i));
    end

    for i = (q+1):length(time)
        sizeu = length(u);
        soma = 0;
        for j=1:p+1
            soma = soma + a(j)*u(:,i-j) + dt * b(j)*f(time(i-j),u(:,i-1));
        end
        if b_1==0
            u(:,i) = soma;
        else
            Func = @(soma,dt,b_1,f,time,i,x) x-soma-dt*b_1*(f(time(i)+dt,x));
            u(:,i)=Newton_(Func,u(:,i-1),f,dfy,time,i,b_1,dt,soma,max_it,tol);
        end
    end
    
    function [uj_next] = heun(time,j,dt,f,uj)
        uj_next=uj+dt/2*(f(time(j),uj)+f(time(j+1),uj+dt*f(time(j),uj)));
    end

    function [x] = Newton_(Func,x,f,dfy,time,i,b_1,dt,soma,max_it,tol)
        count = 1; erro=10;
        while erro>tol && count < max_it 
            val=(Func(soma,dt,b_1,f,time,i,x))';
            J = eye(size(x))-dt*b_1*dfy(time(i)+dt,x);
            df = (-J/val)';
            x=x+df';
            erro = norm(df,Inf);
            count = count+1;
            if count == max_it
                disp('max iterações, rever params')
            end
        end
    end
end
