clear all
clc
close all

%funcao = @(t) t.^3 - t
f = @(t,y) [y(2)-funcao(y(1)); -y(1)];
t0 = 0;
T = 5;
dt = 0.1;
y_real = @(t);
dfy  = @(t,y) [-3*y(1).^2; 0]'
y0 = [1; 1];
a = [18/11 -9/11 2/11 0 0 0] ;
b = [0 0 0 0 0 0];
b_1 = 6/11;
u_1 = multistep_geral(a,b,b_1,t0,T,dt,f,y0,10000,1e-4,dfy)
ordemconv = give_convergence_ms(f,dfy,y_real,y0,t0,T,a,b,b_1,1e-6,5000)

%plot da curva aproximada e da solução exata

% figure
% plot(t0:dt:T,y_real(t0:dt:T),'Linewidth',2, 'color', 'r')
% hold on
% plot(t0:dt:T,u_1,'--', 'Linewidth', 2, 'color','b')
% title('Adam Moulton and Exact Solution')
% legend('Solução Exata','Aproximação Numérica')
% grid on;
% hold off


%plot do erro

% figure
% plot(t0:dt:T,(y_real(t0:dt:T)-u_1),'Linewidth',2, 'color', 'g')
% title('Erro')
% grid on;


%gráfico de fase - caso van der pol

% for j=1:length(t0:dt:T)
%     y1_dot_aprox(j)=u_1(1,j);
%     y2_dot_aprox(j)=u_1(2,j);
% end
% plot(y1_dot_aprox,y2_dot_aprox,'r')