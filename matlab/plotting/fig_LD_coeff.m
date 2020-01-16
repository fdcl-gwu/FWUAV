clear all;
close all;

N=501;
alpha=linspace(0,pi/2,N);
for k=1:N
    [C_L(k) C_D(k)]=wing_QS_LD_coeff(alpha(k));
end

plot(alpha,C_L,'r', alpha, C_D, 'b');
xlabel('$\alpha$','interpreter','latex');
ylabel('$C_L,C_D$','interpreter','latex');
hl=legend({'$C_D$', '$C_D$'});
set(hl,'interpreter','latex');



