function [pow, E, E_dot, eff] = compute_power(m, t, x, x_dot, tau, Q_R, Q_L, Q_A, W_R, W_L, W_A)
% Computes the power and energy input to the system

N=length(t);
t=reshape(t, [N, 1]);
pow=zeros(3,N);
for k=1:N
    pow(1,k) = tau(4:6,k)'*Q_R(:,:,k)*W_R(:,k); 
    pow(2,k) = tau(7:9,k)'*Q_L(:,:,k)*W_L(:,k); 
    pow(3,k) = tau(10:12,k)'*Q_A(:,:,k)*W_A(:,k); 
    E(k) = 0.5*m*x_dot(:,k)'*x_dot(:,k) - m*9.81*x(3,k);
end

E_dot = diff(E')./diff(t);
E_dot = [E_dot; E_dot(end)];

for k=1:N
    eff(k) = E_dot(k)/sum(pow(:,k));
end

end
