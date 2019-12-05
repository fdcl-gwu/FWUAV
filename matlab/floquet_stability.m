% load('sim_QS_x_hover_stability_data');

B = zeros(6, 6, 1+ix_d);
eig_vals_B = zeros(6, 1+ix_d);
for i=1:(1+ix_d)
    d_F = (F_linear(:, :, i) - F_linear(:, :, i+ix_d))./ F_linear(:, :, i);
    d_F(max(abs(F_linear(:, :, i+ix_d)), abs(F_linear(:, :, i))) < 1e-10) = 0;
    if(~all(abs(d_F) < 1e-2, [1, 2]))
        disp(d_F)
    end
    B(:, :, i) = delta_mat(:, :, i) \ delta_mat(:, :, i+ix_d);
    eig_vals_B(:, i) = sort(eig(B(:, :, i)));
    if(~all(abs(eig_vals_B(:, i) ./ eig_vals_B(:, 1) - 1) < 1e-2))
        disp(eig_vals_B(:, i))
        disp(i)
    end
%     if(rank(delta_mat(:, :, i), 1e-12) < 6)
%         disp(delta_mat(:, :, i))
%     end
end

B = B(:, :, round(end/2));
[e_vecs, rhos] = eig(B);
mus = log(diag(rhos)) * WK.f;

N_vars = size(F_linear, 2);
eig_vals = zeros(N_vars, N);

char_soln_mat = zeros(6, 6, N);
per_val_mat = zeros(6, 6, N);

for i=1:N
    eig_vals(:, i) = eig(F_linear(:, :, i));
    char_soln_mat(:, :, i) = delta_mat(:, :, i) * e_vecs;
    Y_0 = diag(exp(mus*t(i)));
    per_val_mat(:, :, i) = char_soln_mat(:, :, i) * inv(Y_0);
end

err = (trapz(t, sum(eig_vals)) - det(B)) / det(B);
eig_vals = sort(reshape(eig_vals(eig_vals ~= 0), N_vars/2, N), 1);

c_ix = 6;
per_val = vecnorm(squeeze(per_val_mat(:, c_ix, :)));
omega = 2*pi*WK.f;
ft_per=fittype('a*sin(omega*x + b) + c');
fit_val_per=fit(t(1:stop_idx)', per_val(1:stop_idx)', ft_per, ...
    'StartPoint', [1, 1, 1, 1], 'Lower', [-Inf, -Inf, -Inf, omega], 'Upper', [Inf, Inf, Inf, omega]);
fit_per_val = fit_val_per.a * sin(omega*t(1:stop_idx) + fit_val_per.b) + fit_val_per.c;

plot(t, per_val, 'k');
hold on;
plot(t, fit_per_val, 'b');

% per_fft = fft(per_val-fit_val_per.c);
% P2 = abs(per_fft/N);
% P1 = P2(1:round(N/2));
% P1(2:end-1) = 2*P1(2:end-1);
% f = WK.f*(1:round(N/2))/N;
% figure
% plot(f, P1);
% axis auto;

polarplot(complex(eig_vals_B(:, end)), 'r*');
% print('char_multipliers', '-depsc');