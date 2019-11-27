load('sim_QS_x_hover_stability_data');

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

polarplot(complex(eig_vals_B(:, end)), 'r*');
% print('char_multipliers', '-depsc');