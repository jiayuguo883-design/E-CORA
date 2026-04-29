function [E_total, S_j_opt, f_opt, p_opt, x_opt] = upper_layer_simple(T_star, alpha, g_jk, B0, sigma2, kappa, beta, T_tot, f_max_j, w_j)
J = length(alpha);
K = size(g_jk,2);
T_tx = T_tot - T_star;
if T_tx <= 0.05
    E_total = Inf; S_j_opt = []; f_opt = []; p_opt = []; x_opt = [];
    return;
end

S_j_opt = zeros(J,1);
f_opt = zeros(J,1);
p_opt = zeros(J,1);
x_opt = zeros(J,K);
P_max = 0.2;   % 设备最大功率

for j = 1:J
    [~, k] = max(g_jk(j,:));
    g = g_jk(j,k);
    x_opt(j,k) = 1;
    
    alpha_j = alpha(j);
    s_min = max(0, alpha_j - T_tot * f_max_j(j) / beta);
    max_rate = B0 * log2(1 + P_max * g / sigma2);
    s_max = min(alpha_j, T_tx * max_rate);
    if s_max < s_min
        s_max = s_min;
    end
    
    obj = @(s) w_j(j) * ( kappa * (alpha_j - s)^3 * beta^3 / T_tot^2 ...
                + T_tx * sigma2/g * (2^(s/(B0*T_tx)) - 1) );
    s_opt = fminbnd(obj, s_min, s_max);
    S_j_opt(j) = s_opt;
    
    if s_opt > 0
        p_opt(j) = sigma2/g * (2^(s_opt/(B0*T_tx)) - 1);
    else
        p_opt(j) = 0;
    end
    f_opt(j) = (alpha_j - s_opt) * beta / T_tot;
    f_opt(j) = min(f_opt(j), f_max_j(j));
end

E_lc = kappa * (alpha - S_j_opt) .* beta .* (f_opt.^2);
E_tx = p_opt * T_tx;
E_total = sum(w_j .* (E_lc + E_tx));
end