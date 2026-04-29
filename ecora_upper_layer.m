function [E_total, S_opt, f_opt, p_opt] = ecora_upper_layer(Sa, T_star, alpha, g_jk, B0, sigma2, kappa, beta, T_tot, f_max_j, w_j, P_max)
%% E-CORA 上层优化（带 sum(S_j) = Sa 约束）
%  通过拉格朗日对偶法，找到满足 sum(S_j)=Sa 的最优卸载比特分配

J = length(alpha);
T_tx = T_tot - T_star;
if T_tx <= 0.05
    E_total = Inf; S_opt = []; f_opt = []; p_opt = [];
    return;
end
if nargin < 12, P_max = 0.2; end

% 各设备选最优信道
g_best = zeros(J,1);
for j = 1:J
    [~, k] = max(g_jk(j,:));
    g_best(j) = g_jk(j,k);
end

% Sa=0: 纯本地计算
if Sa <= 1e-12
    S_opt = zeros(J,1); f_opt = alpha * beta / T_tot; p_opt = zeros(J,1);
    for j = 1:J
        f_opt(j) = min(f_opt(j), f_max_j(j));
    end
    E_lc = kappa * (alpha - S_opt) .* beta .* (f_opt.^2);
    E_tx = 0;
    E_total = sum(w_j .* (E_lc + E_tx));
    return;
end

% 拉格朗日对偶：二分搜索 lambda 使得 sum(S_j) = Sa
% 对给定 lambda，每个设备最小化：E_j(s) + lambda * s
% lambda → -∞: s_j = alpha_j, sum = sum(alpha)
% lambda → +∞: s_j = 0, sum = 0

% 先估计 lambda 范围
dE0 = zeros(J,1);
for j = 1:J
    dE0(j) = -3*kappa*alpha(j)^2*beta^3/T_tot^2 + sigma2/g_best(j)*log(2)/B0;
end
lambda_low = min(dE0) - 10;   % 肯定让所有设备全卸载
lambda_high = max(-dE0) + 10; % 肯定让所有设备不卸载

% 边界检查
[sum_low, ~] = sum_s_lambda(lambda_low, alpha, g_best, T_tx, T_tot, sigma2, kappa, beta, B0, T_tot, f_max_j, w_j, P_max);
[sum_high, ~] = sum_s_lambda(lambda_high, alpha, g_best, T_tx, T_tot, sigma2, kappa, beta, B0, T_tot, f_max_j, w_j, P_max);

if sum_low < Sa
    % 即使最极端情况也达不到 Sa，取 sum 最大的解
    lambda = lambda_low;
elseif sum_high > Sa
    % 最小卸载量都大于 Sa，取 sum 最小的解
    lambda = lambda_high;
else
    % 二分搜索
    for iter = 1:60
        lambda = (lambda_low + lambda_high) / 2;
        [total_s, ~] = sum_s_lambda(lambda, alpha, g_best, T_tx, T_tot, sigma2, kappa, beta, B0, T_tot, f_max_j, w_j, P_max);
        if abs(total_s - Sa) < 1e-6
            break;
        end
        if total_s > Sa
            lambda_low = lambda;  % 需要更大的 lambda 来减少卸载
        else
            lambda_high = lambda; % 需要更小的 lambda 来增加卸载
        end
    end
end

% 用找到的 lambda 计算最终结果
[~, S_opt, f_opt, p_opt] = sum_s_lambda(lambda, alpha, g_best, T_tx, T_tot, sigma2, kappa, beta, B0, T_tot, f_max_j, w_j, P_max);

E_lc = kappa * (alpha - S_opt) .* beta .* (f_opt.^2);
E_tx = p_opt * T_tx;
E_total = sum(w_j .* (E_lc + E_tx));
end

function [total_s, S_opt, f_opt, p_opt] = sum_s_lambda(lambda, alpha, g_best, T_tx, T_tot, sigma2, kappa, beta, B0, ~, f_max_j, w_j, P_max)
    J = length(alpha);
    S_opt = zeros(J,1);
    f_opt = zeros(J,1);
    p_opt = zeros(J,1);

    for j = 1:J
        alpha_j = alpha(j);
        g = g_best(j);
        w = w_j(j);

        s_min = max(0, alpha_j - T_tot * f_max_j(j) / beta);
        max_rate = B0 * log2(1 + P_max * g / sigma2);
        s_max = min(alpha_j, T_tx * max_rate);
        if s_max < s_min, s_max = s_min; end

        % 带拉格朗日项的目标函数
        obj = @(s) w * (kappa*(alpha_j - s)^3*beta^3/T_tot^2 ...
                     + T_tx*sigma2/g*(2^(s/(B0*T_tx))-1) + lambda*s);
        s_opt = fminbnd(obj, s_min, s_max);
        S_opt(j) = s_opt;

        if s_opt > 0
            p_opt(j) = sigma2/g * (2^(s_opt/(B0*T_tx)) - 1);
            if p_opt(j) > P_max, p_opt(j) = P_max; end
        end
        f_opt(j) = (alpha_j - s_opt) * beta / T_tot;
        f_opt(j) = min(f_opt(j), f_max_j(j));
    end

    total_s = sum(S_opt);
end
