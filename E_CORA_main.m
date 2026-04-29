%% E-CORA 主程序（最终正确版：sigma2=1e-12）
clear; clc; close all;

J = 5; N = 3; K = 16; M = 1;
B0 = 5e6; B1 = 100e6;
sigma2 = 1e-12;          % 关键：噪声功率与论文仿真一致
Q = 100;
rng(42);
alpha = 0.5e6 + 2.5e6 * rand(J,1);
T_tot = 0.6;
kappa = 1e-26;           % 使用原文值
beta = 120;
f_cpu_sat = 8e9;
f_max_j = 1e9 * ones(J,1);
w_j = ones(J,1);
H_n = 500 + 1000 * rand(N,1);

fprintf('===== 参数初始化 =====\n');
fprintf('噪声功率 sigma^2 = %.1e\n', sigma2);

% ---- 地面信道（标准路径损耗）----
g_jk = zeros(J,K);
for j = 1:J
    for k = 1:K
        d = 100 + 400*rand();
        PL_dB = 128.1 + 37.6*log10(d/1000);
        PL = 10^(-PL_dB/10);
        g_jk(j,k) = PL * abs((randn+1i*randn)/sqrt(2))^2;
    end
end

% ---- 卫星信道（同前）----
h_nm = zeros(N,M); G_nn = zeros(N,M); G_nn_prime = zeros(N,M,N);
for n = 1:N
    for m = 1:M
        d_n = H_n(n)*1e3;
        PL_sat_dB = 20*log10(4*pi*d_n/(3e8/30e9));
        PL_sat = 10^(-PL_sat_dB/10);
        L_atm = 10^(-5.2/10); L_pol = 10^(-0.1/10); L_align = 10^(-0.35/10);
        h_sat = sqrt(7/8) + sqrt(1/8)*(randn+1i*randn)/sqrt(2);
        h_nm(n,m) = PL_sat * L_atm * L_pol * L_align * abs(h_sat)^2 * 10^(30/10);
        G_nn(n,m) = 10^(30/10);
        for np = 1:N
            if np ~= n
                phi = 5 + 10*rand();
                G_nn_prime(n,m,np) = G_nn(n,m) * (min(phi,0.5)/phi)^2;
            end
        end
    end
end

% ---- 线性搜索 S_a ----
S_a_values = linspace(0, sum(alpha), 30);
E_vals = inf(size(S_a_values));
T_star_vals = zeros(size(S_a_values));
best_E = inf; best_Sa = 0;

fprintf('\n===== 开始搜索 S_a =====\n');
for idx = 1:length(S_a_values)
    Sa = S_a_values(idx);
    T_star = lower_layer_fmincon(Sa, Q, h_nm, G_nn, G_nn_prime, sigma2, B1, f_cpu_sat, H_n, beta);
    if T_star >= T_tot - 0.15
        continue;
    end
    [E, ~, ~, ~, ~] = upper_layer_simple(T_star, alpha, g_jk, B0, sigma2, kappa, beta, T_tot, f_max_j, w_j);
    if ~isfinite(E) || E <= 0
        continue;
    end
    E_vals(idx) = E;
    T_star_vals(idx) = T_star;
    fprintf('Sa=%.2f Mb, T*=%.3fs, Ttx=%.3fs, E=%.4f J\n', Sa/1e6, T_star, T_tot-T_star, E);
    if E < best_E
        best_E = E; best_Sa = Sa;
    end
end

valid = find(isfinite(E_vals) & E_vals>0);
if isempty(valid)
    error('无有效数据，请检查卫星链路参数');
end

fprintf('\n最优 Sa* = %.2f Mbits, 最小能耗 = %.4f J\n', best_Sa/1e6, best_E);

figure;
plot(S_a_values(valid)/1e6, E_vals(valid), 'b-','LineWidth',2); hold on;
plot(best_Sa/1e6, best_E, 'ro','MarkerSize',8,'MarkerFaceColor','r');
xlabel('S_a (Mbits)'); ylabel('总能耗 (J)');
title('E-CORA: 能耗 vs S_a'); grid on;

save('E_CORA_results.mat', 'best_Sa', 'best_E', 'alpha', 'J', 'N');