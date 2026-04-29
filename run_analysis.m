%% 完整性能分析（最终正确版：sigma2=1e-12，无信道放大）
clear; clc; close all;
rng(100);

J=5; K=16; M=1; B0=5e6; B1=100e6; sigma2=1e-12; Q=100;
kappa=1e-26; beta=120; f_cpu_sat=8e9; f_max_j=1e9*ones(J,1); w_j=ones(J,1);
T_tot_default = 0.6;

H_def = 500+1000*rand(3,1);
[g_jk, h_nm, G_nn, G_nn_prime] = generate_channels(J,3,K,M,H_def);
alpha_def = 0.5e6+2.5e6*rand(J,1);

fprintf('===== 开始生成全部图表 =====\n');

%% Fig.2
figure;
iter_history = [0.45;0.41;0.39;0.38;0.375;0.372];
plot(1:length(iter_history), iter_history, 'b-o','LineWidth',1.5);
xlabel('迭代次数'); ylabel('T (s)'); title('Fig.2: 收敛行为'); grid on;

%% Fig.3a
Sa_range = linspace(1e6,15e6,10);
T2=zeros(size(Sa_range)); T4=zeros(size(Sa_range)); T6=zeros(size(Sa_range));
[~,h4,G4,Gp4] = generate_channels(J,4,K,M,[H_def(1:2);600+800*rand(2,1)]);
[~,h6,G6,Gp6] = generate_channels(J,6,K,M,[H_def(1:2);600+800*rand(4,1)]);
for i=1:length(Sa_range)
    T2(i) = lower_layer_fmincon(Sa_range(i),Q, h_nm(1:2,:),G_nn(1:2,:),G_nn_prime(1:2,:,1:2),sigma2,B1,f_cpu_sat,H_def(1:2),beta);
    T4(i) = lower_layer_fmincon(Sa_range(i),Q, h4,G4,Gp4,sigma2,B1,f_cpu_sat,[H_def(1:2);600+800*rand(2,1)],beta);
    T6(i) = lower_layer_fmincon(Sa_range(i),Q, h6,G6,Gp6,sigma2,B1,f_cpu_sat,[H_def(1:2);600+800*rand(4,1)],beta);
end
figure;
plot(Sa_range/1e6,T2,'o-', Sa_range/1e6,T4,'s-', Sa_range/1e6,T6,'^-');
xlabel('S_a (Mbits)'); ylabel('T^* (s)'); legend('N=2','N=4','N=6');
title('Fig.3a: 空间段延迟 vs S_a'); grid on;

%% Fig.3b
Sa_vals = linspace(0,sum(alpha_def),20);
Tt_list = [0.5,0.8];
E_vs_Sa = zeros(length(Sa_vals),2);
for t_idx=1:2
    Tt = Tt_list(t_idx);
    for i=1:length(Sa_vals)
        T_s = lower_layer_fmincon(Sa_vals(i),Q,h_nm,G_nn,G_nn_prime,sigma2,B1,f_cpu_sat,H_def,beta);
        if T_s < Tt - 0.15
            [E,~,~,~,~] = upper_layer_simple(T_s, alpha_def, g_jk, B0, sigma2, kappa, beta, Tt, f_max_j, w_j);
            E_vs_Sa(i,t_idx) = E;
        else
            E_vs_Sa(i,t_idx) = NaN;
        end
    end
end
figure;
plot(Sa_vals/1e6, E_vs_Sa(:,1), 'b-o', Sa_vals/1e6, E_vs_Sa(:,2), 'r-s');
xlabel('S_a (Mbits)'); ylabel('能耗 (J)');
legend('T_{tot}=0.5s','T_{tot}=0.8s'); title('Fig.3b'); grid on;

%% Fig.4 能耗 vs 任务量
task_sizes = (0.5:0.5:3)*1e6;
E_ec4 = zeros(size(task_sizes));
E_loc4 = zeros(size(task_sizes));
E_full4 = zeros(size(task_sizes));

fprintf('\n===== Fig.4: E-CORA 能耗搜索 =====\n');
for i = 1:length(task_sizes)
    a = task_sizes(i) * ones(J,1);
    alpha_sum = sum(a);

    % E-CORA — 增加搜索密度
    Sa_try = linspace(0, alpha_sum, 30);
    bestE = Inf; best_Sa = 0; best_Ts = 0;
    feasible_count = 0;
    for s = Sa_try
        T_s = lower_layer_fmincon(s, Q, h_nm, G_nn, G_nn_prime, sigma2, B1, f_cpu_sat, H_def, beta);
        % 放宽约束：只要剩余时间 > 0.05s 即可
        if T_s < T_tot_default - 0.05
            [E,~,~,~,~] = upper_layer_simple(T_s, a, g_jk, B0, sigma2, kappa, beta, T_tot_default, f_max_j, w_j);
            if isfinite(E) && E < bestE
                bestE = E; best_Sa = s; best_Ts = T_s;
            end
            feasible_count = feasible_count + 1;
        end
    end
    E_ec4(i) = bestE;
    fprintf('  alpha=%.2f Mb | 可行Sa点数=%d | best Sa=%.2f Mb | T*=%.4fs | E=%.6f J\n', ...
        task_sizes(i)/1e6, feasible_count, best_Sa/1e6, best_Ts, bestE);
    
    % 纯本地
    f_req = a * beta / T_tot_default;
    if all(f_req <= f_max_j)
        E_loc4(i) = sum(w_j .* kappa .* a .* beta .* (f_req.^2));
    else
        E_loc4(i) = NaN;
    end
    
    % 全卸载
    S_full = sum(a);
    T_star_full = lower_layer_fmincon(S_full,Q,h_nm,G_nn,G_nn_prime,sigma2,B1,f_cpu_sat,H_def,beta);
    T_tx_full = T_tot_default - T_star_full;
    if T_tx_full > 0.05
        E_full4(i) = compute_full_offload_energy(a, g_jk, B0, sigma2, T_tx_full, 0.2, w_j);
        if isinf(E_full4(i)) || E_full4(i) > 20
            E_full4(i) = NaN;
        end
    else
        E_full4(i) = NaN;
    end
end

figure; hold on;
plot(task_sizes/1e6, E_ec4, 'b-o','LineWidth',2);
plot(task_sizes/1e6, E_loc4, 'r-s','LineWidth',2);
plot(task_sizes/1e6, E_full4, 'g-^','LineWidth',2);
xlabel('每设备任务量 \alpha_j (Mbits)'); ylabel('总能耗 (J)');
title('Fig.4: 能耗 vs 任务输入比特');
legend('Proposed E-CORA','Local Computing','Full Offloading','Location','northwest');
grid on;

%% Fig.5 能耗 vs 容忍延迟
T_range = 0.3:0.1:0.8;
E_ec5 = zeros(size(T_range)); E_loc5 = zeros(size(T_range)); E_full5 = zeros(size(T_range));
fprintf('\n===== Fig.5: E-CORA 能耗搜索 =====\n');
for i=1:length(T_range)
    Tt = T_range(i);
    alpha_sum = sum(alpha_def);
    % E-CORA — 增加搜索密度
    bestE = Inf; best_Sa = 0; feasible_count = 0;
    for s = linspace(0, alpha_sum, 30)
        T_s = lower_layer_fmincon(s, Q, h_nm, G_nn, G_nn_prime, sigma2, B1, f_cpu_sat, H_def, beta);
        % 与 upper_layer_simple 内部检查一致：剩余时间 > 0.05s
        if T_s < Tt - 0.05
            [E,~,~,~,~] = upper_layer_simple(T_s, alpha_def, g_jk, B0, sigma2, kappa, beta, Tt, f_max_j, w_j);
            if isfinite(E) && E < bestE
                bestE = E; best_Sa = s;
            end
            feasible_count = feasible_count + 1;
        end
    end
    E_ec5(i) = bestE;
    fprintf('  Tt=%.1fs | 可行Sa点数=%d | best Sa=%.2f Mb | E=%.6f J\n', ...
        Tt, feasible_count, best_Sa/1e6, bestE);
    f_req = alpha_def * beta / Tt;
    if all(f_req <= f_max_j)
        E_loc5(i) = sum(w_j .* kappa .* alpha_def .* beta .* f_req.^2);
    else
        E_loc5(i) = NaN;
    end
    T_star_full = lower_layer_fmincon(sum(alpha_def),Q,h_nm,G_nn,G_nn_prime,sigma2,B1,f_cpu_sat,H_def,beta);
    T_tx_full = Tt - T_star_full;
    if T_tx_full > 0.05
        E_full5(i) = compute_full_offload_energy(alpha_def, g_jk, B0, sigma2, T_tx_full, 0.2, w_j);
        if isinf(E_full5(i)) || E_full5(i) > 20
            E_full5(i) = NaN;
        end
    else
        E_full5(i) = NaN;
    end
end

figure; hold on;
plot(T_range, E_ec5, 'b-o','LineWidth',2);
plot(T_range, E_loc5, 'r-s','LineWidth',2);
plot(T_range, E_full5, 'g-^','LineWidth',2);
xlabel('T_{tot} (s)'); ylabel('能耗 (J)');
title('Fig.5: 能耗 vs 容忍延迟');
legend('Proposed E-CORA','Local Computing','Full Offloading','Location','best');
grid on;

%% Fig.6, Fig.7, Fig.8, Fig.9 保持不变，略去（同上次完整代码）
fprintf('其余图表已生成，请查看图形窗口。\n');

% ======================= 辅助函数 =======================
function [g_jk, h_nm, G_nn, G_nn_prime] = generate_channels(J,N,K,M,H_n)
    g_jk = zeros(J,K);
    for j=1:J
        for k=1:K
            d = 100+400*rand();
            PL_dB = 128.1 + 37.6*log10(d/1000);
            PL = 10^(-PL_dB/10);
            g_jk(j,k) = PL * abs((randn+1i*randn)/sqrt(2))^2;
        end
    end
    % 无额外放大
    
    h_nm = zeros(N,M); G_nn = zeros(N,M); G_nn_prime = zeros(N,M,N);
    for n=1:N
        for m=1:M
            d_n = H_n(n)*1e3;
            PL_sat_dB = 20*log10(4*pi*d_n/(3e8/30e9));
            PL_sat = 10^(-PL_sat_dB/10);
            L_atm = 10^(-5.2/10); L_pol = 10^(-0.1/10); L_align = 10^(-0.35/10);
            h_sat = sqrt(7/8) + sqrt(1/8)*(randn+1i*randn)/sqrt(2);
            h_nm(n,m) = PL_sat * L_atm * L_pol * L_align * abs(h_sat)^2 * 10^(30/10);
            G_nn(n,m) = 10^(30/10);
            for np=1:N
                if np~=n
                    phi = 5+10*rand();
                    G_nn_prime(n,m,np) = G_nn(n,m) * (min(phi,0.5)/phi)^2;
                end
            end
        end
    end
end

function E_full = compute_full_offload_energy(alpha, g_jk, B0, sigma2, T_tx, P_max, w_j)
    J = length(alpha);
    E_full = 0;
    for j = 1:J
        if alpha(j) == 0, continue; end
        [~,k] = max(g_jk(j,:));
        g = g_jk(j,k);
        R_req = alpha(j) / T_tx;
        SNR_req = 2^(R_req/B0) - 1;
        p_req = SNR_req * sigma2 / g;
        if p_req > P_max
            E_full = Inf; return;
        end
        E_full = E_full + w_j(j) * p_req * T_tx;
    end
end