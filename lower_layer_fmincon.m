function T_star = lower_layer_fmincon(Sa, Q, h_nm, G_nn, G_nn_prime, sigma2, B1, f_cpu_sat, H_n, beta)
[N, M] = size(h_nm);
c = 3e8; T_trip = 2 * H_n * 1e3 / c;
q0 = ones(N, M) * (Q / (N * M));
obj = @(q) compute_T(q, Sa, N, M, h_nm, G_nn, G_nn_prime, sigma2, B1, f_cpu_sat, T_trip, beta);
options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
[~, T_star] = fmincon(obj, q0(:), ones(1,N*M), Q, [],[], zeros(N*M,1), [],[], options);
end

function T = compute_T(q_vec, Sa, N, M, h, G, Gp, sigma2, B1, f_cpu, T_trip, beta)
q = reshape(q_vec, N, M);
R = zeros(N,1);
for n = 1:N
    for m = 1:M
        interf = 0;
        for np = 1:N
            if np ~= n
                interf = interf + q(np,m) * Gp(n,m,np) * h(np,m);
            end
        end
        snr = q(n,m) * G(n,m) * h(n,m) / (interf + sigma2);
        R(n) = R(n) + B1 * log2(1 + max(snr,1e-15));
    end
end
sumR = sum(R);
if sumR < 1e-9, T = 1e10; return; end
T_s = Sa / sumR;
T_comp = max(R .* Sa * beta ./ (f_cpu * sumR) + T_trip);
T = T_s + T_comp;
end