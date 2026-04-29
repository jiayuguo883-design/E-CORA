# E-CORA 项目文档

## 项目概述

E-CORA（Energy-efficient Computation Offloading and Resource Allocation）是发表于 IEEE Internet of Things Journal 的论文的 MATLAB 仿真实现。

**论文信息：**
- 标题：Energy-Efficient Multiaccess Edge Computing for Terrestrial-Satellite Internet of Things
- 作者：Zhengyu Song, Yuanyuan Hao, Yuanwei Liu, Xin Sun
- 期刊：IEEE Internet of Things Journal, Vol. 8, No. 18, September 15, 2021
- DOI：10.1109/JIOT.2021.3068141

## 系统模型

### 网络架构

- **IMD（IoT Mobile Device）**：J = 5 个物联网移动设备，每个有 alpha_j 比特任务输入
- **TST（Terrestrial-Satellite Terminal）**：1 个地面卫星终端，作为接入点，无计算能力
- **LEO 卫星**：N = 3 颗低轨卫星，部署 MEC 服务器，提供边缘计算能力
- **频段**：C 波段（IMD→TST，地面段），Ka 波段（TST→卫星，空间段）

### 任务卸载流程（两阶段）

1. **地面段**：IMD 通过 OFDMA 将任务比特卸载到 TST
   - 占用时间 T_tx，带宽 B0 = 5 MHz，K = 16 个子载波
   - S_j = 设备 j 卸载到 TST 的比特数
   - Sa = sum(S_j) = 总卸载比特数
2. **空间段**：TST 将收集的比特通过 Ka 波段上传到多颗 LEO 卫星
   - 占用时间 T_star = T(Sa)，带宽 B1 = 100 MHz，M = 1 个子载波
   - 卫星 MEC 服务器并行执行边缘计算

### 本地计算

未卸载的比特（alpha_j - S_j）在 IMD 本地执行，利用 DVFS 调频。

## 双层优化框架

### 下层问题 —— 最小化空间段延迟（lower_layer_fmincon.m）

```
min_q  T(Sa) = Ts + max_n (TnC + Tntrip)
s.t.   sum(q) <= Q
```

- Ts = Sa / sum(RnT) —— TST 到卫星的传输时间
- TnC = RnT \* Ts \* beta / fnC —— 卫星 n 的边缘计算时间
- Tntrip = 2 \* Hn / c —— 往返传播延迟
- q = TST 发射功率分配，通过 fmincon（interior-point 算法）求解
- **优化变量**：qnm（TST 到卫星 n 在子载波 m 上的发射功率）
- **约束**：总功率不超过 Q = 100

### 上层问题 —— 最小化设备加权和能耗（ecora_upper_layer.m）

```
min_{f,s,p}  sum(w_j * (E_j_LC + E_j_TX))
s.t.         sum(S_j) = Sa
             T_tx + T(Sa) = T_tot
             0 <= S_j <= alpha_j
```

- E_j_LC = kappa \* (alpha_j - S_j) \* beta \* f_j^2 —— 本地计算能耗
- E_j_TX = p_j \* T_tx —— 卸载传输能耗
- **求解方法**：拉格朗日对偶法（Lagrangian dual decomposition）
- 通过二分搜索拉格朗日乘子 lambda，使 sum(S_j(lambda)) = Sa
- 每个设备用 fminbnd 在 [s_min, s_max] 内求解

### 线性搜索

由于 T(Sa) 无法解析表达，在 [0, sum(alpha)] 范围内对 Sa 进行线性搜索（30 个点），对每个 Sa 求解上下层问题，取使能耗最小的 Sa。

## 文件结构

| 文件 | 功能 |
|------|------|
| **E_CORA_main.m** | 主程序。生成信道，线性搜索 Sa，绘制 E vs Sa 曲线 |
| **lower_layer_fmincon.m** | 下层优化。用 fmincon 求解最小空间段延迟 T(Sa)。内含 compute_T 辅助函数 |
| **ecora_upper_layer.m** | 上层优化[新]。带 sum(S_j)=Sa 约束的拉格朗日对偶法。内含 sum_s_lambda 辅助函数 |
| **upper_layer_simple.m** | 上层优化[旧]。独立优化各设备 S_j（无 sum 约束），已不被主脚本调用，仅作参考 |
| **run_analysis.m** | 完整性能分析。生成 Fig.2~Fig.9 对比图。内含 generate_channels 和 compute_full_offload_energy 辅助函数 |
| **E_CORA_results.mat** | 保存的仿真结果 |
| **AGENTS.md** | 本文档 |

## 信道模型

### 地面信道（C 波段）

- 大尺度衰落：标准路径损耗 PL_dB = 128.1 + 37.6\*log10(d/1000)
- 小尺度衰落：瑞利衰落
- 距离 d 在 [100, 500]m 间随机

### 卫星信道（Ka 波段）

- 自由空间路径损耗：PL_sat_dB = 20\*log10(4\*pi\*Hn / lambda)
- 大气损耗：5.2 dB
- 极化损耗：0.1 dB
- 指向损耗：0.35 dB
- 小尺度衰落：Shadowed-Rician 衰落（K 因子 = 7）
- 卫星高度 Hn 在 [500, 1500]km 间随机
- 频率 f = 30 GHz

## 关键参数

| 参数 | 值 | 说明 |
|------|-----|------|
| J | 5 | 设备数 |
| N | 3 | 卫星数 |
| K | 16 | 地面子载波数 |
| M | 1 | 卫星子载波数 |
| B0 | 5 MHz | C 波段带宽 |
| B1 | 100 MHz | Ka 波段带宽 |
| sigma2 | 1e-12 | 噪声功率 |
| Q | 100 | TST 最大发射功率 |
| kappa | 1e-26 | 芯片架构系数 |
| beta | 120 cycles/bit | 计算强度 |
| f_cpu_sat | 8e9 cycles/s | 卫星 CPU 频率 |
| f_max_j | 1e9 cycles/s | 设备最大 CPU 频率 |
| P_max | 0.2 W | 设备最大发射功率 |
| T_tot | 0.6 s | 任务容忍延迟（默认） |

## 仿真输出图表

| 图 | 内容 |
|----|------|
| **Fig.2** | SFP 算法收敛行为 |
| **Fig.3a** | 空间段最小延迟 vs Sa（N=2,4,6 卫星对比） |
| **Fig.3b** | 能耗 vs Sa（Ttot=0.5s 和 0.8s），U 形曲线 |
| **Fig.4** | 能耗 vs 任务输入比特（E-CORA vs Local vs Full Offloading） |
| **Fig.5** | 能耗 vs 容忍延迟（E-CORA vs Local vs Full Offloading） |
| **Fig.6** | 能耗和卸载比特数 vs TST 最大发射功率 Q |
| **Fig.7** | 卸载比特在卫星间的分布（四种情形） |
| **Fig.8** | 能耗 vs 卫星数量 |
| **Fig.9** | 卸载比例 vs 卫星数量 |

## 论文预期结果

1. **Fig.3b**：能耗 vs Sa 呈 U 形，存在最优 Sa（Ttot=0.6s 时 Sa\* ≈ 8.048 Mbits）
2. **Fig.4**：E-CORA < Local Computing < Full Offloading
3. **Fig.5**：E-CORA < Local Computing < Full Offloading，Ttot 越小 E-CORA 优势越明显

## 状态记录

- 2026-04-29：修复了上层优化中缺少 sum(S_j)=Sa 约束的关键 bug
- 2026-04-29：新增 ecora_upper_layer.m 使用拉格朗日对偶法
- 2026-04-29：将 Sa 搜索点数从 10 增加到 30
- 2026-04-29：约束条件从 T_s < Ttot-0.15 放宽到 T_s < Ttot-0.05
- 2026-04-29：项目上传至 GitHub（jiayuguo883-design/E-CORA）
