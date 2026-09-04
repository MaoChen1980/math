# -*- coding: utf-8 -*-
"""工具 -> 新命题: 推导 + 数值检验.
P1: Mahler/2-kernel 生成函数 => v2 序列 DFT 闭式 (S-E2 升级)
P2: p-adic 赋值/LTE => L 链长 = 1 + v3(r+1) (S-H1 机制)
P3: Ramanujan 和/类和分解 => sigma 谱梯幅值闭式 (q=4 情形)
P4: LTE => 2^K-1 可证窗口延长 v2(3^K-1)=2+v2(K)
附: sigma 残差自相关衰减 (Wiener 连续谱判据支持)
"""
import math
import numpy as np

ok = lambda name, cond, detail="": print(
    f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")

memo = {1: 0}


def sg(n):
    path = []
    x = n
    while x not in memo:
        path.append(x)
        x = x >> 1 if x % 2 == 0 else 3 * x + 1
    v = memo[x]
    for y in reversed(path):
        v += 1
        memo[y] = v
    return memo[n]


# ============ P1: v2 序列 DFT 闭式 ============
# g(j)=v2(3j+2); 2-kernel: g(2m)=1+h(m), g(2m+1)=0; h(2m)=0, h(2m+1)=1+g(m)
# 展开: G[k] = sum_{i=1}^{K} 1[(N/2^i)|k]*(N/2^i)*omega^{k*e_i}
#        + omega^{k*e_base}
#   e_i = (2*4^{i-1}-2)/3 (分形指数 0,2,10,42,...), e_base = (2^{2K+1}-2)/3
N, K = 256, 8
g = np.array([(lambda m: (m & -m).bit_length() - 1)(3 * j + 2)
              for j in range(N)], float)
G_np = np.fft.fft(g)


def closed(k):
    s = 0j
    for i in range(1, K + 1):
        q = N // 2 ** i
        if k % q == 0:
            e = ((k % N) * (((2 * 4 ** (i - 1) - 2) // 3) % N)) % N
            s += q * np.exp(-2j * np.pi * e / N)
    e = ((k % N) * (((2 * 4 ** K - 2) // 3) % N)) % N
    s += np.exp(-2j * np.pi * e / N)
    return s


diff = max(abs(G_np[k] - closed(k)) for k in range(N))
ok("P1 v2-DFT 闭式 vs numpy.fft (N=256, 全 k)", diff < 1e-8,
   f"max_diff={diff:.2e}; G[N/2]=N? {G_np[N//2].real:.1f} vs {N}")

# ============ P2: L 链长 = 1 + v3(r+1) ============
# 恒等式: L(m)+1 = 2(m+1)/3 => v3 每链一步恰减 1
# 注: 测真实链长(不去重); M20 采样器的 seen 去重会截断链, 不属定理对象
nodes, seen = [], set()
roots_used = []
r = 2
while len(nodes) < 20000:
    m, ln = r, 1
    while (2 * m - 1) % 3 == 0:
        m = (2 * m - 1) // 3
        ln += 1
    roots_used.append((r, ln))
    # 采样器视角: 沿链加入未访问节点
    m2 = r
    while True:
        if m2 in seen:
            break
        seen.add(m2)
        nodes.append(m2)
        if (2 * m2 - 1) % 3 != 0:
            break
        m2 = (2 * m2 - 1) // 3
    r += 3


def v3(x):
    v = 0
    while x % 3 == 0:
        x //= 3
        v += 1
    return v


bad = [(r, ln) for r, ln in roots_used[:-1] if ln != 1 + v3(r + 1)]
ok("P2 真实L链长=1+v3(r+1) (全部根)", not bad,
   f"roots={len(roots_used)-1} 不符={len(bad)}")

# ============ P3': sigma 谱梯 q=4 分解 ============
# S0(4m)=2q+Σσ(m); S1(4m-3)=3q-3+Σσ(3m-2); S2(4m-2)=3q-3+Σσ(3m-1)
#   [S1/S2 的 m=1 项撞轨道终点 n=1, 恒等式例外, 各差 3]
# S3(4m-1)=不闭合 (奇链: 4m-1->6m-1->9m-1->...)
# X[N/4] = [s32 + 3q - 3 - S3] + i*[s_m - s31 - q + 3]
# 命题内容: 幅值 = 闭式部分 + 奇类残差 S3 (3-膨胀困难在每个谱梯阶级重现, 呼应 S-F3)
M = 20000
sig_arr = np.array([sg(n) for n in range(1, M + 1)], float)
q = M // 4
s_m = sum(sg(m) for m in range(1, q + 1))
s_31 = sum(sg(3 * m - 1) for m in range(1, q + 1))
s_32 = sum(sg(3 * m - 2) for m in range(1, q + 1))
S3 = sum(sg(4 * m - 1) for m in range(1, q + 1))
CF = (s_32 + 3 * q - 3 - S3) + 1j * (s_m - s_31 - q + 3)
X = np.fft.fft(sig_arr)[M // 4]
closed_share = (s_32 + 3 * q - 3) / (s_32 + 3 * q - 3 - S3)
ok("P3' sigma X[N/4] = 闭式部分 + 奇类残差 S3 (q=4 分解恒等式)",
   abs(X - CF) < 1e-6,
   f"dft={X:.1f} closed={CF:.1f}; S3 占实部 |S3|/|实部|="
   f"{abs(S3/(s_32 + 3*q - 3 - S3)):.3f}")

# ============ P4: LTE 可证窗口延长 ============
Kb = 10000
v2_3k = ((3 ** Kb - 1) & -(3 ** Kb - 1)).bit_length() - 1
pred = 2 + (Kb & -Kb).bit_length() - 1  # K 偶: v2(3^K-1)=2+v2(K)
ok("P4 LTE: v2(3^K-1)=2+v2(K)", v2_3k == pred, f"v2={v2_3k} pred={pred}")
# 直接迭代验证: 2K 步后到 3^K-1, 再 v2 次减半
x = 2 ** Kb - 1
for _ in range(2 * Kb):
    x = x >> 1 if x % 2 == 0 else 3 * x + 1
ok("P4 迭代 2K 步后 = 3^K-1", x == 3 ** Kb - 1)
for _ in range(v2_3k):
    x = x >> 1
ok("P4 再 v2 步减半 = (3^K-1)/2^v2", x == (3 ** Kb - 1) >> v2_3k)

# ============ 附: sigma 残差自相关 (长记忆发现) ============
t = np.arange(M)
sl, ic = np.polyfit(t, sig_arr, 1)
res = sig_arr - (sl * t + ic)
ks = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
rhos = [float(np.corrcoef(res[:-k], res[k:])[0, 1]) for k in ks]
# 发现性检验: 自相关不衰减 => 长程相关 => Wiener 快判据不适用
ok("附 残差长记忆: rho(k) 不衰减 (新发现, 谱测度类型开放)",
   abs(rhos[-1]) > 0.05 and abs(rhos[-1]) > 0.3 * abs(rhos[0]),
   f"rho(1)={rhos[0]:.4f} rho(64)={rhos[6]:.4f} rho(4096)={rhos[-1]:.4f}")
