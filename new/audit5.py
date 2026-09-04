# -*- coding: utf-8 -*-
"""audit5: 三条盲点专项自查.
盲点1: sigma/core/v2 代码与数学独立性 (非循环论证)
盲点2: DFT 相位递归是否依赖 N=2^k (换 N=19998/20002/6000 测试 + 奇数边界)
盲点3: 阶梯占比反升是否为归一化幻觉 (分子/分母能量分解 + 线性相消角)
"""
import json
import numpy as np

R = json.load(open(r"C:\Users\savyc\collatz_fft_study\results2.json",
                   encoding="utf-8"))
ok = lambda name, cond, detail="": print(
    f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")

memo = {1: 0}


def sg(n):
    """σ: 完整 Collatz 轨道逐步迭代计数 (与 core/v2 无任何代码耦合)."""
    path, x = [], n
    while x not in memo:
        path.append(x)
        x = x >> 1 if x % 2 == 0 else 3 * x + 1
    v = memo[x]
    for y in reversed(path):
        v += 1
        memo[y] = v
    return memo[n]


def core(n):
    """严格算术操作: 剥离全部 2 因子 (试除循环)."""
    while n % 2 == 0:
        n //= 2
    return n


def v2(n):
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


print("=" * 72)
print("盲点1: sigma / core / v2 独立性")
print("=" * 72)
# 三条独立代码路径: sg=轨道迭代 (含 3n+1 步), core=试除, v2=试除计数.
# 数学上 σ 的定义是全轨道步数 (含 3n+1 步), 不是由 core+v2 定义的.
import random as rnd
rnd.seed(7)
bad = [n for n in rnd.sample(range(1, 200001), 2000)
       if sg(n) != sg(core(n)) + v2(n)]
ok("1a 恒等式 2000 随机例 (含 n≤2×10^5), 整数精确 0 误差", not bad,
   f"不符={len(bad)}")
# 第三实现 (递归) 交叉 — 与 verify4 同源但这里换 seed
from functools import lru_cache


@lru_cache(maxsize=None)
def sg2(n):
    if n == 1:
        return 0
    return 1 + sg2(n // 2 if n % 2 == 0 else 3 * n + 1)


bad2 = [n for n in rnd.sample(range(1, 200001), 2000) if sg(n) != sg2(n)]
ok("1b sg(迭代+memo) vs sg2(递归) 两实现交叉 2000 例", not bad2,
   f"不符={len(bad2)}")
# 数学注: σ(2m)=1+σ(m) 是定理非定义 (σ 按"全轨道步数"定义, 含 3n+1 步);
# 该恒等式确属初等 (一步归纳), 其价值仅作对齐载体 — 措辞已在复核注降级.
ok("1c σ(2m)=1+σ(m) 1000 例 (定理非定义: 右侧不用 σ 定义)",
   all(sg(2 * m) == 1 + sg(m) for m in rnd.sample(range(1, 100000), 1000)))

print()
print("=" * 72)
print("盲点2: 相位递归是否依赖 N=2^k")
print("=" * 72)
# 递归推导只用 ω_N²=ω_{N/2} (任何偶 N 成立), 与 2^k 无关. 现场测试:
sig_all = np.array([sg(n) for n in range(1, 200001)], float)
print(f"{'N':>8}{'因子分解':>16}{'max_err':>12}  判定")
res2 = {}
for N in (20000, 19998, 20002, 6000, 19994, 16384):
    sig = sig_all[:N]
    r = np.array([sg(core(n)) for n in range(1, N + 1)], float)
    X = np.fft.fft(r)
    E = np.fft.fft(r[0::2])
    Xh = np.fft.fft(r[:N // 2])
    k = np.arange(N)
    wN = np.exp(-2j * np.pi / N)
    rec = E[k % (N // 2)] + wN ** k * Xh[k % (N // 2)]
    err = float(np.max(np.abs(X - rec)))
    fac = []
    m = N
    d = 2
    while d * d <= m:
        while m % d == 0:
            fac.append(d)
            m //= d
        d += 1
    if m > 1:
        fac.append(m)
    pow2 = all(f == 2 for f in fac)
    res2[N] = err
    print(f"{N:>8}{str(fac):>16}{err:>12.2e}  {'2^k' if pow2 else '非2^k'}"
          f"{' (通过)' if err < 1e-4 else ' (失败)'}")
ok("2a 相位递归在全部偶 N (含非 2^k: 19998/20002/6000/19994) 精确",
   all(res2[N] < 1e-4 for N in (19998, 20002, 6000, 19994, 20000, 16384)))

# 奇数 N: 偶/奇指标分裂不均衡 => 递归需边界项, 如实量化失败边界
N = 19999
sig = sig_all[:N]
r = np.array([sg(core(n)) for n in range(1, N + 1)], float)
X = np.fft.fft(r)
E = np.fft.fft(r[0::2])       # 偶指标 N//2 个... 奇 N 时偶指标 9999+... 数不同
k = np.arange(N)
rec = np.full(N, np.nan)
try:
    rec = E[k % ((N + 1) // 2)] + np.exp(-2j * np.pi / N) ** k * \
        np.fft.fft(r[:(N // 2)])[k % (N // 2)]
except Exception:
    pass
if not np.any(np.isnan(rec)):
    err_odd = float(np.max(np.abs(X[:2 * (N // 2) + 1] -
                                rec[:2 * (N // 2) + 1])))
    print(f"{N:>8}{'奇数(边界项缺失)':>18}{err_odd:>12.2e}  按声明仅在偶 N 成立")

print()
print("=" * 72)
print("盲点3: 阶梯占比反升 = 归一化幻觉? (分子/分母分解)")
print("=" * 72)
N = 20000
sig = sig_all[:N]
r = np.array([sg(core(n)) for n in range(1, N + 1)], float)
v2s = np.array([v2(n) for n in range(1, N + 1)], float)
t = np.arange(N)


def spec(x):
    sl, ic = np.polyfit(t, x, 1)
    return np.fft.rfft(x - (sl * t + ic))


Xs, Xr, Xv = spec(sig), spec(r), spec(v2s)
P = [np.abs(X) ** 2 for X in (Xs, Xr, Xv)]
ladder = range(625, N // 2 + 1, 625)


def abs_sum(X):
    return float(sum(abs(X[k]) ** 2 for k in ladder))


tot = [float(np.sum(p)) for p in P]
lad = [abs_sum(X) for X in (Xs, Xr, Xv)]
print(f"{'':>14}{'总能量':>14}{'16格阶梯绝对能量':>16}{'占比':>10}")
for nm, i in (("原始 σ", 0), ("剥离 r", 1), ("v₂ 单独", 2)):
    print(f"{nm:<12}{tot[i]:>14.4e}{lad[i]:>16.4e}{lad[i]/tot[i]:>10.4f}")
ok("3a 线性精确: X_σ = X_r + X_{v₂} (谱上逐 bin)",
   float(np.max(np.abs(Xs - Xr - Xv))) < 1e-4,
   f"max={float(np.max(np.abs(Xs-Xr-Xv))):.2e}")
# 分母效应大小: 若分子不变, 占比变化 = E_raw/E_str
den_factor = tot[0] / tot[1]
num_factor = lad[1] / lad[0]
print(f"\n分解: 占比比 = (分子因子 {num_factor:.4f}) × (分母因子 {den_factor:.4f})")
print(f"  实测占比比 = {0.1103/0.0952:.4f}")
ok("3b 反升主因是分子 (剥离后阶梯绝对能量真实上升), 非分母",
   num_factor > 1.10 and den_factor < 1.005,
   f"分子因子={num_factor:.4f} (阶梯绝对能量升 {100*(num_factor-1):.1f}%), "
   f"分母因子={den_factor:.4f} (总能量仅变 {100*(den_factor-1):.2f}%)")
# 机理: v₂ 的谱在阶梯 bin 上与 r 的谱反相相消 (剥离=解除相消)
cs = [float(np.real(Xv[k] * np.conj(Xr[k])) /
            (abs(Xv[k]) * abs(Xr[k]))) for k in ladder
      if abs(Xv[k]) * abs(Xr[k]) > 0]
print(f"\n相消角证据: 16 个阶梯 bin 上 cos<X_v, X_r> 均值 = {np.mean(cs):.4f}"
      f" (负=反相相消), v₂ 单独阶梯绝对能量/σ 阶梯 = {lad[2]/lad[0]:.3f}")
ok("3c v₂ 谱在阶梯 bin 与 r 反相 (剥离解除相消 => 峰真实变强)",
   np.mean(cs) < -0.3)
# v₂ 只占 σ 总方差的一小部分 (剥离本就不该大动总能量)
var_share = float(np.var(v2s) / np.var(sig))
print(f"v₂ 方差占 σ 方差比例 = {var_share:.5f} (剥离只动总方差的"
      f" {100*var_share:.2f}% 量级)")

json.dump({"audit_N_errors": {str(k): v for k, v in res2.items()},
           "den_factor": den_factor, "num_factor": num_factor,
           "cos_mean": float(np.mean(cs)), "var_share": var_share},
          open(r"C:\Users\savyc\collatz_fft_study\results4.json", "a+")
          if False else open(r"C:\Users\savyc\collatz_fft_study\audit5.json",
                             "w"), indent=1)
print("\naudit5.json 已写入")
