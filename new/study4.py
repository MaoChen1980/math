# -*- coding: utf-8 -*-
"""Round5: p-adic 处理(对齐结构) -> FFT 找规律.

信号 = results2.json 已验证采样 (24/24), 不重算. 管线特征导向:
  Stage A (p-adic 对齐):
    A1 2-adic 奇核剥离: core(n)=n>>v2(n); 恒等式 σ(n)=σ(core)+v₂ (精确)
       => r[j]=σ(core(j+1)): 所有偶数被对齐到奇核, 可证 dyadic 部分被剥掉
    A2 v₂ 分层: 层计数 c_v = floor(N/2^v)-floor(N/2^{v+1}) (精确)
    A3 3-adic 分层: 奇核信号按 n mod 6 分层 (3-结构 = 困难所在, S-F3/S-J2)
  Stage B (FFT on aligned):
    B1 DFT 相位递归: X_r[k]=E[k]+ω^k·X_r^{(N/2)}[k mod N/2] (精确)
       => σ 谱的全部内容 = 奇数 n 上 σ 值的各级 DFT + 可计算相位
    B2 阶梯定位: 原谱 625-格离散峰 vs 剥离后 E 谱 (S-F3 的谱学定位)
    B3 层均值恒等式: E[σ|v₂=v] = v + mean(σ over odd m ≤ N/2^v) (S-E1 机制定量化)
    B4 mod-6 分层谱浓度 + M01 随机对照
"""
import json
import math
import numpy as np

R = json.load(open(r"C:\Users\savyc\collatz_fft_study\results2.json",
                   encoding="utf-8"))
N = 20000
ok = lambda name, cond, detail="": print(
    f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")

memo = {1: 0}


def sg(n):
    path, x = [], n
    while x not in memo:
        path.append(x)
        x = x >> 1 if x % 2 == 0 else 3 * x + 1
    v = memo[x]
    for y in reversed(path):
        v += 1
        memo[y] = v
    return memo[n]


def v2(n):
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def core(n):
    while n % 2 == 0:
        n //= 2
    return n


print("=" * 72)
print("Stage A: p-adic 对齐")
print("=" * 72)

sig16 = np.array(R["M16"]["entries"][0]["signal"], float)  # σ(n), n=1..20000
# --- A1 奇核剥离 (精确恒等式) ---
r = np.array([sg(core(n)) for n in range(1, N + 1)], float)
v2s = np.array([v2(n) for n in range(1, N + 1)], float)
err = float(np.max(np.abs(sig16 - (r + v2s))))
ok("A1 剥离恒等式 σ(n)=σ(core)+v₂ 全窗精确 (N=20000)", err == 0.0,
   f"max_err={err:.2e}")
ok("A1b 剥离不变性 r[2m-1]=r[m-1] (j 奇指标=半指标, 0 基) 精确",
   all(r[2 * m + 1] == r[m] for m in range(N // 2)))

# --- A2 v₂ 层计数 (精确) ---
cnts = [int((v2s == v).sum()) for v in range(15)]
pred = [N // 2 ** v - N // 2 ** (v + 1) for v in range(15)]
ok("A2 层计数 c_v=floor(N/2^v)-floor(N/2^{v+1}) 精确", cnts == pred,
   f"前6层={cnts[:6]} (几何减半)")

print()
print("=" * 72)
print("Stage B: FFT on aligned")
print("=" * 72)

# --- B1 DFT 相位递归 (可证, 全 k 机器验证) ---
X = np.fft.fft(r)
E = np.fft.fft(r[0::2])            # 奇数 n 的 σ 值 DFT (10000 点)
Xh = np.fft.fft(r[:N // 2])        # r 前半 (10000 点)
omega = np.exp(-2j * np.pi / (N // 2))
rec = E + omega ** np.arange(N // 2) * Xh
err = float(np.max(np.abs(X - np.concatenate([rec, rec]))))
# ω_{N}^{k} 对 k 与 k+N/2 同 (ω_{N/2}), 且 E, Xh 以 N/2 为周期 => X 前后半相同?
# 注意: 实为 X[k] = E[k mod N/2] + ω_N^k·Xh[k mod N/2], ω_N^{k+N/2} = -ω_N^k
wN = np.exp(-2j * np.pi / N)
k = np.arange(N)
rec = E[k % (N // 2)] + wN ** k * Xh[k % (N // 2)]
err = float(np.max(np.abs(X - rec)))
ok("B1 相位递归 X[k]=E[k'] + ω^k·Xh[k'], k'=k mod N/2 全 k 精确", err < 1e-6,
   f"max_err={err:.2e}")
# 两级展开: Xh[k'] = E[k''] + ω_{N/2}^{k'}·Xq[k''], k''=k' mod N/4
Xq = np.fft.fft(r[:N // 4])
E2 = np.fft.fft(r[:N // 2][0::2])
wN2 = np.exp(-2j * np.pi / (N // 2))
rec2 = E[k % (N // 2)] + wN ** k * (
    E2[(k % (N // 2)) % (N // 4)]
    + wN2 ** (k % (N // 2)) * Xq[(k % (N // 2)) % (N // 4)])
err = float(np.max(np.abs(X - rec2)))
ok("B1b 两级展开 (E, E2, Xq + 复合相位) 精确", err < 1e-6, f"max_err={err:.2e}")

# --- B2 阶梯定位: 原谱 vs 剥离谱的 625-格离散峰 ---
t = np.arange(N)
sl, ic = np.polyfit(t, sig16, 1)
res_s = sig16 - (sl * t + ic)
sl2, ic2 = np.polyfit(t, r, 1)
res_r = r - (sl2 * t + ic2)
P_raw = np.abs(np.fft.rfft(res_s))[1:] ** 2
P_str = np.abs(np.fft.rfft(res_r))[1:] ** 2
med_raw, med_str = float(np.median(P_raw)), float(np.median(P_str))
peaks_raw = [kk for kk in (10000, 5000, 7500, 1250) if P_raw[kk - 1] > 10 * med_raw]
peaks_str = [kk for kk in (10000, 5000, 7500, 1250) if P_str[kk - 1] > 10 * med_str]
fr_raw = float(sum(P_raw[kk - 1] for kk in range(625, N // 2 + 1, 625)) / P_raw.sum())
fr_str = float(sum(P_str[kk - 1] for kk in range(625, N // 2 + 1, 625)) / P_str.sum())
ok("B2 剥离前: 625-格 4 离散峰 (>10x 中位数)", len(peaks_raw) == 4,
   f"峰={peaks_raw} 格能量占比={fr_raw:.4f}")
ok("B2b 否证预设: 阶梯在剥离后仍存在 => 阶梯来自奇核盖印几何, 非 +v₂ 部分",
   len(peaks_str) == 4 and fr_str > 0.5 * fr_raw,
   f"剥离后峰={peaks_str} 格占比 {fr_raw:.4f}->{fr_str:.4f}")
# B2d 奇子信号 e[j']=σ(2j'+1) 自身谱 (10000 点格): 阶梯的最终居所?
e = sig16[0::2]
sle, ice = np.polyfit(np.arange(N // 2), e, 1)
P_e = np.abs(np.fft.rfft(e - (sle * np.arange(N // 2) + ice)))[1:] ** 2
med_e = float(np.median(P_e))
latt_e = [kk for kk in (5000, 2500, 1250, 625) if P_e[kk - 1] > 10 * med_e]
fr_e = float(sum(P_e[kk - 1] for kk in range(625, N // 4 + 1, 625)) / P_e.sum())
ok("B2d 奇子信号 e 自身谱测量 (10000 格, 多重数 625/1250/2500/5000)",
   True,
   f">10x中位数峰={latt_e}; 625-格占比={fr_e:.4f}; top峰k="
   f"{[int(kk) + 1 for kk in np.argsort(P_e)[::-1][:5]]}")


def top_share(P, n=32):
    return float(np.sort(P)[::-1][:n].sum() / P.sum())


def lagcorr(x, lag):
    return float(np.corrcoef(x[:-lag], x[lag:])[0, 1])


print(f"\nB2c 对照表 (原始σ vs 剥离r vs 奇子信号e, 去趋势):")
res_e = e - (sle * np.arange(N // 2) + ice)
print(f"{'信号':<10}{'top32占比':>10}{'625格占比':>10}{'ρ(4096)':>10}")
for nm, P, res in (("原始 σ", P_raw, res_s), ("剥离 r", P_str, res_r),
                   ("奇子信号 e", P_e, res_e)):
    fr = {"原始 σ": fr_raw, "剥离 r": fr_str, "奇子信号 e": fr_e}[nm]
    print(f"{nm:<10}{top_share(P):>10.4f}{fr:>10.4f}"
          f"{lagcorr(res, 4096):>10.4f}")

# --- B3 层均值恒等式 (S-E1 机制定量化) ---
bad = []
for v in range(0, 8):
    mask = v2s == v
    obs = float(sig16[mask].mean())
    cv = int(mask.sum())
    th = v + float(np.mean([sg(m) for m in range(1, N // 2 ** v + 1, 2)]))
    if abs(obs - th) > 1e-9:
        bad.append(v)
ok("B3 E[σ|v₂=v] = v + mean(σ over odd m ≤ N/2^v) 全层精确 (v<8)", not bad,
   f"不符层={bad}")
print("  层均值示例: " + ", ".join(
    f"v={v}:{sig16[v2s == v].mean():.2f}" for v in range(6)))

# --- B4 3-adic 分层: 奇核信号按 n mod 6 分层谱浓度 + 对照 ---
print("\nB4 奇核信号 r 按 n mod 6 分层 (层内去趋势) top-10 谱浓度:")
for a in (1, 3, 5):
    idx = np.arange(N) % 6 == a
    x = r[idx]
    tt = np.arange(len(x))
    slp, icp = np.polyfit(tt, x, 1)
    P = np.abs(np.fft.rfft(x - (slp * tt + icp)))[1:] ** 2
    print(f"  n≡{a} (mod 6): 样本={len(x)}, top-10 占比={top_share(P, 10):.4f}")
rng = np.random.default_rng(11)
x = rng.normal(size=N // 2)
P = np.abs(np.fft.rfft(x))[1:] ** 2
print(f"  高斯对照:            top-10 占比={top_share(P, 10):.4f}")

# --- B5 逐方法剥离 (值=σ(n) 且 n 可恢复的方法; 其余如实排除) ---
# 剥离 = 方法自身采样值减 v₂(被采样整数): σ(n)=σ(core(n))+v₂(n) 逐值应用
print("\nB5 逐方法: 剥离前 -> 剥离后 (去趋势, 625-格占比 / top-32 占比):")
print(f"{'方法':<6}{'raw格':>8}{'strip格':>9}{'raw32':>8}{'strip32':>9}")
comp5 = {}
NS = {"M03": [RES[i // 2500] + 16 * (i % 2500) for i in range(N)]
      if (RES := [1, 3, 5, 7, 9, 11, 13, 15]) else None,
      "M04": [1 + 3 * i for i in range(N)],
      "M16": list(range(1, N + 1))}
for mid, ns in NS.items():
    s = np.array(R[mid]["entries"][0]["signal"], float)[:N]
    tt = np.arange(N)
    slp, icp = np.polyfit(tt, s, 1)
    s_d = s - (slp * tt + icp)
    rr = s - np.array([float(v2(n)) for n in ns])
    slp2, icp2 = np.polyfit(tt, rr, 1)
    rr_d = rr - (slp2 * tt + icp2)
    P1 = np.abs(np.fft.rfft(s_d))[1:] ** 2
    P2 = np.abs(np.fft.rfft(rr_d))[1:] ** 2
    f1 = float(sum(P1[kk - 1] for kk in range(625, N // 2 + 1, 625)) / P1.sum())
    f2 = float(sum(P2[kk - 1] for kk in range(625, N // 2 + 1, 625)) / P2.sum())
    comp5[mid] = {"raw_ladder": f1, "strip_ladder": f2,
                  "raw32": top_share(P1), "strip32": top_share(P2)}
    print(f"{mid:<6}{f1:>8.4f}{f2:>9.4f}{top_share(P1):>8.4f}"
          f"{top_share(P2):>9.4f}")
# M02: O(n) 的核心不变性定理 (O(2m)=O(m) => O(core)=O, 剥离是恒等)
s02 = np.array(R["M02"]["entries"][0]["signal"], float)[:N]
ok("B5b M02: O(core)=O 全窗精确 (O 本身 v₂ 不变, 阶梯=盖印几何非 v₂ 部分)",
   all(s02[n - 1] == s02[core(n) - 1] for n in range(1, N + 1)))
print("(排除: M01 n 未存档 / M18 值为 log-比 / M05,M13,M17,M20 非格点序信号;"
      " 对这些方法硬套剥离无语义)")

json.dump({"B1_maxerr": err, "ladder": {"raw": fr_raw, "strip": fr_str},
           "peaks_raw": peaks_raw, "peaks_strip": peaks_str,
           "per_method": comp5},
          open(r"C:\Users\savyc\collatz_fft_study\results4.json", "w"),
          indent=1)
print("\nresults4.json 已写入")
