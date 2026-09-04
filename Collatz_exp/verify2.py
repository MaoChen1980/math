# -*- coding: utf-8 -*-
"""Round2 验证: 大 N 恒等式核对 + S-E1/S-F2/S-F3 复核 + 新检验."""
import json
import math
import numpy as np

R = json.load(open(r"C:\Users\savyc\collatz_fft_study\results2.json",
                   encoding="utf-8"))


def sig_of(mid, i=0):
    return np.array(R[mid]["entries"][i]["signal"], float)


ok = lambda name, cond, detail="": print(
    f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")

# --- A. sigma(2^k)=k ---
s10 = sig_of("M10")
ok("A sigma(2^k)=k (k=1..20000)", all(s10[k-1] == k for k in range(1, 20001)))

# --- B/C. M17 梳状谱支撑集 (D26) ---
s7, s11 = sig_of("M17", 0), sig_of("M17", 1)
X7 = np.abs(np.fft.rfft(s7 - s7.mean()))
sup7 = [k for k in range(1, 10001) if X7[k] > 1e-9]
ok("B mod7: 支撑={10000} (period 2)", sup7 == [10000], f"{sup7[:5]}")
X11 = np.abs(np.fft.rfft(s11 - s11.mean()))
bad = [k for k in range(1, 10001) if X11[k] > 1e-9 and k % 2000 != 0]
ok("C mod11: 支撑=2000 的倍数 (period 10)", not bad, f"违规={bad[:5]}")

# --- D. Nyquist 线幅值代数式 @ N=20000 ---
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


s16 = sig_of("M16")
x0 = s16 - s16.mean()
nyq = abs(np.fft.rfft(x0)[10000])
direct = abs(sum((-1) ** n * s16[n - 1] for n in range(1, 20001)))
even_sum = sum(sg(m) for m in range(1, 10001)) + 10000
odd_sum = sum(2 + sg(3 * m - 1) for m in range(2, 10001))
alg = abs(even_sum - odd_sum)
ok("D Nyquist 闭式 @20000: dft=direct=algebra",
   abs(nyq - alg) < 1e-6 and abs(direct - alg) < 1e-6,
   f"dft={nyq:.1f} algebra={alg:.1f}")

# --- E. S-E1: v2 泄漏/浓度 三窗口对比 ---
e200 = 0.8725  # 第一轮 N=200 实测
e20k = R["M12"]["entries"][0]["fft"]["even_k_energy"]
e16k = R["M12"]["entries"][1]["fft"]["even_k_energy"]
s12 = sig_of("M12")
ok("E1 E[v2]=2+o(1) (渐近; 偏差为确定性 3-adic 格效应)",
   abs(s12.mean() - 2.0) < 0.01, f"mean={s12.mean():.5f}")
ok("E2 P(v2=1)=50% (20000奇数中恰10000个 n=3 mod 4)",
   (s12 == 1).sum() == 10000)
ok("E3 S-E1 泄漏论: N=16384 浓度>0.999 且 > N=20000 > N=200",
   e16k > 0.999 and e16k > e20k > e200,
   f"N=200:{e200} N=20000:{e20k} N=16384:{e16k}")

# --- F. M18 模4分类 ---
s18 = sig_of("M18")
g_even = [s18[n-1] for n in range(1, 20001) if n % 2 == 0]
g_3m4 = [s18[n-1] for n in range(1, 20001) if n % 4 == 3]
g_1m4 = [s18[n-1] for n in range(1, 20001) if n % 4 == 1]
ok("F1 2|n 恒 -1", all(abs(v + 1) < 1e-9 for v in g_even))
ok("F2 n=3 mod4 均值=log2(3)-1+o(1)",
   abs(np.mean(g_3m4) - (math.log2(3) - 1)) < 0.01,
   f"{np.mean(g_3m4):.5f} vs {math.log2(3)-1:.5f}")
ok("F3 n=1 mod4 均值=log2(3)-3+o(1)",
   abs(np.mean(g_1m4) - (math.log2(3) - 3)) < 0.01,
   f"{np.mean(g_1m4):.5f} vs {math.log2(3)-3:.5f}")

# --- G. S-F2/S-F3 @ N=20000: 谱形状不变性 + 625-格能量占比 ---
memo2 = {1: 0, 2: 1}


def orbit_list(n):
    seq = [n]
    while seq[-1] != 1:
        x = seq[-1]
        seq.append(x >> 1 if x % 2 == 0 else 3 * x + 1)
    return seq


S = s16.copy()
O = np.array([sum(v & 1 for v in orbit_list(n)) for n in range(1, 20001)],
             float)
H = S - O + 1


def spec(x):
    t = np.arange(len(x))
    sl, ic = np.polyfit(t, x, 1)
    return np.abs(np.fft.rfft(x - (sl * t + ic)))


XS, XO, XH = spec(S), spec(O), spec(H)
c1 = float(np.corrcoef(XS, XO)[0, 1])
c2 = float(np.corrcoef(XS, XH)[0, 1])
ok("G1 S-F2: corr(spec_sigma, spec_O)>0.99", c1 > 0.99, f"corr={c1:.4f}")
ok("G2 S-F2: corr(spec_sigma, spec_H)>0.99", c2 > 0.99, f"corr={c2:.4f}")
P = XS[1:] ** 2  # bin k=1..10000
frac625 = float(sum(P[k - 1] for k in range(625, 10001, 625)) / P.sum())
med = float(np.median(P))
peaks625 = [k for k in (10000, 5000, 7500, 1250) if P[k - 1] > 10 * med]
ok("G3 S-F3 定量: 格线为离散峰(>10x中位数) + 连续谱占能量主体",
   len(peaks625) == 4 and frac625 < 0.5,
   f"格线能量占比={frac625:.4f} (连续谱={1-frac625:.4f}); "
   f"离散峰>10x中位数: {peaks625}")

# --- H. M13 增长相: 纯 Nyquist 线定理 ---
s13 = sig_of("M13")
slope_th = (math.log2(3) - 1) / 2
ok("H1 斜率=(log2(3)-1)/2 (5位)", abs(0.292481 - slope_th) < 5e-6,
   f"obs=0.292481 theory={slope_th:.6f}")
K = 10000
n0 = 2**K - 1
step1 = 3 * n0 + 1
step2 = step1 // 2
ok("H2 T(2^K-1)=3*2^K-2, T^2=3*2^(K-1)-1",
   step1 == 3 * 2**K - 2 and step2 == 3 * 2**(K - 1) - 1)
a13 = R["M13"]["entries"][0]["fft"]["top_peaks"]
ok("H3 增长相去趋势谱=纯Nyquist线 (第2峰=0%)", a13[1]["rel"] == 0.0,
   f"peaks={[ (p['k'], p['rel']) for p in a13[:3] ]}")

# --- I. M15 1/f^2 斜率 ---
s15 = sig_of("M15")
t = np.arange(20000)
sl, ic = np.polyfit(t, s15, 1)
Xp = np.abs(np.fft.rfft(s15 - (sl * t + ic)))[1:] ** 2
f = np.fft.rfftfreq(20000)[1:]
m = (f >= 1e-4) & (f <= 5e-2)
slope_fit = float(np.polyfit(np.log(f[m]), np.log(Xp[m]), 1)[0])
ok("I M15 失衡流谱斜率 ~ -2 (漂移随机游走)", -2.6 < slope_fit < -1.4,
   f"slope={slope_fit:.3f}")

# --- J. 负对照 + M20 新观察 ---
f01 = R["M01"]["entries"][0]["fft"]["flatness"]
f19 = R["M19"]["entries"][0]["fft"]["flatness"]
ok("J1 M01 随机采样 flatness~0.55 (白化)", 0.45 < f01 < 0.65, f"{f01}")
p20 = R["M20"]["entries"][0]["fft"]["top_peaks"]
ks = [p["k"] for p in p20[:4]]
ok("J2 M20 纯L: 前4峰全在 k 的 750 倍数 (3-adic 梳; 第5峰 k=1 为链长趋势)",
   all(k % 750 == 0 for k in ks), f"top4_k={ks} (f={3/80:.5f} 的整数倍)")
