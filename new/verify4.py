# -*- coding: utf-8 -*-
"""Round5 独立复核: 与 study4.py 不同代码路径重算关键结论.
独立性: 递归 lru_cache sigma / 试除 v2 / 直和 DFT (位运算相位) / 独立阶梯指标."""
import json
import math
import random
from functools import lru_cache

import numpy as np

R = json.load(open(r"C:\Users\savyc\collatz_fft_study\results2.json",
                   encoding="utf-8"))
R4 = json.load(open(r"C:\Users\savyc\collatz_fft_study\results4.json",
                    encoding="utf-8"))
N = 20000
ok = lambda name, cond, detail="": print(
    f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")


@lru_cache(maxsize=None)
def sg2(n):
    if n == 1:
        return 0
    return 1 + sg2(n // 2 if n % 2 == 0 else 3 * n + 1)


def v2b(n):
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def coreb(n):
    while n % 2 == 0:
        n //= 2
    return n


sig16 = np.array(R["M16"]["entries"][0]["signal"], float)

# ---------- V1. 剥离恒等式: 递归 sigma + 试除 v2 (300 随机例) ----------
random.seed(2026)
bad = [n for n in random.sample(range(1, N + 1), 300)
       if sig16[n - 1] != sg2(coreb(n)) + v2b(n)]
ok("V1 σ(n)=σ(core)+v₂ 独立重算 (300例)", not bad, f"不符={len(bad)}")

# ---------- V2. 相位递归: 直和 DFT 独立重算 (选 5 个 k) ----------
r = np.array([sg2(coreb(n)) for n in range(1, N + 1)], float)
e = r[0::2]
Xh = np.fft.fft(r[:N // 2])
E = np.fft.fft(e)
Nh = N // 2
bad = []
for k in (0, 1, 625, 5000, 12345):
    j = np.arange(N)
    ph = np.exp(-2j * np.pi * ((k % N) * j % N) / N)
    direct = complex(np.dot(r, ph))
    kh = k % Nh
    rhs = E[kh] + np.exp(-2j * np.pi * k / N) * Xh[kh]
    if abs(direct - rhs) > 1e-4:
        bad.append(k)
ok("V2 相位递归 X[k]=E[k']+ω^k·Xh[k'] 直和重算 (k=0,1,625,5000,12345)",
   not bad, f"不符={bad}")

# ---------- V3. 两级展开直和 (k=625, 12345) ----------
Xq = np.fft.fft(r[:N // 4])
E2 = np.fft.fft(r[:N // 2][0::2])
bad = []
for k in (625, 12345):
    j = np.arange(N)
    ph = np.exp(-2j * np.pi * ((k % N) * j % N) / N)
    direct = complex(np.dot(r, ph))
    kh, kq = k % Nh, k % (N // 4)
    rhs = E[kh] + np.exp(-2j * np.pi * k / N) * (
        E2[kq] + np.exp(-2j * np.pi * kh / Nh) * Xq[kq])
    if abs(direct - rhs) > 1e-4:
        bad.append(k)
ok("V3 两级展开直和 (k=625,12345)", not bad, f"不符={bad}")

# ---------- V4. 层均值恒等式 (v=2,3 独立计数) ----------
bad = []
for v in (2, 3):
    members = [n for n in range(1, N + 1) if v2b(n) == v]
    obs = float(np.mean([sig16[n - 1] for n in members]))
    th = v + float(np.mean([sg2(m) for m in range(1, N // 2 ** v + 1, 2)]))
    if abs(obs - th) > 1e-9:
        bad.append(v)
ok("V4 层均值恒等式独立重算 (v=2,3)", not bad, f"不符层={bad}")

# ---------- V5. 阶梯指标复算 (独立实现, 全量循环) ----------
t = np.arange(N)
sl, ic = np.polyfit(t, sig16, 1)
res_s = sig16 - (sl * t + ic)
sl2, ic2 = np.polyfit(t, r, 1)
res_r = r - (sl2 * t + ic2)
Xs = np.fft.rfft(res_s)
Xr = np.fft.rfft(res_r)
# 注意 bin 约定: 625-格 = DFT bin k 的 625 倍数 (Xs[k]), 与 study4 的 P[k-1] ([1:] 位移) 同一
fr1 = sum(float(Xs[k].real ** 2 + Xs[k].imag ** 2)
          for k in range(625, N // 2 + 1, 625)) / float(np.sum(np.abs(Xs) ** 2))
fr2 = sum(float(Xr[k].real ** 2 + Xr[k].imag ** 2)
          for k in range(625, N // 2 + 1, 625)) / float(np.sum(np.abs(Xr) ** 2))
ok("V5 阶梯占比复算: raw/strip 与 results4 一致",
   abs(fr1 - R4["ladder"]["raw"]) < 1e-9 and
   abs(fr2 - R4["ladder"]["strip"]) < 1e-9,
   f"raw={fr1:.4f} strip={fr2:.4f}")
ok("V5b 否证复认: 剥离后阶梯占比不降反微升", fr2 > fr1,
   f"{fr1:.4f} -> {fr2:.4f}")

# ---------- V6. O 核心不变性抽查 ----------
s02 = np.array(R["M02"]["entries"][0]["signal"], float)[:N]
idx = random.sample(range(1, N + 1), 300)
ok("V6 O(core)=O 独立抽查 (300例)",
   all(s02[n - 1] == s02[coreb(n) - 1] for n in idx))

# ---------- V7. 奇子信号 e 谱峰复算 ----------
re_ = e - np.polyval(np.polyfit(np.arange(Nh), e, 1), np.arange(Nh))
Pe = np.abs(np.fft.rfft(re_))[1:] ** 2
med = float(np.median(Pe))
pk = [kk for kk in (5000, 2500, 1250, 625) if Pe[kk - 1] > 10 * med]
ok("V7 奇子信号 e 谱: 625/1250/2500/5000 全部 >10x 中位数", pk == [5000, 2500, 1250, 625],
   f"峰={pk}")

print("\nRound5 独立复核完毕")
