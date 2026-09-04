# -*- coding: utf-8 -*-
"""Round4 独立复核: 与 study3.py 不同代码路径重算关键结论.
独立性: 显式 chi 矩阵 WHT / pywt.wavedec / 位运算直和 / pow 模运算 /
递归 lru_cache sigma."""
import json
import math
from functools import lru_cache

import numpy as np
import pywt

R3 = json.load(open(r"C:\Users\savyc\collatz_fft_study\results3.json",
                    encoding="utf-8"))
N, K = 16384, 14
ok = lambda name, cond, detail="": print(
    f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")


@lru_cache(maxsize=None)
def sg2(n):
    if n == 1:
        return 0
    return 1 + sg2(n // 2 if n % 2 == 0 else 3 * n + 1)


# ---------- V1. fwht 正确性: 显式 chi 矩阵 (N=64) + 对合性 (N=16384) ----------
def fwht(x):
    a = np.asarray(x, float).copy()
    n, h = a.size, 1
    while h < n:
        a = a.reshape(-1, 2, h)
        u, v = a[:, 0, :].copy(), a[:, 1, :].copy()
        a = np.stack([u + v, u - v], axis=1).reshape(-1)
        h *= 2
    return a


def chi_matrix(n):
    j = np.arange(n)
    bits = ((j[:, None] >> np.arange(n.bit_length() - 1)) & 1)
    return (-1.0) ** (bits @ bits.T)  # W[q,j]


n64 = 64
x64 = np.arange(1, n64 + 1, dtype=float)
ok("V1a fwht vs 显式 chi 矩阵 (N=64, x=j+1)",
   np.max(np.abs(fwht(x64) - chi_matrix(n64) @ x64)) < 1e-9)
rng = np.random.default_rng(7)
big = rng.normal(size=N)
ok("V1b fwht 对合性: fwht(fwht(x))=N·x (N=16384)",
   np.max(np.abs(fwht(fwht(big)) - N * big)) < 1e-8)

# ---------- V2. ramp 闭式: 位运算直和独立重算 (选 6 个 q) ----------
a = np.array([sg2(n) for n in range(1, N + 1)], float)  # = [1..N] ramp
j = np.arange(N)
qsel = [0, 1, 2, 3, 1024, 4096, 12345]
bad = []
for q in qsel:
    qb = ((q >> np.arange(K)) & 1)
    chi = (-1.0) ** (((j[:, None] >> np.arange(K)) & 1) @ qb)
    direct = float(np.sum((j + 1) * chi))
    closed = N * (N + 1) / 2 if q == 0 else (
        -N * 2 ** (q.bit_length() - 2) if q & (q - 1) == 0 else 0.0)
    if abs(direct - closed) > 1e-6:
        bad.append(q)
ok("V2 ramp WHT 闭式独立直和 (q=0,1,2,3,1024,4096,12345)", not bad,
   f"不符={bad}")

# ---------- V3. 折叠公式: 小 N 显式矩阵 + 真 σ 值 ----------
n8 = 64
a8 = np.array([sg2(n) for n in range(1, n8 + 1)], float)
M = chi_matrix(n8)
We, Wo = M[:0], None
e8, o8 = a8[0::2], a8[1::2]
Mh = chi_matrix(n8 // 2)
rec = np.empty(n8)
for q in range(n8):
    qh, q0 = q >> 1, q & 1
    rec[q] = (Mh @ e8)[qh] + ((Mh @ o8)[qh] if q0 == 0 else -(Mh @ o8)[qh])
ok("V3 折叠公式 vs 显式矩阵 @ N=64 真 σ",
   np.max(np.abs(rec - M @ a8)) < 1e-9)

# ---------- V4. B1 细节恒等式 @16384 (递归 sigma 独立重算) ----------
d1 = (a[0::2] - a[1::2]) / math.sqrt(2)
rhs = np.array([(sg2(2 * k + 1) - sg2(k + 1) - 1) / math.sqrt(2)
                for k in range(N // 2)])
ok("V4 d1[k]=(σ(2k+1)-σ(k+1)-1)/√2 全窗独立重算",
   float(np.max(np.abs(d1 - rhs))) < 1e-9)

# ---------- V5. DWT 细节谱: pywt.wavedec 独立路径 (对象 = ramp [1..N]) ----------
ramp = np.arange(1, N + 1, dtype=float)
# wavedec 细节顺序: coef[1]=最粗(len 1) ... coef[K]=最细(pairs) => 反转对齐层1..K
coef = pywt.wavedec(ramp, "haar", mode="periodization")  # [approx, dK..d1]
E_wd = [float(np.sum(c * c)) for c in coef[1:]][::-1]
pred = [N * 2 ** (2 * lev - 4) for lev in range(1, K + 1)]
rel = max(abs(x - y) / y for x, y in zip(E_wd, pred))
ok("V5 ramp DWT 细节谱 E(l)=N·2^(2l-4) via pywt.wavedec", rel < 1e-9,
   f"相对误差={rel:.2e}")

# ---------- V6. W5 表抽查: M16 的 WHT 三个 q 位运算直和 ----------
s16 = np.array(json.load(open(
    r"C:\Users\savyc\collatz_fft_study\results2.json", encoding="utf-8"))
    ["M16"]["entries"][1]["signal"], float)[:N]
W16 = fwht(s16)
bad = []
for q in (0, 1, 4096, 8192):
    qb = ((q >> np.arange(K)) & 1)
    chi = (-1.0) ** (((np.arange(N)[:, None] >> np.arange(K)) & 1) @ qb)
    if abs(float(np.sum(s16 * chi)) - W16[q]) > 1e-3:
        bad.append(q)
ok("V6 M16 WHT 抽查直和 (q=0,1,4096,8192)", not bad, f"不符={bad}")

# ---------- V7. C1 p-adic: 循环除法 v2 独立 ----------
s12 = np.array(json.load(open(
    r"C:\Users\savyc\collatz_fft_study\results2.json", encoding="utf-8"))
    ["M12"]["entries"][1]["signal"], float)
idx = np.random.default_rng(3).integers(0, 8192, 300)
bad = []
for i in idx:
    m, v = 3 * (2 * int(i) + 1) + 1, 0
    while m % 2 == 0:
        m //= 2
        v += 1
    if s12[int(i)] != v:
        bad.append(int(i))
ok("V7 v2 循环除法独立重算 (300例)", not bad, f"不符={len(bad)}")
ok("V7b P(v2=1)=8192/8192 精确", int((s12 == 1).sum()) == 8192)

# ---------- V8. C4 LTE: pow 模运算独立 (v2=16 => 3^K ≡ 1+2^16 mod 2^17) ----------
r17 = pow(3, 16384, 1 << 17)
ok("V8 LTE 独立: 3^16384 ≡ 1+2^16 (mod 2^17)", r17 == 1 + (1 << 16),
   f"pow mod 2^17 = {r17}")

# ---------- V9. W5 JSON 数据与文中一致性 ----------
wc = R3["wht_compare"]
ok("V9 results3 表: M13 支撑=8 且 WHT占比=1.0",
   wc["M13"]["support"] == 8 and abs(wc["M13"]["wht32"] - 1.0) < 1e-9,
   f"M13: {wc['M13']}")
ok("V9b 20法均值 WHT≈DFT (0.2325 vs 0.2357)",
   abs(np.mean([c["wht32"] for c in wc.values()]) - 0.2325) < 1e-3)

print("\nRound4 独立复核完毕")
