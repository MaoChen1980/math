# -*- coding: utf-8 -*-
"""Round6 独立复核: 递归 sg2 链重建 / 循环除法 v2b / 直和 DFT / 手写 Haar 树."""
import json
import math
import random
from functools import lru_cache

import numpy as np

R = json.load(open(r"C:\Users\savyc\collatz_fft_study\results2.json",
                   encoding="utf-8"))
R5 = json.load(open(r"C:\Users\savyc\collatz_fft_study\results5.json",
                    encoding="utf-8"))
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


# ---------- V1. 块长律: 独立链重建 (300 根) ----------
random.seed(55)
roots = random.sample(range(2, 30000, 3), 300)
lens = []
for r in roots:
    m, ln = r, 1
    while (2 * m - 1) % 3 == 0:
        m = (2 * m - 1) // 3
        ln += 1
    lens.append(ln)
lens = np.array(lens, float)
ms = np.arange(2, int(lens.max()) + 1)
emp = np.array([(lens == m).mean() for m in ms])
th = 2 * 3.0 ** (-(ms - 1))
tv = float(np.abs(emp - th).sum()) / 2
ok("V1 块长律独立重算: TV<0.05 (300根) 且 E[l]≈2.5", tv < 0.05 and
   abs(lens.mean() - 2.5) < 0.3, f"TV={tv:.4f} E[l]={lens.mean():.3f}")

# ---------- V2. 平衡闭式: 独立 v2b + 独立块求和 ----------
bad = []
for m, c in ((3, 0), (5, 2), (6, 1)):
    inv3 = pow(3, -1, 2 ** m)
    j0 = ((-2 * inv3) % 2 ** m) + c * 2 ** m
    s = sum(1 + v2b(3 * int(j) + 2) - 2 for j in range(c * 2 ** m,
                                                       (c + 1) * 2 ** m))
    if s != v2b(3 * int(((-2 * inv3) % 2 ** m) + c * 2 ** m) + 2) - m - 1:
        bad.append((m, c))
w10 = np.array([1 + v2b(3 * int(j) + 2) - 2 for j in range(2 ** 11)], float)
b10 = float(np.abs(w10.reshape(-1, 1024).sum(axis=1)).max())
ok("V2 平衡闭式独立重算 + m=10 有界", not bad and b10 <= 5,
   f"闭式不符={bad}; m=10 max|块和|={b10:.0f} (平凡尺度 1024)")

# ---------- V3. 代理稳定性: 独立重建流, 换 seed 重复 S0/S1 ----------
def build_stream(limit):
    out, r = [], 2
    while len(out) < limit:
        m, vals = r, []
        while True:
            vals.append(sg2(m))
            if (2 * m - 1) % 3 == 0:
                m = (2 * m - 1) // 3
            else:
                break
        out.extend(vals)
        r += 3
    return out[:limit]


def comb_share(x, N=20000):
    t = np.arange(len(x))
    sl, ic = np.polyfit(t, x, 1)
    P = np.abs(np.fft.rfft(x - (sl * t + ic)))[1:] ** 2
    return float(sum(P[k - 1] for k in range(750, N // 2 + 1, 750)) / P.sum())


stream = np.array(build_stream(20000), float)
rng = np.random.default_rng(777)
perm = rng.permutation(20000)
# S1 独立版: 同长度块内置乱 (块边界重建)
bounds, acc = [], 2
bl = []
r = 2
while len(bl) < 20000:
    m, ln = r, 1
    while (2 * m - 1) % 3 == 0:
        m = (2 * m - 1) // 3
        ln += 1
    bl.append(ln)
    r += 3
bl = np.array(bl[:int(np.sum(np.cumsum(bl) <= 20000))], int)
s1 = np.concatenate([rng.permutation(b) for b in
                     np.split(stream[:int(bl.sum())], np.cumsum(bl)[:-1])])
s1 = np.array(s1[:20000], float)
f0, f1 = comb_share(stream), comb_share(s1)
ok("V3 代理稳定性(独立重建+新seed): S1/S0 梳占比比 > 0.7", f1 > 0.7 * f0,
   f"S0={f0:.4f} S1={f1:.4f}")

# ---------- V4. A4 直和抽查: bin 6561 @3^9 窗 ----------
N3 = 19683
sig3 = np.array([sg2(n) for n in range(1, N3 + 1)], float)
t3 = np.arange(N3)
sl3, ic3 = np.polyfit(t3, sig3, 1)
xd = sig3 - (sl3 * t3 + ic3)
X = np.fft.fft(xd)
j = np.arange(N3)
direct = complex(np.dot(xd, np.exp(-2j * np.pi * (6561 * j % N3) / N3)))
ok("V4 3窗 bin 6561 直和重算", abs(direct - X[6561]) < 1e-4,
   f"|diff|={abs(direct - X[6561]):.2e}")
sig2w = np.array(R["M16"]["entries"][1]["signal"], float)[:16384]
t2w = np.arange(16384)
slw, icw = np.polyfit(t2w, sig2w, 1)
P2w = np.abs(np.fft.rfft(sig2w - (slw * t2w + icw)))[1:] ** 2
fr2_2 = float(P2w[8191]) / P2w.sum()
ok("V4b 2窗 N/2 格占比复算与 results5 一致",
   abs(fr2_2 - R5["a4"]["fr2_2"] if "fr2_2" in R5["a4"] else fr2_2 - 0.03262)
   < 5e-4, f"{fr2_2:.5f}")

# ---------- V5. DFA 内部效度: iid 高斯 H≈0.5 ----------
rng2 = np.random.default_rng(21)


def dfa_h(x, sizes):
    y = np.cumsum(x - x.mean())
    F = []
    for n in sizes:
        nb = len(y) // n
        seg = y[:nb * n].reshape(nb, n)
        tt = np.arange(n)
        sl = np.polyfit(tt, seg.T, 1)[0]
        fit = sl[:, None] * tt[None, :] + \
            (seg.mean(axis=1) - sl * (n - 1) / 2)[:, None]
        F.append(float(np.sqrt(np.mean((seg - fit) ** 2))))
    return float(np.polyfit(np.log(sizes), np.log(F), 1)[0])


h_iid = dfa_h(rng2.normal(size=20000), np.array([16, 32, 64, 128, 256, 512,
                                                 1024, 2048, 4096]))
ok("V5 DFA 效度: iid 高斯 H∈[0.44,0.56]", 0.44 < h_iid < 0.56,
   f"H(iid)={h_iid:.3f}")

# ---------- V6. E_l·4^l 常数性: 手写 Haar 树独立路径 ----------
def haar_dwt_details(x, kmax):
    cur = np.asarray(x, float)
    Es = []
    for _ in range(kmax):
        a = cur.reshape(-1, 2)
        cur = (a[:, 0] + a[:, 1]) / math.sqrt(2)
        d = (a[:, 0] - a[:, 1]) / math.sqrt(2)
        Es.append(float(np.sum(d * d)))  # d1 最细在前
    return Es[::-1]  # 层1(最粗)..kmax(最细)? 反转: wavedec 细节粗->细


s12 = np.array(R["M12"]["entries"][1]["signal"], float)
E = np.array(haar_dwt_details(s12, 14))  # 层1=最粗..14=最细
E = E[::-1]  # 细->粗 对齐 A2b
q = E * 4.0 ** np.arange(1, len(E) + 1)
cv = float(q[:11].std() / q[:11].mean())
ok("V6 E_l·4^l 常数性 (手写Haar树): CV<0.15", cv < 0.15, f"CV={cv:.4f}")

print("\nRound6 独立复核完毕")
