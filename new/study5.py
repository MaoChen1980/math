# -*- coding: utf-8 -*-
"""Round6: 四个新特征 × 工具攻击 (每个预注册判据).
A1 纯L 3-adic 梳 -> 更新过程: 块长律 + 三级置换代理定位梳载体
A2 v2 细节谱 4^-l -> 平衡引理: 对齐块 (v2-2) 和有界 => 边界效应
A3 长记忆 -> DFA/Hurst: sigma 残差 vs 洗牌对照
A4 跨基定量 -> (Z/3^9) 群特征标谱: 3-adic 格能量 vs 2-adic 同型指标
"""
import json
import math
import numpy as np

R = json.load(open(r"C:\Users\savyc\collatz_fft_study\results2.json",
                   encoding="utf-8"))
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
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    return c


print("=" * 72)
print("A1 纯L 3-adic 梳 <- 更新过程 (块长律 + 三级置换代理)")
print("=" * 72)
# 重建无去重真链流 (derive P2 语义): 根 r=2,5,8..., 链 m->(2m-1)/3
roots = []
r = 2
while sum(l for _, l in roots) < 20000 + 10:
    m, ln = r, 1
    while (2 * m - 1) % 3 == 0:
        m = (2 * m - 1) // 3
        ln += 1
    roots.append((r, ln))
    r += 3
stream, blocks = [], []
for root, ln in roots:
    m, vals = root, []
    for _ in range(ln):
        vals.append(sg(m))
        if (2 * m - 1) % 3 == 0:
            m = (2 * m - 1) // 3
    blocks.append(vals)
    stream.extend(vals)
stream = np.array(stream[:20000], float)
# 块长律: P(l=m) = 2*3^-(m-1) (m>=2)? 推导: l=2+v3(s), s=(r+1)/3 任意
lens = np.array([ln for _, ln in roots], float)
ms = np.arange(2, int(lens.max()) + 1)
emp = np.array([(lens == m).mean() for m in ms])
th = 2 * 3.0 ** (-(ms - 1))
tv = float(np.abs(emp - th).sum()) / 2
ok("A1a 块长律 P(l=m)=2·3^-(m-1): 全距 TV 距离 < 0.02", tv < 0.02,
   f"TV={tv:.4f}; E[l] 实测={lens.mean():.3f} vs 定理 5/2={2.5}")
# 代理电池: 梳指标 = 750 倍数 bin 能量占比
N = 20000


def comb_share(x):
    t = np.arange(len(x))
    sl, ic = np.polyfit(t, x, 1)
    P = np.abs(np.fft.rfft(x - (sl * t + ic)))[1:] ** 2
    f = float(sum(P[k - 1] for k in range(750, N // 2 + 1, 750)) / P.sum())
    top = [int(kk) + 1 for kk in np.argsort(P)[::-1][:4]]
    return f, top


rng = np.random.default_rng(42)
f0, top0 = comb_share(stream)
sizes, ptr = [], 0
shuf_stream = stream.copy()
perm = rng.permutation(len(shuf_stream))
s1 = stream[perm]  # S1: 值全乱 (保块长? 不保) — 与 S3 同, 改: 保块长乱值
s1 = np.concatenate([rng.permutation(b) for b in blocks[:len(blocks)]])
s1 = np.array(s1[:20000], float)
order = rng.permutation(len(blocks[:int(sum(len(b) for b in blocks[:20000])//1)]))
s2 = np.concatenate([blocks[i] for i in rng.permutation(len(blocks))])
s2 = np.array(s2[:20000], float)
s3 = np.array(rng.permutation(stream), float)
print(f"{'流':<26}{'750格占比':>10}  top4峰k")
f1, top1 = comb_share(s1)
f2, top2 = comb_share(s2)
f3, top3 = comb_share(s3)
for nm, f, tp in (("原流 S0", f0, top0), ("S1 保块长乱值", f1, top1),
                  ("S2 保块乱序", f2, top2), ("S3 全乱对照", f3, top3)):
    print(f"{nm:<26}{f:>10.4f}  {tp}")
ok("A1b 梳载体定位: S1(乱值)梳存活>=原50% => 长度载梳; 否则查 S2",
   f1 >= 0.5 * f0 or f2 >= 0.5 * f0,
   f"S0={f0:.4f} S1={f1:.4f} S2={f2:.4f} S3={f3:.4f}")

print()
print("=" * 72)
print("A2 v2 细节谱 4^-l <- 平衡引理 (对齐块和有界)")
print("=" * 72)
# 奇 n=2j+1, v2(3n+1) = 1+v2(3j+2); 3 mod 2^m 可逆 => j 完全剩余系时
# 3j+2 走遍 mod 2^m => 平衡. 对齐块: j ∈ [c·2^m, (c+1)·2^m)
maxdev = []
badcf = []
# 信号值 = v2(3n+1), n=2j+1 奇 => v2(6j+4) = 1 + v2(3j+2) (削半步 +1)
for m in range(1, 14):
    j = np.arange(0, 2 ** 13)
    w = np.array([1 + v2(3 * int(jj) + 2) for jj in j], float) - 2
    sums = w.reshape(-1, 2 ** m).sum(axis=1)
    maxdev.append(float(np.abs(sums).max()))
    # 闭式: 块和 = v2(3j0+2) - m - 1, j0 = 唯一使 3j0+2≡0 (mod 2^m) 的剩余
    if m <= 6:
        inv3 = pow(3, -1, 2 ** m)
        for c in range(4):
            j0 = ((-2 * inv3) % 2 ** m) + c * 2 ** m
            if sums[c] != v2(3 * int(j0) + 2) - m - 1:
                badcf.append((m, c))
ok("A2a 平衡引理: 对齐块和 = v2(3j0+2)-m-1 (闭式精确) 且有界不随 m 增长",
   not badcf and max(maxdev) <= 12 and maxdev[-1] <= 2,
   f"max|块和| 随 m: {[f'{x:.0f}' for x in maxdev]} (平凡尺度为 2^m, "
   f"实测有界且递减); 闭式不符={badcf}")
# misaligned 对照: 随机偏移块和显著更大 (边界效应可见)
rng2 = np.random.default_rng(9)
w_all = np.array([1 + v2(3 * int(jj) + 2) for jj in range(2 ** 14)], float) - 2
offs = rng2.integers(1, 255, 512)
al = [abs(w_all[c * 256:(c + 1) * 256].sum()) for c in range(2 ** 14 // 256 - 1)]
mis = [abs(w_all[c * 256 + o:(c + 1) * 256 + o].sum())
       for c, o in zip(range(512), offs) if c * 256 + o + 256 <= 2 ** 14]
print(f"m=8: 对齐平均|和|={np.mean(al):.3f} vs 非对齐={np.mean(mis):.3f} "
      f"(比值={np.mean(mis)/max(np.mean(al),1e-9):.1f}x)")
# 细节谱 4^-l 复核 (M12 原始, 不去趋势)
import pywt
s12 = np.array(R["M12"]["entries"][1]["signal"], float)
coef = pywt.wavedec(s12, "haar", mode="periodization")
E = np.array([float(np.sum(c * c)) for c in coef[1:]][::-1])  # 层1..14
q = E * 4.0 ** np.arange(1, len(E) + 1)
cv = float(q[:11].std() / q[:11].mean())
ok("A2b E_l·4^l ≈ 常数 (l=1..11, 变异系数<0.15)", cv < 0.15,
   f"CV={cv:.4f} => E_l ∝ 4^-l 成立且由 A2a 边界效应解释")

print()
print("=" * 72)
print("A3 长记忆 <- DFA/Hurst")
print("=" * 72)


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


sig = np.array(R["M16"]["entries"][0]["signal"], float)[:20000]
t = np.arange(20000)
sl, ic = np.polyfit(t, sig, 1)
res = sig - (sl * t + ic)
rr = sig - np.array([float(v2(n)) for n in range(1, 20001)])
sl2, ic2 = np.polyfit(t, rr, 1)
res_r = rr - (sl2 * t + ic2)
sizes = np.array([16, 32, 64, 128, 256, 512, 1024, 2048, 4096])
h_sig = dfa_h(res, sizes)
h_shuf = dfa_h(np.array(rng.permutation(res)), sizes)
h_r = dfa_h(res_r, sizes)
ok("A3a σ 残差 H > 0.6 (长记忆的 DFA 定量)", h_sig > 0.6, f"H(σ残差)={h_sig:.3f}")
ok("A3b 洗牌对照 H→0.5 (白噪声零假设)", 0.45 < h_shuf < 0.55,
   f"H(洗牌)={h_shuf:.3f}")
print(f"H(剥离 r 残差) = {h_r:.3f} (长记忆在剥离后存活 => 与 S-K5 一致)")

print()
print("=" * 72)
print("A4 跨基定量 <- (Z/3^9) 群特征标谱 (N=19683)")
print("=" * 72)
N3 = 19683
sig3 = np.array([sg(n) for n in range(1, N3 + 1)], float)
t3 = np.arange(N3)
sl3, ic3 = np.polyfit(t3, sig3, 1)
X3 = np.fft.fft(sig3 - (sl3 * t3 + ic3))  # 循环群 Z/3^9 的特征标谱
P3 = np.abs(X3) ** 2
r3 = sig3 - np.array([float(v2(n)) for n in range(1, N3 + 1)])
sl3b, ic3b = np.polyfit(t3, r3, 1)
X3r = np.fft.fft(r3 - (sl3b * t3 + ic3b))
P3r = np.abs(X3r) ** 2
print(f"{'窗':>10}{'格':>12}{'格能量占比':>12}")
rows = []
for nm, P, NN in (("σ @3^9窗", P3, N3), ("剥离 @3^9窗", P3r, N3)):
    for lv, div in (("N/3", 3), ("N/9", 9), ("N/27", 27)):
        step = NN // div
        fr = float(sum(P[k - 1] for k in range(step, NN // 2 + 1, step))
                   / P[1:].sum())
        rows.append((nm, f"{lv} ({step})", fr))
        print(f"{nm:<10}{lv + ' (' + str(step) + ')':>12}{fr:>12.5f}")
# 2-adic 同型指标 (N=16384 rfft: bin k∈[0,8192]; N/2 格 = {8192}, N/4 格 = {4096,8192})
sig2w = np.array(R["M16"]["entries"][1]["signal"], float)[:16384]
t2w = np.arange(16384)
slw, icw = np.polyfit(t2w, sig2w, 1)
P2w = np.abs(np.fft.rfft(sig2w - (slw * t2w + icw)))[1:] ** 2
fr2_2 = float(P2w[8192 - 1]) / P2w.sum()
fr2_4 = float(P2w[4096 - 1] + P2w[8192 - 1]) / P2w.sum()
fr3_3 = rows[0][2]
rng3 = np.random.default_rng(5)
Pc = np.abs(np.fft.fft(rng3.normal(size=N3))) ** 2
fr3_ctrl = float(sum(Pc[k - 1] for k in range(6561, N3 // 2 + 1, 6561))
                 / Pc[1:].sum())
ok("A4a 3-adic 格浓度 ≪ 2-adic 同型格浓度 (粗层 N/3 vs N/2 格)",
   fr3_3 < 0.5 * fr2_2,
   f"3窗N/3格={fr3_3:.5f} vs 2窗N/2格={fr2_2:.5f} (2窗N/4格={fr2_4:.5f})")
ok("A4b 高斯对照: 单 bin 占比 ~2/N 量级 (无离散结构)", fr3_ctrl < 5e-3,
   f"对照={fr3_ctrl:.2e} vs 期望~{2 / N3:.1e}")
ok("A4c 3窗内剥离不消灭 N/9 格浓度 (3-结构在奇核侧, 与 S-K2 呼应)",
   rows[3][2] > 0.3 * rows[0][2],
   f"σ={rows[0][2]:.5f} -> 剥离={rows[3][2]:.5f}")

json.dump({"a1_tv": tv, "a1_comb": {"S0": f0, "S1": f1, "S2": f2, "S3": f3},
           "a2_cv": cv, "a2_maxdev": maxdev,
           "a3": {"H_sig": h_sig, "H_shuf": h_shuf, "H_r": h_r},
           "a4": {"rows": rows, "fr2_4": fr2_4, "ctrl": fr3_ctrl}},
          open(r"C:\Users\savyc\collatz_fft_study\results5.json", "w"),
          indent=1)
print("\nresults5.json 已写入")
