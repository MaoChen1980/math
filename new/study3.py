# -*- coding: utf-8 -*-
"""Round4: Walsh-Hadamard + 小波包 + p-adic 特征分析 (N=16384 dyadic 窗).

信号 = results2.json 已验证采样 (24/24), 不重算. 全部计算特征导向:
  A. WHT   -- 沃尔什函数 = (Z/2)^K 特征标 => 本质是 2-adic Fourier 分析
  B. 小波包 -- Haar 尺度树, 细节系数携带精确奇偶恒等式
  C. p-adic -- 2-kernel / 2-automatic / 3-adic 类均值衰减
"""
import json
import math
import numpy as np
import pywt

R = json.load(open(r"C:\Users\savyc\collatz_fft_study\results2.json",
                   encoding="utf-8"))
N = 16384
K = 14
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


def sig(mid, i=0):
    return np.array(R[mid]["entries"][i]["signal"], float)


def sig_pow2(mid, i=0):
    """dyadic 窗: 原生 16384 用 entries[1]; 否则取不超过信号长的最大 2^K 前缀."""
    s = sig(mid, i)
    m = 1 << (len(s).bit_length() - 1)
    if mid in ("M12", "M16"):
        s = sig(mid, 1)
        m = N
    return s[:m], m


def sig16384(mid, i=0):
    return sig_pow2(mid, i)[0]


# ---------- 快速 WHT (Hadamard 自然序, 非归一) ----------
def fwht(x):
    a = np.asarray(x, float).copy()
    n, h = a.size, 1
    while h < n:
        a = a.reshape(-1, 2, h)
        u, v = a[:, 0, :].copy(), a[:, 1, :].copy()
        a = np.stack([u + v, u - v], axis=1).reshape(-1)
        h *= 2
    return a


# ---------- Haar 小波包树 (正交, Parseval 精确; 手写以便恒等式核对) ----------
def haar_packet(x, kmax):
    cur, levels = [np.asarray(x, float)], []
    for _ in range(kmax):
        nxt = []
        for arr in cur:
            a = arr.reshape(-1, 2)
            nxt.append((a[:, 0] + a[:, 1]) / math.sqrt(2))
            nxt.append((a[:, 0] - a[:, 1]) / math.sqrt(2))
        cur = nxt
        levels.append(cur)
    return levels  # levels[j-1] = 第 j 层 (2^j 节点, j=1 最粗)


def level_energy(levels):
    return [sum(float(np.sum(n * n)) for n in nodes) for nodes in levels]


def odd_chain_share(levels):
    """全 detail 子带链 (每层取 detail 子节点) 能量占比."""
    tot = sum(level_energy(levels))
    e = 0.0
    # 逐层 detail 链: 索引 2^lev-1 (全 1 的二进制路径)
    for lev, nodes in enumerate(levels, 1):
        idx = 2 ** lev - 1
        if idx < len(nodes):
            e += float(np.sum(nodes[idx] * nodes[idx]))
    return e / tot


print("=" * 72)
print("A. Walsh-Hadamard 变换 (沃尔什函数 = 2-adic 特征标)")
print("=" * 72)

# --- W1 可证: 分块常值 => WHT 支撑 ⊆ 1024 的倍数 ---
# 沃尔什函数 W_q 在高 m 位相同的位置上取同值 ⟺ q 被 2^{K-m} 整除.
# 检验对象: mod-16 分层 (M03 采样设计) 的分块常值信号 (M03 存储信号块内是
# 变化的 σ 值, 不属本定理对象 —— 第一版误用, 已修正).
x16 = np.repeat(np.arange(16, dtype=float), N // 16)
W16 = fwht(x16)
sup = sorted(np.flatnonzero(np.abs(W16) > 1e-9 * np.max(np.abs(W16))).tolist())
ok("W1 分块常值信号 (mod16 分层): 支撑 ⊆ {q: 1024|q}",
   all(q % 1024 == 0 for q in sup) and len(sup) <= 16,
   f"支撑={len(sup)}条: {[q // 1024 for q in sup]}")

# --- W2 可证: σ(2^k)=k ramp 的 WHT 闭式 ---
# x[j]=j+1: W[q] = N(N+1)/2 (q=0); -N*2^(b-1) (q=2^b); 0 其余
m10 = sig16384("M10")  # σ(2^k)=k, k=1..16384 => x[j]=j+1
W10 = fwht(m10)
closed = np.zeros(N)
closed[0] = N * (N + 1) / 2
for b in range(K):
    closed[2 ** b] = -N * 2 ** (b - 1)
err = float(np.max(np.abs(W10 - closed)))
ok("W2 ramp 闭式 W[q]=N(N+1)/2·δq0 - N·2^(b-1)·δ(q=2^b)", err < 1e-4,
   f"max_err={err:.2e}; 支撑恰={np.count_nonzero(np.abs(W10)>1e-6)} 条")

# --- W3 可证: 奇偶递归 => WHT 域折叠公式 (S-F3 的 Walsh 版能量分解) ---
# a[j]=σ(j+1); 实测 halves e[j']=a[2j'], o[j']=a[2j'+1]:
#   W(a)[q] = W(e)[q>>1] + (-1)^{q0}·W(o)[q>>1]   (N/2 点原始 WHT)
# 闭式核对: o[j']=a[j']+1 精确; e[j']=2+σ(3j'+2) 对 j'≥1 (j'=0 撞 n=1 终点,
#   与第一轮发现的 m=1 例外同一处)
m16 = sig16384("M16", 1)  # 原生 16384
Wa = fwht(m16)
e_part = m16[0::2].copy()
o_part = m16[1::2].copy()
Wo, We = fwht(o_part), fwht(e_part)
rec = np.empty(N)
for q in range(N):
    qh, q0 = q >> 1, q & 1
    rec[q] = We[qh] + (Wo[qh] if q0 == 0 else -Wo[qh])
err = float(np.max(np.abs(rec - Wa)))
ok("W3 WHT 折叠公式 W(a)[q]=W(e)[q>>1]±W(o)[q>>1] (原始WHT)", err < 1e-4,
   f"max_err={err:.2e}")
ok("W3a 闭式 o[j']=a[j']+1 精确",
   float(np.max(np.abs(o_part - m16[:N // 2] - 1))) == 0.0)
ecf = np.array([2 + sg(3 * j + 2) for j in range(1, N // 2)], float)
ok("W3a 闭式 e[j']=2+σ(3j'+2) (j'≥1)",
   float(np.max(np.abs(ecf - e_part[1:]))) == 0.0)
# Parseval 能量分解 (按各自点数归一): e 半 = 奇数 n (3-膨胀域)
Ea = float(np.sum(Wa ** 2)) / N
Eo = float(np.sum(Wo ** 2)) / (N // 2)
Ee = float(np.sum(We ** 2)) / (N // 2)
ok("W3b Parseval: Σa²=Σe²+Σo²", abs(Ea - (Ee + Eo)) / Ea < 1e-10,
   f"偏差={abs(Ea-(Ee+Eo))/Ea:.2e}; 奇数n半(3-膨胀域)能量占比={Ee/(Ee+Eo):.4f}")

# --- W4 对照: M01 随机采样白化 ---
m01 = sig16384("M01")
W01 = fwht(m01 - m01.mean())
share01 = float(np.sort(np.abs(W01))[::-1][:32].sum() / np.sum(np.abs(W01)))
m01dft = np.abs(np.fft.rfft(m01 - m01.mean()))
sh_d01 = float(np.sort(m01dft)[::-1][:32].sum() / m01dft.sum())
ok("W4 负对照 M01: WHT top-32 能量占比 < 0.5 (白化)", share01 < 0.5,
   f"WHT={share01:.4f} DFT={sh_d01:.4f}")

# --- W5 全 20 方法: WHT vs DFT 压缩性对比 (top-32 能量占比) ---
print("\nW5 逐方法对比 (top-32 能量占比, 16384 窗, 去均值/去趋势):")
print(f"{'方法':<6}{'WHT':>9}{'DFT':>9}   WHT支撑条数")
comp = {}
for mid in ["M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08", "M09",
            "M10", "M11", "M12", "M13", "M14", "M15", "M16", "M17", "M18",
            "M19", "M20"]:
    i = 1 if mid in ("M12", "M16") else 0
    x = sig16384(mid, i)
    t = np.arange(len(x))
    sl, ic = np.polyfit(t, x, 1)
    x = x - (sl * t + ic)
    W = fwht(x)
    D = np.abs(np.fft.rfft(x))
    w32 = float(np.sort(np.abs(W[1:]))[::-1][:32].sum() / np.sum(np.abs(W[1:])))
    d32 = float(np.sort(D[1:])[::-1][:32].sum() / D[1:].sum())
    nsup = int(np.count_nonzero(np.abs(W[1:]) > 1e-6 * np.max(np.abs(W[1:]))))
    comp[mid] = {"wht32": w32, "dft32": d32, "support": nsup}
    print(f"{mid:<6}{w32:>9.4f}{d32:>9.4f}   {nsup}")
wavg = np.mean([c["wht32"] for c in comp.values()])
davg = np.mean([c["dft32"] for c in comp.values()])
print(f"20法均值: WHT={wavg:.4f} vs DFT={davg:.4f}")

print()
print("=" * 72)
print("B. Haar 小波包分解 (尺度树)")
print("=" * 72)

# --- pywt 交叉确认 (用户指定小波包) + Parseval ---
wp = pywt.WaveletPacket(m16, wavelet="haar", mode="periodization",
                        maxlevel=3)
mine = haar_packet(m16, 3)
agree = True
for lev in range(1, 4):
    nodes = wp.get_level(lev, "natural")
    idx = 0
    for nd in nodes:
        path_bits = [0 if c == "a" else 1 for c in nd.path]
        i = 0
        for b in path_bits:
            i = 2 * i + b
        if not np.allclose(nd.data, mine[lev - 1][i], atol=1e-9):
            agree = False
ok("B0 pywt(haar,periodization) 与手写包树逐节点一致(1~3层)", agree)

# --- B1 可证: 细节系数精确恒等式 ---
# d1[k]=(σ(2k+1)-σ(2k+2))/√2 = (σ(2k+1)-σ(k+1)-1)/√2
d1 = mine[0][1]
rhs = np.array([(sg(2 * k + 1) - sg(k + 1) - 1) / math.sqrt(2)
                for k in range(N // 2)])
err = float(np.max(np.abs(d1 - rhs)))
ok("B1 第1层 detail: d1[k]=(σ(2k+1)-σ(k+1)-1)/√2 精确", err < 1e-9,
   f"max_err={err:.2e}")

# --- B2 可证+特征: Haar DWT 细节谱 (node index 1 = a^(l-1)d 路) ---
# 注: 包树第 l 层全部节点能量之和恒等 Σx² (每层都是完备正交基) —— 无信息量;
#     有信息量的是 a^(l-1)d 节点 = 经典 Haar DWT 尺度 l 细节能量.
# ramp 闭式: E(l) = N·2^(2l-4), l=1..K (块 2^l, 半差 h², 归一 2^(l/2))
levels10 = haar_packet(m10, K)
E10 = [float(np.sum(levels10[lev - 1][1] ** 2)) for lev in range(1, K + 1)]
pred10 = [N * 2 ** (2 * lev - 4) for lev in range(1, K + 1)]
rel = max(abs(a - b) / b for a, b in zip(E10, pred10))
ok("B2 ramp DWT 细节谱闭式 E(l)=N·2^(2l-4)", rel < 1e-9, f"相对误差={rel:.2e}")

print("\nB2b DWT 细节谱 E(l)/Σ (层1=最细相邻对...层14=最粗, %; 解析信号用原始值,"
      " 其余去趋势):")
print(f"{'方法':<6}{'E1':>7}{'E2':>7}{'E3':>7}{'E4':>7}{'E6':>7}{'E8':>7}"
      f"{'E10':>7}{'E12':>7}{'E14':>7}  lnE-lnl 斜率(层1..6)")
prof = {}
for mid in ["M10", "M16", "M02", "M04", "M12", "M18", "M20", "M01", "M11",
            "M15", "M13"]:
    if mid in ("M10", "M13"):  # 解析信号: 趋势即内容, 不去趋势
        x = m10 if mid == "M10" else sig16384("M13")
    else:
        x = sig16384(mid, 1 if mid in ("M12", "M16") else 0)
        t = np.arange(len(x))
        sl, ic = np.polyfit(t, x, 1)
        x = x - (sl * t + ic)
    lv = haar_packet(x, len(x).bit_length() - 1)
    E = np.array([float(np.sum(nodes[1] ** 2)) for nodes in lv])
    E = E / E.sum()
    slp = float(np.polyfit(np.arange(1, 7), np.log(E[:6]), 1)[0])
    prof[mid] = {"E": E.tolist(), "slope": slp}
    print(f"{mid:<6}" + "".join(f"{E[l-1]*100:>7.2f}" for l in
                                (1, 2, 3, 4, 6, 8, 10, 12, 14)) +
          f"  {slp:>7.3f}")
print("(锚点: 原始 ramp 斜率=+1.386=ln4 (E∝4^l, 细->粗); 白噪声=0; "
      "M12 实测 −1.387=−ln4 => E(l)∝4^(−l), 精确反比律, 待严格化)")

# --- B3 特征: 全-detail 链 (奇指标路径) 能量占比 ---
print("\nB3 全-detail 子带链能量占比 (奇指标路径, 越细越集中=3-adic 局域):")
oddS = {}
for mid in ["M16", "M02", "M12", "M20", "M01"]:
    i = 1 if mid in ("M12", "M16") else 0
    x = sig16384(mid, i)
    sh = odd_chain_share(haar_packet(x, K))
    oddS[mid] = sh
    print(f"  {mid}: {sh:.5f}")

print()
print("=" * 72)
print("C. p-adic 分析")
print("=" * 72)

# --- C1 2-adic: v2 序列 2-kernel 复验 @16384 (S-E2) ---
s12d = sig("M12", 1)
ok("C1 E[v2]=2 (16384 窗, 3-adic 格偏差<0.01)", abs(s12d.mean() - 2) < 0.01,
   f"mean={s12d.mean():.5f}")
ok("C1b P(v2=1)=50% 精确", int((s12d == 1).sum()) == 8192,
   f"{int((s12d == 1).sum())}/8192")

# --- C2 2-automatic 检验: σ(n) mod 2 与 O(n) mod 2 的 WHT 支撑 ---
par = np.array([sg(n) % 2 for n in range(1, N + 1)], float)


def oddcnt(n):
    c, x = 0, n
    while x != 1:
        c += x & 1
        x = x >> 1 if x % 2 == 0 else 3 * x + 1
    return c


O = np.array([oddcnt(n) for n in range(1, N + 1)], float)
for name, v in (("σ mod 2", par), ("O mod 2", O % 2)):
    Wv = fwht(v - v.mean())
    nsup = int(np.count_nonzero(np.abs(Wv) > 1e-6 * np.max(np.abs(Wv))))
    top = 32
    share = float(np.sort(np.abs(Wv))[::-1][:top].sum() / np.abs(Wv).sum())
    print(f"  {name}: WHT支撑={nsup}/{N} top-32占比={share:.4f}")

# σ mod 2 恒等式: σ(2m)=σ(m)+1 => 偶指标翻转; 自相似 => 2-regular 候选
flip_ok = all(sg(2 * m) % 2 != sg(m) % 2 for m in range(1, N // 2 + 1))
ok("C2a σ(2m)≡σ(m)+1 (mod 2) 全窗精确", flip_ok)

# --- C3 3-adic: 类均值极差 —— 否证 p-adic 连续性 (新发现) ---
M3 = 200000
sl3 = np.array([sg(n) for n in range(1, M3 + 1)], float)
print("\nC3 σ 按 n ≡ a (mod 3^j) 的类均值极差 (n≤200000):")
decay, amaxs = [], []
for j in range(1, 6):
    p = 3 ** j
    means = np.array([sl3[np.arange(M3) % p == a].mean() for a in range(p)])
    rng = float(means.max() - means.min())
    decay.append(rng)
    amaxs.append(int(means.argmax()))
    print(f"  3^{j}: 极差={rng:.3f} argmax类={amaxs[-1]}"
          f" (3^j-1={p - 1})")
slp3 = float(np.polyfit(range(1, 6), np.log(decay), 1)[0])
sd = float(sl3.std())
null5 = sd * math.sqrt(2 / (M3 / 3 ** 5))  # 3^5 类均值差的抽样噪声尺度
ok("C3 否证: 类均值极差随 3^j 增长 => σ 非 3-adic 连续",
   slp3 > 0.5 and decay[-1] > 3 * null5,
   f"斜率={slp3:.3f} (γ={-slp3:.3f}<0); 3^5极差={decay[-1]:.2f} vs "
   f"抽样噪声~{null5:.2f}; argmax类漂移={amaxs} (无固定极限点)")

# --- C4 可证回归: LTE v2(3^K-1)=2+v2(K) (K=16384=2^14 => 16) ---
v = ((3 ** 16384 - 1) & -(3 ** 16384 - 1)).bit_length() - 1
ok("C4 LTE v2(3^16384-1)=2+14=16", v == 16, f"v2={v}")

# ---------- 存档 ----------
out = {
    "note": "Round4: WHT+WPD+p-adic, 16384 dyadic 窗, 信号=results2 已验证复用",
    "wht_compare": comp,
    "wpd_profile": prof,
    "odd_chain": oddS,
    "c3_decay": decay,
    "c3_slope": slp3,
}
json.dump(out, open(r"C:\Users\savyc\collatz_fft_study\results3.json", "w"),
          indent=1)
print("\nresults3.json 已写入")
