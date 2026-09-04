# -*- coding: utf-8 -*-
"""
Collatz FFT Sampling Study — Round 2: N=20000
20 采样方式 (沿用第一轮定义, 含大 N 适配) + N=16384 dyadic 对照窗
输出: results2.json + summary2.txt
"""
import json
import math
import random
import time
from collections import deque

import numpy as np

LOG2_3 = math.log2(3.0)
N_SAMPLES = 20000
N_DYAD = 16384  # 2^14 对照窗

_t0 = time.time()


def log(msg):
    print(f"[{time.time()-_t0:6.1f}s] {msg}", flush=True)


# ---------------- 核心 Collatz 库 ----------------

_memo = {1: 0}


def sigma_memo(n):
    path = []
    x = n
    while x not in _memo:
        path.append(x)
        x = x >> 1 if (x & 1) == 0 else 3 * x + 1
    s = _memo[x]
    for y in reversed(path):
        s += 1
        _memo[y] = s
    return _memo[n]


def sigma_direct(n):
    s, x = 0, n
    while x != 1:
        x = x >> 1 if (x & 1) == 0 else 3 * x + 1
        s += 1
    return s


def next_x(x):
    return x >> 1 if (x & 1) == 0 else 3 * x + 1


def stopping_orbit(n):
    """Terras 停时轨道: 到首个 < n 的值为止."""
    if n == 1:
        return [1]
    seq, x = [n], n
    while True:
        x = next_x(x)
        if x < n:
            break
        seq.append(x)
    return seq


def log2_big(x):
    if x <= 0:
        return 0.0
    bl = x.bit_length()
    return (bl - 1) + math.log2(x / (1 << (bl - 1)))


# ---------------- 20 种采样方式 (N=20000) ----------------

def m01_uniform_random():
    random.seed(42)
    ns = [random.randint(1, 2_000_000) for _ in range(N_SAMPLES)]
    return [float(sigma_direct(n)) for n in ns], {"n_range": "U[1,2e6] seed=42"}


def m02_lattice_odd_count():
    return [float(_odd_count(n)) for n in range(1, N_SAMPLES + 1)], \
        {"identity": "O(2m)=O(m)"}


def _odd_count(n):
    x, c = n, 0
    while True:
        if x & 1:
            c += 1
        if x == 1:
            break
        x = next_x(x)
    return c


def m03_mod16_stratified():
    sig = []
    for r in [1, 3, 5, 7, 9, 11, 13, 15]:
        for j in range(2500):
            sig.append(float(sigma_memo(r + 16 * j)))
    return sig, {"classes": 8, "per_class": 2500, "n_max": 39999}


def m04_mod3_condition():
    return [float(sigma_memo(1 + 3 * j)) for j in range(N_SAMPLES)], \
        {"cond": "n=1 mod 3", "n_max": 59998}


def m05_inverse_tree_bfs():
    q = deque([1])
    seen = {1}
    nodes = []
    while q and len(nodes) < N_SAMPLES:
        m = q.popleft()
        nodes.append(m)
        preds = [2 * m]
        if m % 6 == 4:
            preds.append((m - 1) // 3)
        for p in preds:
            if p not in seen:
                seen.add(p)
                q.append(p)
    return [float(sigma_memo(m)) for m in nodes], {"root": 1, "nodes": len(nodes)}


def m06_pow2_plusminus():
    """适配: k<=1250, 交错 -> 2500 样本 (大整数轨道成本截断)."""
    K = 1250
    sig = []
    for k in range(1, K + 1):
        sig.append(float(sigma_direct(2 ** k - 1)))
        sig.append(float(sigma_direct(2 ** k + 1)))
    return sig, {"k": f"1..{K}", "adapted": "大整数轨道成本截断 2500/20000"}


def m07_extreme_sigma():
    top = 20000
    dom = 200000
    sigmas = [sigma_direct(n) for n in range(1, dom + 1)]
    idx = sorted(np.argsort(sigmas)[::-1][:top].tolist())
    records, best = [], -1
    for n, s in enumerate(sigmas, 1):
        if s > best:
            best = s
            records.append((n, s))
    return [float(sigmas[i]) for i in idx], {
        "top_range": f"n<={dom}", "n_records": len(records),
        "last_records": records[-5:]}


def m08_level_sets():
    """适配: 域 n<=100000; bins k=1..20000, sigma 支撑 ~[1,350] 其余 0."""
    dom = 100000
    counts = np.zeros(N_SAMPLES)
    for n in range(1, dom + 1):
        s = sigma_direct(n)
        if 1 <= s <= N_SAMPLES:  # s=0 (n=1) 不属于任何 k>=1 档
            counts[s - 1] += 1
    return counts.tolist(), {"domain": f"n<={dom}",
                             "support": "k in [1,~350], 其余 0 填充"}


def m09_primes_composites():
    limit = 110000
    sieve = bytearray([1]) * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    primes = []
    for i in range(2, limit + 1):
        if sieve[i]:
            primes.append(i)
            sieve[i * i:: i] = b"\x00" * len(sieve[i * i:: i])
        if len(primes) >= 10000:
            break
    sig = [float(sigma_memo(p)) for p in primes[:10000]]
    sig += [float(sigma_memo(4 * i)) for i in range(1, 10001)]
    return sig, {"block1": "前10000素数", "block2": "4,8,...,40000"}


def m10_trivial_ramp():
    return [float(k) for k in range(1, N_SAMPLES + 1)], {
        "identity": "sigma(2^k)=k (可证 ramp)"}


def _stream(min_len):
    """Terras 停时轨道拼接流: n=2,3,4,... (排除共享尾部)."""
    vals, odd = [], 0
    n = 2
    while len(vals) < min_len:
        for v in stopping_orbit(n):
            vals.append(v)
            odd += v & 1
        n += 1
    return vals, odd, n - 1


def m11_parity_stream():
    vals, _, n_used = _stream(N_SAMPLES)
    return [1.0 if v & 1 else -1.0 for v in vals[:N_SAMPLES]], {
        "stream": f"停时轨道拼接 n=2..{n_used}",
        "adapted": "单轨道 sigma>=20000 不存在 (需 ~3000 位整数)"}


def m12_v2_sequence(N=N_SAMPLES):
    sig = []
    n = 1
    while len(sig) < N:
        m = 3 * n + 1
        sig.append(float((m & -m).bit_length() - 1))
        n += 2
    return sig, {"domain": f"奇 n=1..{2*N-1}", "model": "P2 几何分布"}


def m13_growth_phase_ramp():
    """n=2^10000-1 已证增长相: T^{2m}=3^m*2^{K-m}-1, T^{2m+1}=3^{m+1}*2^{K-m}-2
    (无条件代数: (a*2^j-1) 经奇+偶两步变 (3a, j-1); 奇步记录 3n+1 本身).
    K=10000 覆盖前 20000 步."""
    K = 10000
    sig = []
    p3, p2 = 1, 1 << K  # 3^m, 2^(K-m)
    for j in range(N_SAMPLES):
        m = j >> 1
        if j % 2 == 0:
            x = p3 * p2 - 1
        else:
            x = 3 * p3 * p2 - 2
        sig.append(log2_big(x))
        if j % 2 == 1:  # 完成一对, 进入下一个 m
            p3 *= 3
            p2 >>= 1
    return sig, {"start": "2^10000-1", "phase": "已证增长相 (斜率=(log2(3)-1)/2)"}


def m14_so_ratio_stream():
    vals, _, n_used = _stream(N_SAMPLES)
    sig, odd = [], 0
    for k in range(1, len(vals) + 1):
        odd += vals[k - 1] & 1
        if odd >= 1:
            sig.append(k / odd)
        if len(sig) == N_SAMPLES:
            break
    return sig, {"stream": f"n=2..{n_used}", "ref": "XF3"}


def m15_imbalance_stream():
    vals, _, n_used = _stream(N_SAMPLES)
    sig, odd = [], 0
    for k in range(1, N_SAMPLES + 1):
        odd += vals[k - 1] & 1
        sig.append(k - LOG2_3 * odd)
    return sig, {"stream": f"n=2..{n_used}",
                 "model": "漂移随机游走 -> 预测 1/f^2 谱"}


def m16_sigma_function(N=N_SAMPLES):
    return [float(sigma_memo(n)) for n in range(1, N + 1)], {
        "identity": "sigma(2m)=sigma(m)+1"}


def m17_rmin_mod_p():
    def gen(p):
        r, out = 0, []
        for _ in range(N_SAMPLES):
            r = (6 * r + 1) % p
            out.append(float(r))
        return out
    return gen(7), gen(11), {"periods": "ord_7(6)=2, ord_11(6)=10 (D26)"}


def m18_map_ratio():
    sig = []
    for n in range(1, N_SAMPLES + 1):
        if n % 2 == 0:
            t = n >> 1
        else:
            m = 3 * n + 1
            t = m >> ((m & -m).bit_length() - 1)
        sig.append(log2_big(t) - log2_big(n))
    return sig, {"even": -1.0, "odd": "log2(3)-v2(+o(1))"}


def m19_cross_section():
    sig = []
    for i in range(1, N_SAMPLES + 1):
        x = 100 * i
        for _ in range(60):
            if x == 1:
                break
            x = next_x(x)
        sig.append(log2_big(x))
    return sig, {"fixed_k": 60, "n": "100,200,...,2e6"}


def m20_pure_L_branches():
    nodes, seen = [], set()
    r = 2
    while len(nodes) < N_SAMPLES:
        m = r
        while True:
            if m in seen:
                break
            seen.add(m)
            nodes.append(m)
            if (2 * m - 1) % 3 != 0:
                break
            m = (2 * m - 1) // 3
        r += 3
    return [float(sigma_memo(m)) for m in nodes[:N_SAMPLES]], {
        "roots": "m=2 mod 3 递增", "O3": "纯L分支"}


# ---------------- FFT ----------------

def analyze(sig):
    x = np.asarray(sig, dtype=float)
    n = x.size
    t = np.arange(n)
    slope, icpt = np.polyfit(t, x, 1)
    resid = x - (slope * t + icpt)
    mag = np.abs(np.fft.rfft(resid))
    if mag.size > 1 and mag[1:].max() > 0:
        order = np.argsort(mag[1:])[::-1][:5] + 1
        top = [(int(k), round(k / n, 6),
                round(float(mag[k] / mag[order[0]]), 3)) for k in order]
        p = mag[1:] ** 2 + 1e-300
        flat = float(np.exp(np.mean(np.log(p))) / np.mean(p))
    else:
        top, flat = [], 0.0
    even_frac = float(mag[2::2].__pow__(2).sum() / (mag[1:] ** 2).sum()) \
        if mag.size > 2 else 0.0
    return {
        "mean": round(float(x.mean()), 4), "std": round(float(x.std()), 4),
        "slope": round(float(slope), 6),
        "top_peaks": [{"k": k, "freq": f, "rel": r} for k, f, r in top],
        "flatness": round(flat, 5),
        "even_k_energy": round(even_frac, 5),
    }


METHODS = [
    ("M01", "均匀随机 sigma", m01_uniform_random),
    ("M02", "格点 O(n)", m02_lattice_odd_count),
    ("M03", "模16分层", m03_mod16_stratified),
    ("M04", "模3条件", m04_mod3_condition),
    ("M05", "逆树BFS", m05_inverse_tree_bfs),
    ("M06", "2^k±1 [2500适配]", m06_pow2_plusminus),
    ("M07", "极端sigma top20000", m07_extreme_sigma),
    ("M08", "水平集 [0填充适配]", m08_level_sets),
    ("M09", "素数/合数 x10000", m09_primes_composites),
    ("M10", "平凡环 ramp", m10_trivial_ramp),
    ("M11", "奇偶±1拼接流", m11_parity_stream),
    ("M12", "v2 序列", m12_v2_sequence),
    ("M13", "2^k-1 已证增长相", m13_growth_phase_ramp),
    ("M14", "S/O 累积比拼接流", m14_so_ratio_stream),
    ("M15", "失衡拼接流(1/f²?)", m15_imbalance_stream),
    ("M16", "sigma(n) 函数[主选]", m16_sigma_function),
    ("M17", "R_min mod 7/11", m17_rmin_mod_p),
    ("M18", "log2(T(n)/n)", m18_map_ratio),
    ("M19", "横截面 T^60", m19_cross_section),
    ("M20", "纯L分支 O3", m20_pure_L_branches),
]


def main():
    results = {}
    for mid, desc, fn in METHODS:
        out = fn()
        entries = []
        if mid == "M17":
            s7, s11, extra = out
            for name, s, per in (("M17a R_min mod 7", s7, 2),
                                 ("M17b R_min mod 11", s11, 10)):
                a = analyze(s)
                a["theory_period"] = per
                entries.append({"name": name, "signal": s, "fft": a})
        elif mid == "M12":
            s, extra = out
            sd, extrad = m12_v2_sequence(N_DYAD)
            entries.append({"name": "M12 v2 N=20000", "signal": s,
                            **extra, "fft": analyze(s)})
            entries.append({"name": "M12d v2 N=16384(dyadic窗)",
                            "signal": sd, **extrad, "fft": analyze(sd)})
        elif mid == "M16":
            s, extra = out
            sd, extrad = m16_sigma_function(N_DYAD)
            entries.append({"name": "M16 sigma N=20000", "signal": s,
                            **extra, "fft": analyze(s)})
            entries.append({"name": "M16d sigma N=16384(dyadic窗)",
                            "signal": sd, **extrad, "fft": analyze(sd)})
        else:
            s, extra = out
            entries.append({"name": mid, "signal": s, **extra,
                            "fft": analyze(s)})
        results[mid] = {"desc": desc, "entries": entries}
        log(f"{mid} {desc} ok")

    results["_meta"] = {"N": N_SAMPLES, "N_dyad": N_DYAD,
                        "log2_3": LOG2_3, "date": "2026-09-04",
                        "note": "Round2: 适配见各条目 adapted 字段"}

    with open(r"C:\Users\savyc\collatz_fft_study\results2.json", "w",
              encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, default=float)

    lines = [f"Collatz FFT Study Round 2 — N={N_SAMPLES} (对照窗 {N_DYAD})",
             "=" * 64, ""]
    for mid, desc, _ in METHODS:
        for e in results[mid]["entries"]:
            a = e["fft"]
            pk = " | ".join(f"k={p['k']}(f={p['freq']:.5f},{p['rel']*100:.0f}%)"
                            for p in a["top_peaks"][:4])
            lines.append(f"{e['name']}  [{results[mid]['desc']}]")
            lines.append(f"  mean={a['mean']} std={a['std']} slope={a['slope']} "
                         f"flatness={a['flatness']} even_k_E={a['even_k_energy']}")
            lines.append(f"  peaks: {pk}")
            for k in e:
                if k not in ("name", "signal", "fft"):
                    lines.append(f"  {k}: {e[k]}")
            lines.append("")
    with open(r"C:\Users\savyc\collatz_fft_study\summary2.txt", "w",
              encoding="utf-8") as f:
        f.write("\n".join(lines))
    log("DONE -> summary2.txt / results2.json")


if __name__ == "__main__":
    main()
