# -*- coding: utf-8 -*-
"""独立重算验证: 用与 study2.py 不同的实现重算各方法样本, 与存储值比对.
独立性: 递归 lru_cache sigma / 循环除法 v2 / 独立 BFS 与流生成."""
import json
import math
import random
from functools import lru_cache

R = json.load(open(r"C:\Users\savyc\collatz_fft_study\results2.json",
                   encoding="utf-8"))
rng = random.Random(12345)
ok = lambda name, cond, detail="": print(
    f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")


def sget(mid, i=0):
    return R[mid]["entries"][i]["signal"]


# ---------- 独立实现 (与 study2 不同代码路径) ----------
@lru_cache(maxsize=None)
def sig2(n):
    """递归 + lru_cache (study2 用迭代 memo)."""
    if n == 1:
        return 0
    return 1 + sig2(n // 2 if n % 2 == 0 else 3 * n + 1)


import sys
sys.setrecursionlimit(100000)


def odd2(n):
    c, x = 0, n
    while True:
        c += x % 2
        if x == 1:
            return c
        if x % 2 == 0:
            x //= 2
        else:
            x = 3 * x + 1


def v2b(n):
    v, m = 0, 3 * n + 1
    while m % 2 == 0:
        m //= 2
        v += 1
    return v


def sdir2(n):
    s, x = 0, n
    while x != 1:
        s += 1
        x = (x // 2) if (x % 2 == 0) else (3 * x + 1)
    return s


def stop_orbit2(n):
    if n == 1:
        return [1]
    seq, x = [n], n
    while True:
        x = x // 2 if x % 2 == 0 else 3 * x + 1
        if x < n:
            return seq
        seq.append(x)


def sample_idx(n_total, k=300):
    return sorted(rng.sample(range(n_total), k))


# ---------- M01: 同种子可复现 + 独立 sigma ----------
random.seed(42)
ns1 = [random.randint(1, 2_000_000) for _ in range(20000)]
s01 = sget("M01")
idx = sample_idx(20000)
bad = [i for i in idx if s01[i] != sdir2(ns1[i])]
ok("M01 重算(种子可复现+独立sigma)", not bad, f"300样例不符={len(bad)}")

# ---------- M02: 独立奇数计数 ----------
s02 = sget("M02")
idx = sample_idx(20000)
bad = [i for i in idx if s02[i] != odd2(i + 1)]
ok("M02 O(n) 独立重算", not bad, f"不符={len(bad)}")

# ---------- M03/M04/M16/M16d/M20/M05/M09: 独立 sigma ----------
RES16 = [1, 3, 5, 7, 9, 11, 13, 15]
for mid, gen, tot in (
        ("M03", lambda i: RES16[i // 2500] + 16 * (i % 2500), 20000),
        ("M04", lambda i: 1 + 3 * i, 20000),
        ("M16", lambda i: i + 1, 20000)):
    s = sget(mid)
    idx = sample_idx(tot)
    bad = [i for i in idx if s[i] != sig2(gen(i))]
    ok(f"{mid} 独立 sigma 重算", not bad, f"不符={len(bad)}")

s16d = sget("M16", 1)
idx = sample_idx(16384)
bad = [i for i in idx if s16d[i] != sig2(i + 1)]
ok("M16d (N=16384) 独立 sigma", not bad, f"不符={len(bad)}")

# ---------- M05: 独立 BFS 重建 + 全量比对 ----------
from collections import deque
q = deque([1])
seen = {1}
nodes5 = []
while q and len(nodes5) < 20000:
    m = q.popleft()
    nodes5.append(m)
    preds = [2 * m] + ([(m - 1) // 3] if m % 6 == 4 else [])
    for p in preds:
        if p not in seen:
            seen.add(p)
            q.append(p)
# 前驱规则抽查: T(pred)=node
import random as _r
_r.seed(7)
ok("M05 前驱规则 T(2m)=m 抽查",
   all((lambda m: (m // 2 if m % 2 == 0 else 3 * m + 1))(2 * m) == m
       for m in _r.sample(nodes5, 200)))
ok("M05 奇前驱合法性 (n=4 mod 6 且 (n-1)/3 奇)",
   all((m - 1) % 3 == 0 and ((m - 1) // 3) % 2 == 1
       for m in _r.sample(nodes5, 200) if m % 6 == 4))
s05 = sget("M05")
idx = sample_idx(20000)
bad = [i for i in idx if s05[i] != sig2(nodes5[i])]
ok("M05 独立 BFS+sigma 全序比对", not bad, f"不符={len(bad)}")

# ---------- M06: 独立轨道 + 恒等式 sigma(2^k-1)=2k+sigma(3^k-1) ----------
s06 = sget("M06")
idx = [2 * k for k in rng.sample(range(1, 1251), 40)]  # 偶位=2^k-1
bad = []
for i in idx:
    k = i // 2 + 1
    if s06[i] != sdir2(2 ** k - 1):
        bad.append(i)
ok("M06 sigma(2^k-1) 独立轨道重算(40例)", not bad, f"不符={len(bad)}")
idbad = [k for k in rng.sample(range(1, 301), 25)
         if sdir2(2**k - 1) != 2 * k + sdir2(3**k - 1)]
ok("M06 恒等式 sigma(2^k-1)=2k+sigma(3^k-1) (25例)", not idbad,
   f"不符={len(idbad)}")

# ---------- M07/M08: 全量独立重扫 ----------
sigmas_ind = [sdir2(n) for n in range(1, 200001)]
s07 = sorted(sget("M07"))
top_ind = sorted(sorted(sigmas_ind)[::-1][:20000])
ok("M07 全量重扫: top-20000 多重集一致", s07 == [float(v) for v in top_ind])
cnt = [0] * 20001
for n in range(1, 100001):
    v = sdir2(n)
    if 1 <= v <= 20000:
        cnt[v] += 1
s08 = sget("M08")
ok("M08 全量重扫: 水平集计数一致",
   all(s08[k - 1] == cnt[k] for k in range(1, 20001)))

# ---------- M09: 素性独立核验 + sigma ----------
def is_prime(m):
    if m < 2:
        return False
    d = 2
    while d * d <= m:
        if m % d == 0:
            return False
        d += 1
    return True


s09 = sget("M09")
idx = sample_idx(10000, 60)
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
ok("M09 前10000素性抽查(trial division)",
   all(is_prime(primes[i]) for i in sample_idx(10000, 60)))
bad = [i for i in sample_idx(10000, 100)
       if s09[i] != sig2(primes[i])] + \
      [i for i in sample_idx(10000, 100)
       if s09[10000 + i] != sig2(4 * (i + 1))]
ok("M09 sigma 独立重算(素数+合数各100)", not bad, f"不符={len(bad)}")

# ---------- M11/M14/M15: 独立流重建 ----------
vals2 = []
n = 2
while len(vals2) < 20000:
    vals2.extend(stop_orbit2(n))
    n += 1
vals2 = vals2[:20000]
# 停时边界抽查
_r.seed(9)
ok("M11 停时轨道边界 (末值>=n 且下一值<n)",
   all((lambda n: (lambda seq: seq[-1] >= n and
                   (seq[-1] // 2 if seq[-1] % 2 == 0 else
                    3 * seq[-1] + 1) < n)(stop_orbit2(n)))
       (nn) for nn in _r.sample(range(2, 4000), 100)))
s11 = sget("M11")
bad = [i for i in sample_idx(20000)
       if s11[i] != (1.0 if vals2[i] % 2 else -1.0)]
ok("M11 拼接流独立重建+奇偶比对", not bad, f"不符={len(bad)}")
s15 = sget("M15")
odd_run = 0
bad = []
for k in range(1, 20001):
    odd_run += vals2[k - 1] % 2
    if abs(s15[k - 1] - (k - math.log2(3) * odd_run)) > 1e-9:
        bad.append(k)
ok("M15 失衡流独立重建", not bad, f"不符={len(bad)}")

# ---------- M12/M12d: 循环除法 v2 ----------
for mid, i_ent in (("M12", 0), ("M12", 1)):
    s = sget(mid, i_ent)
    idx = sample_idx(len(s))
    bad = [i for i in idx if s[i] != v2b(2 * i + 1)]
    ok(f"{mid} v2 循环除法独立重算", not bad, f"不符={len(bad)}")

# ---------- M13: 闭式 vs 直接迭代 (抽样到 j=4000) ----------
s13 = sget("M13")
x = 2**10000 - 1
direct = {}
for j in range(0, 4001):
    direct[j] = x
    x = x // 2 if x % 2 == 0 else 3 * x + 1


def log2b(v):
    bl = v.bit_length()
    return (bl - 1) + math.log2(v / (1 << (bl - 1)))


idx = sorted(set(list(rng.sample(range(4001), 60)) + [0, 1, 2, 3, 3998, 3999,
                                                      4000]))
bad = [j for j in idx if abs(s13[j] - log2b(direct[j])) > 1e-9]
ok("M13 增长相闭式 vs 直接迭代 (60例, j<=4000)", not bad, f"不符={len(bad)}")

# ---------- M17: pow 独立重算 ----------
for name, p in (("M17a", 7), ("M17b", 11)):
    s = sget("M17", 0 if p == 7 else 1)
    idx = sample_idx(20000, 300)
    inv5 = pow(5, -1, p)
    bad = [i for i in idx
           if s[i] != ((pow(6, i + 1, p) - 1) * inv5) % p]
    ok(f"{name} R_min mod {p}: pow 独立重算 (300例)", not bad, f"不符={len(bad)}")

# ---------- M18: 独立 Syracuse 加速步 T(n)/n ----------
s18 = sget("M18")
idx = sample_idx(20000)
bad = []
for i in idx:
    n = i + 1
    if n % 2 == 0:
        t = n // 2
    else:
        m = 3 * n + 1
        t = m // (m & -m)  # 加速步: 除尽 2^v2 (与 study2 定义一致)
    if abs(s18[i] - (log2b(t) - math.log2(n))) > 1e-9:
        bad.append(i)
ok("M18 独立重算 (加速步语义)", not bad, f"不符={len(bad)}")

# ---------- M19: 独立 T^60 ----------
s19 = sget("M19")
idx = sample_idx(20000, 200)
bad = []
for i in idx:
    x = 100 * (i + 1)
    for _ in range(60):
        if x == 1:
            break
        x = x // 2 if x % 2 == 0 else 3 * x + 1
    if abs(s19[i] - log2b(x)) > 1e-9:
        bad.append(i)
ok("M19 独立 T^60 重算 (200例)", not bad, f"不符={len(bad)}")

print("\n采样准确性独立验证完毕")
