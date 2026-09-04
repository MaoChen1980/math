---
title: "Correction and Conditional Reformulation of the Mod-16 Expansion Bound in Collatz Dynamics"
author: "Chen Mao"
email: "savy.ch@gmail.com"
date: "September 5, 2026"
abstract: "We correct a significant arithmetic error in the previously claimed absolute upper bound τ(n) ≤ C log n for Collatz trajectories. The error arises from counting four-step blocks as two-step blocks, leading to an undercount of the worst-case amplification. We show that the correct four-step net amplification is (3/2)^3 × (1/2) = 27/16 > 1. Consequently, the constant-C absolute bound is invalid. We reformulate P36 as a conditional bound that depends on the choice of the deterministic window length K, and we provide the corrected framework for future analysis of state-machine expansion envelopes."
---

## 1. The Original Claim

In earlier versions of this work, the following proposition was stated as rigorously proven:

**Original P36 (erroneous):** There exists a constant C = 4 / log₂(4/3) ≈ 9.64 such that for all n, τ(n) ≤ C log₂ n, where τ(n) is the number of steps until the trajectory first falls below n.

**Proof attempt (flawed):** Since the mod-16 deterministic transition matrix guarantees at least one division by 2 in every 4-step block, and the maximum single-step amplification is 3/2, the net amplification over 4 steps is at most:

$$
\frac{3}{2} \times \frac{1}{2} = \frac{3}{4} < 1
$$

This is **arithmetically incorrect**.

## 2. The Correction

A 4-step block consists of **four individual transformations**. If the worst case is three odd steps (amplification by approximately 3/2 each) and one even step (division by 2), the net amplification is:

$$
\left(\frac{3}{2}\right)^3 \times \frac{1}{2} = \frac{27}{16} = 1.6875 > 1
$$

The original calculation treated the 4-step block as if it had only two steps: one amplification and one division. This is a factor of 4 error in step counting.

## 3. Implications for the Constant-C Bound

Since 27/16 > 1, a 4-step block can have net amplification. Therefore:

- There is **no immediate guarantee** that the trajectory decreases after any fixed number of steps.
- The claimed constant \(C = 4/\log_2(4/3)\) is invalid.
- The proposition "there exists an absolute constant C such that τ(n) ≤ C log₂ n for all n" is **not established** by this argument.

## 4. Reformulation as a Conditional Bound (P36 Revised)

**Revised P36 (conditional):** For any chosen deterministic window length K (i.e., working modulo 2^K), the maximal possible amplification over K steps is bounded by the maximum product of local transition factors along the state machine's longest path with net growth > 1. Let G_K be the maximum over all K-step paths of:

$$
\prod_{i=1}^K \frac{3}{2^{v_i}} \quad \text{where } v_i \ge 1
$$

If G_K < 1 for some K, then:

$$
\tau(n) \le K \cdot \lceil \log_{1/G_K} n \rceil = \frac{K}{\log_2(1/G_K)} \log_2 n
$$

**Status:** Conditional upon determining G_K for a finite K. For K=4, G_4 = 27/16 > 1, so this K is insufficient.

**Future direction:** Compute G_K for K=5, 6, ... (i.e., mod 32, mod 64, etc.) until the first K with G_K < 1 appears. If such a K exists, an absolute bound can be recovered. If no such K exists, then no constant-C absolute bound is possible from finite-state arguments alone.

## 5. Connection to Broader Lessons

This correction is significant for three reasons:

1. **Methodological honesty:** It demonstrates that the framework is self-correcting and falsifiable. A flawed constant was identified, corrected, and the proposition was downgraded rather than hidden.

2. **Boundary awareness:** The m=1 boundary case in σ(2m−1) = 2 + σ(3m−1) was also corrected (σ(1)=0, not 3). Together, these corrections show that rigorous Collatz analysis requires meticulous attention to the **base cases** and **worst-case path enumeration**.

3. **State machine limits:** The correction reveals that the mod-16 state machine alone is insufficient to prove a global boundedness constant. Higher moduli (mod 2^K, K≥5) must be examined, or entirely different methods must be employed.

## 6. Conclusion

The corrected P36 no longer claims an absolute constant-C bound. Instead, it provides a **conditional framework**: if a sufficiently long deterministic window with net contraction G_K < 1 exists, then a bound follows; otherwise, finite-state methods are insufficient for global boundedness. This honest reformulation strengthens the overall framework by clarifying its assumptions and boundaries.

---

## References

1. Conway, J. H. (1972). Unpredictable iterations. *Proc. 1972 Number Theory Conf.*
2. Lagarias, J. C. (1985). The 3x+1 problem and its generalizations. *Amer. Math. Monthly.*
3. Tao, T. (2019). Almost all orbits of the Collatz map attain almost all values. *arXiv:1908.08104.*
4. Allouche, J.-P., & Shallit, J. (2003). Automatic Sequences. *Cambridge UP.*