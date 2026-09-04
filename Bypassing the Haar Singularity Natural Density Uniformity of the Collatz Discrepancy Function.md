---
title: "Bypassing the Haar Singularity: Natural Density Uniformity of the Collatz Discrepancy Function"
author: "Chen Mao"
email: "savy.ch@gmail.com"
date: "September 5, 2026"
abstract: "The Haar measure on Z_2, although preserved by the Syracuse map, is singular with respect to the counting measure on N. This prohibits direct transfer of 2-adic almost-everywhere results to integer orbit properties. We introduce M7, which circumvents this obstruction by working entirely in the integer counting measure (natural density). We define the discrepancy function Δ(O) = R - R_min(O) and conjecture that Δ(O) mod p becomes uniformly distributed as O → ∞. This provides a probabilistically valid channel for finite-field sieving that is independent of Haar measure, opening the door for probabilistic number theory methods in Collatz cycle exclusion."
---

## 1. The Haar Obstruction

A beautiful and well-established result in Collatz theory (Kontorovich-Lagarias 2009, M1) states that the Syracuse map T preserves the Haar measure μ_2 on Z_2:

$$
\mu_2(T^{-1}(A)) = \mu_2(A) \quad \forall A \subset Z_2
$$

This has led to powerful "almost everywhere" results: Tao (2019) proves boundedness for μ_2-almost all 2-adic integers, and Terras (1976) proves eventual descent for μ_2-almost all integers.

However, as established in M2 and M4:

$$
\mu_2 \perp \#_N
$$

where \(\#_N\) is the counting measure on N. The two measures are mutually singular. Consequently:

- μ_2-almost-everywhere does **not** imply natural-density-almost-everywhere in N.
- 2-adic spectral properties cannot be directly pulled back to integer frequency spectra.
- This is the fundamental reason why distributional results (T1/T2) do not imply pointwise convergence (T3).

## 2. The Natural Density Framework

Let A ⊂ N. The natural density (or asymptotic density) is defined as:

$$
d(A) = \lim_{N \to \infty} \frac{|A \cap \{1,\dots,N\}|}{N}
$$

if the limit exists. This is the appropriate measure for integer-orbit statements: "almost all n" in Collatz discourse typically means natural density 1.

Our goal is to formulate statements about Collatz dynamics directly in this measure, without reference to Z_2 or Haar measure.

## 3. The Discrepancy Function

Recall the cycle accumulator decomposition:

$$
R = R_{\min}(O) + \Delta(O)
$$

where:

$$
R_{\min}(O) = \frac{6^O - 1}{5}, \quad \Delta(O) = \sum_{j=0}^{O-1} 3^j (2^{e_j} - 2^j)
$$

The distribution of Δ(O) mod p is determined by the distribution of e_j mod ord_p(2). From P2, the e_j are geometrically distributed with P(v_2 = k) = 2^{-k} in the Haar sense. The key question: does this extend to natural density?

## 4. The Uniformity Conjecture (M7)

**Conjecture (M7, Δ-uniformity):** For any prime p, the sequence Δ(O) mod p is uniformly distributed on F_p under the natural counting density:

$$
\lim_{N \to \infty} \frac{1}{N} \sum_{O=1}^N \mathbf{1}_{\{\Delta(O) \equiv a \pmod{p}\}} = \frac{1}{p} \quad \forall a \in \mathbb{F}_p
$$

**Supporting heuristic:** Since P2 implies that v_2(3n+1) is geometrically distributed, we expect that e_j mod ord_p(2) behaves like independent random variables with uniform distribution on the cyclic group. Then Δ(O) mod p is a sum of O independent, identically distributed random variables, and the law of large numbers on finite groups gives uniformity as O → ∞.

**Technical requirement:** Rigorous proof requires non-trivial estimates of 2-adic exponential sums:

$$
\sum_{n \le N} \chi(n) e^{2\pi i a e_j(n) / \operatorname{ord}_p(2)}
$$

which may be approached via 2-adic harmonic analysis or mixed character sum bounds.

## 5. Application to Cycle Exclusion

If M7 holds, then for any fixed prime p and any target residue c ∈ F_p, the set of O for which:

$$
\Delta(O) \equiv c - R_{\min}(O) \pmod{p}
$$

has natural density 1/p. In particular, for a large ensemble of primes {p_1, ..., p_m}, the density of O satisfying the full sieve condition decays as p^{-m}.

Combining with P34, the number of (O, S) pairs passing the modular sieve for m primes is:

$$
O \cdot \left(\prod_{i=1}^m \frac{1}{p_i}\right) \cdot (\# \text{candidate S per O})
$$

Since XF3 constrains S to an interval of length O/log 2, the total number of surviving pairs grows sub-quadratically, and the natural density of surviving cycles tends to 0 as m → ∞.

## 6. Comparison with Haar-Based Approaches

| Approach | Domain | Valid Statements | Limitations |
|----------|--------|------------------|-------------|
| Haar measure on Z_2 | 2-adic integers | Almost-everywhere in μ_2 | Cannot pull back to N |
| Natural density (M7) | Integers N | Almost-everywhere in counting measure | Requires uniformity proof |
| Direct pointwise | Each n ∈ N | For all n | Unproven (Collatz itself) |

M7 occupies a middle ground: it is stronger than Haar results for integer applications, but weaker than full pointwise proof. If proven, it would be the first rigorous connection between 2-adic probabilistic structure and integer asymptotic density.

## 7. Conclusion

M7 offers a rigorous path to bypass the Haar-counting singularity: work directly in the integer counting measure from the outset. The conjecture that Δ(O) mod p is uniformly distributed is plausible from P2's geometric distribution and provides a probabilistic sieve for cycle exclusion. Its proof would open the door for probabilistic number theory methods (e.g., Erdős–Kac type theorems) to contribute to Collatz cycle analysis without invoking the non-transferable Haar measure.

---

## References

1. Kontorovich, A. V., & Lagarias, J. C. (2009). Odds and ends of the Collatz problem. *arXiv:0910.1944.*
2. Tao, T. (2019). Almost all orbits of the Collatz map attain almost all values. *arXiv:1908.08104.*
3. Terras, R. (1976). A stopping time problem on the positive integers. *Acta Arith.*
4. Tenenbaum, G. (1995). Introduction to Analytic and Probabilistic Number Theory. *Cambridge UP.*
5. Matthews, K. (2012). The 3x+1 problem: a probabilistic approach. *J. Integer Seq.*