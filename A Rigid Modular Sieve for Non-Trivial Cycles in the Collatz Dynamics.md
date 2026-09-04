---
title: "A Rigid Modular Sieve for Non-Trivial Cycles in the Collatz Dynamics"
author: "Chen Mao"
email: "savy.ch@gmail.com"
date: "September 5, 2026"
abstract: "We present a novel framework for excluding non-trivial cycles in the Collatz (3n+1) problem by coupling global algebraic constraints with local finite field structures. By decomposing the loop accumulator R into a structural minimum R_min(O) and an orbit-dependent discrepancy Δ(O), we project the cyclic trajectory onto a rigid congruence lattice over F_p. This yields a sieve that systematically eliminates parameter pairs (O, S) based on modular incompatibility, independent of numerical bounds. The method bridges the traditional gap between real-valued logarithmic estimates and 2-adic local analyses, offering a new arithmetic invariant for cycle exclusion."
---

## 1. Introduction

The Collatz map T: N → N is defined for odd integers as T(n) = (3n+1)/2^{v_2(3n+1)}. Despite decades of study, the existence of non-trivial cycles remains open. Traditional exclusion methods fall into two disconnected categories:

- **Real-domain bounds** (Lagarias-Weiss 1992, Eliahou 1993): use logarithmic estimates to show that any cycle must have enormous length (O > 10^{20}).
- **2-adic local constraints** (Kontorovich-Lagarias 2009): analyze the map on Z_2 but suffer from the Haar-counting measure singularity (M2).

These approaches are arithmetic-incompatible: real bounds say nothing about modular residues, and 2-adic constraints say nothing about global cycle length. This paper introduces P34, a **rigid modular sieve** that bridges this gap by projecting the global cycle identity onto a finite field F_p while retaining the full algebraic structure of the orbit.

## 2. The Classical Cycle Identity

For a cycle of odd length O and total step length S (including all divisions by 2), the starting point n_0 satisfies the exact algebraic identity:

$$
2^S n_0 = 3^O n_0 + R
$$

where

$$
R = \sum_{j=0}^{O-1} 3^j \cdot 2^{e_j}, \quad 0 = e_0 < e_1 < \cdots < e_{O-1} < S
$$

and e_j is the cumulative number of divisions by 2 up to the j-th odd step. From this, the divisibility constraint (XD3) follows directly:

$$
n_0 \mid (2^S - 3^O)
$$

and equivalently:

$$
n_0(2^S - 3^O) = R
$$

This is the foundation of all subsequent sieve construction.

## 3. Decomposition of the Accumulator

We introduce a canonical decomposition of R into a structural minimum and an orbit-dependent discrepancy.

Define:

$$
R_{\min}(O) = \frac{6^O - 1}{5}
$$

which is the algebraic lower bound obtained by setting e_j = j (the smallest possible increasing sequence). This bound is established in D12-R1. Although not dynamically attainable for all O, it serves as a universal lower bound.

Define the discrepancy:

$$
\Delta(O) = R - R_{\min}(O) = \sum_{j=0}^{O-1} 3^j (2^{e_j} - 2^j)
$$

Since e_j ≥ j, we have Δ(O) ≥ 0. The cyclic identity becomes:

$$
n_0(2^S - 3^O) = R_{\min}(O) + \Delta(O)
$$

## 4. The Modular Sieve Projector

Projecting onto the finite field F_p for an arbitrary prime p, any valid cycle configuration must satisfy:

$$
n_0 \cdot (2^S - 3^O) \equiv R_{\min}(O) + \Delta(O) \pmod{p}
$$

with the global boundary constraint:

$$
n_0 \le \frac{2^S - 3^O}{R_{\min}(O) + \Delta(O)}
$$

When Δ(O) ≡ 0 (mod p), this collapses to the pure sieve:

$$
n_0 \cdot (2^S - 3^O) \equiv R_{\min}(O) \pmod{p}
$$

## 5. The Cycle-Exclusion Criterion (P34)

**Theorem (P34, Cyclic Congruence Sieve):** If a non-trivial cycle exists, then for every prime p, there exist integers n_0, O, S satisfying:

1. The congruence projection:

$$
n_0(2^S - 3^O) \equiv R_{\min}(O) + \Delta(O) \pmod{p}
$$

2. The inequality constraint:

$$
n_0 \le \frac{2^S - 3^O}{R_{\min}(O) + \Delta(O)}
$$

3. The lower bound:

$$
S/O > \log_2 3
$$

from P3, and the upper bound:

$$
S/O < \log_2 3 + 1/O
$$

from XF3.

**Corollary:** Since R_min(O) mod p is purely periodic with period ord_p(6) (D26), for fixed O and p, the set of possible S is restricted to a rigid congruence lattice. Any (O, S) pair falling outside this lattice is algebraically impossible as a cycle candidate.

## 6. Sieve Efficacy and Computational Implications

The sieve can be applied without numerical iteration over n_0:

- For each O, compute R_min(O) mod p.
- For each candidate S within the logarithmic bounds (P3/XF3), test whether there exists n_0 ∈ N satisfying both congruence and inequality.
- If no n_0 exists, eliminate that (O, S) pair.

This is fundamentally different from traditional methods:

| Method | Domain | Constraint Type | Dependency |
|--------|--------|-----------------|------------|
| Lagarias-Weiss | Real logarithms | Size bounds | Numerical |
| Eliahou | Modular arithmetic | Divisibility | Cycle length only |
| P34 (this work) | Finite field F_p | Congruence + inequality | Algebraic, parameter-free |

## 7. Connection to the Discrepancy Distribution (M7)

The sieve's strength depends on the distribution of Δ(O) mod p. Under the natural density assumption (M7), Δ(O) mod p becomes uniformly distributed as O → ∞. This implies that the fraction of (O, S) pairs passing the sieve decays as O(p^{-m}) when multiple primes are used, providing a probabilistic exclusion of large cycles.

This bypasses the Haar-counting singularity (M2) by working directly on integer counting measure, rather than attempting to pull back from Z_2.

## 8. Conclusion and Outlook

P34 provides the first algebraic sieve that:

1. Uses exact cycle identities, not asymptotic estimates;
2. Operates in finite fields, leveraging D26's periodic structure;
3. Couples global cycle parameters (O, S) with local arithmetic invariants (R_min mod p);
4. Is parameter-free and numerically independent.

Future work includes:
- Quantifying the sieve's density for explicit prime sets;
- Extending to the full discrepancy term Δ(O) without the uniform distribution assumption;
- Implementing the sieve for O up to 10^6 to compare with known numerical bounds.

---

## References

1. Conway, J. H. (1972). Unpredictable iterations. *Proc. 1972 Number Theory Conf.*
2. Eliahou, S. (1993). The 3x+1 problem: new lower bounds for nontrivial cycle lengths. *Discrete Math.*
3. Kontorovich, A. V., & Lagarias, J. C. (2009). Odds and ends of the Collatz problem. *arXiv:0910.1944.*
4. Lagarias, J. C. (1985). The 3x+1 problem and its generalizations. *Amer. Math. Monthly.*
5. Lagarias, J. C., & Weiss, A. (1992). The 3x+1 problem: two stochastic models. *Ann. Appl. Probab.*