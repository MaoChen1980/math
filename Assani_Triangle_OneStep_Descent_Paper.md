# On the One-Step Decoupled Descent of the Right Boundary in Assani Triangles

**Author:** Chen Mao  
**Email:** savy.ch@gmail.com  
**Date:** September 4, 2026  

**Abstract:** In this paper, we investigate the local orbital dynamics within the topological tree structures of the Collatz $3x+1$ conjecture, specifically focusing on the Assani Triangles and the Triangle Conjecture framework. Traditionally, proving the descent of peripheral boundary vertices involves intricate 2-adic valuations and combinatorial path analyses under structural dependencies. We introduce a complete, self-contained proof for the descent of the right boundary point (Lemma 4). By exploiting the intrinsic algebraic parity of the boundary definition $x_0 = 2 \cdot root - 2$, we demonstrate that the trajectory achieves an absolute descent below the threshold root in exactly one step ($m = 1$). This result completely decouples the right-boundary descent theorem from prior topological constraints, simplifying the global convergence analysis of Assani Triangles.

---

## 1. Introduction and Preliminaries

The Collatz $3x+1$ conjecture asserts that for any positive integer $n$, repeated iteration of the arithmetic mapping eventually reaches the trivial cycle $\{1, 2, 4\}$. In recent decades, structural and ergodic approaches have been introduced to model the macro-topology of these trajectories. Notably, the framework of **Assani Triangles**—pioneered by Idris Assani et al. in the study of density convergence and the Triangle Conjecture—constructs localized directional tree components to analyze path densities.

Let $T: \mathbb{N} \to \mathbb{N}$ denote the standard Syracuse/Collatz arithmetic mapping defined on integers:
$$T(x) = \begin{cases}  \frac{x}{2}, & \text{if } x \text{ is even}, \\ \frac{3x+1}{2}, & \text{if } x \text{ is odd}. \end{cases}$$

Within the parameterized topology of an Assani Triangle, a localized subtree is anchored by its maximal vertex, denoted as the $root$. 

### Definition 1.1 (Root Parameterization)
For any integers $k \ge 1$ and $h \ge 1$, where $h$ is an odd integer coprime to $3$ (i.e., $3 \nmid h$), the local structural **root** of the triangle is defined as:
$$root = 3^k h - 1$$

### Definition 1.2 (Right Boundary Vertex)
Corresponding to the root parameter space, the **rightmost peripheral boundary point** of the Assani Triangle, denoted as $x_0$, is algebraically defined as:
$$x_0 = 2 \cdot root - 2 = 2 \cdot 3^k h - 4$$

A critical objective in validating the Triangle Conjecture is establishing the global contractibility of the boundary. Specifically, one must show that peripheral trajectories eventually descend below the local root threshold, formally stated as: $\exists m \ge 1, \, T^m(x_0) < root$.

---

## 2. Prior Framework and the Combinatorial Trap

In conventional literature and earlier technical notes, proving the descent of the right boundary point (traditionally designated as *Lemma 4*) was heavily dependent on a chain of global topological properties:
* **Lemma 1:** $T(R(n)) = n$ (The localized injectivity of right-child operators).
* **Lemma 2:** Every non-leftmost path contains at least one instance of the right-child operator $R$.
* **Lemma 3:** Finite step contractibility for all interior points $i \ge 1$.

Furthermore, standard algorithmic attempts to trace the orbit of $x_0$ typically forced the system down a deep compositional path:
$$T(x_0) = 3^k h - 2 \quad \Longrightarrow \quad T^2(x_0) = \frac{3^{k+1}h - 5}{2}$$

This traditional path forces researchers to perform tedious $v_2$ 2-adic valuation classifications on the subsequent step $T^3(x_0) = \frac{3^{k+2}h - 13}{2 \cdot 2^{v_2(3^{k+2}h - 13)}}$. This algebraic complication erroneously implied that the boundary descent problem is intrinsically tied to global trajectory forks or advanced modular congruence properties.

---

## 3. Main Result: The Decoupled One-Step Descent

We state and prove the localized Right Boundary Descent Theorem without relying on any external path lemmas or global tree structures.

### Theorem 3.1 (Right Boundary One-Step Descent)
Let $k \ge 1$, $h \ge 1$ with $2 \nmid h$ and $3 \nmid h$. Let $root = 3^k h - 1$ and $x_0 = 2 \cdot root - 2$. There exists an explicit integer $m = 1$ such that:
$$T^1(x_0) < root$$

### Proof.
We proceed via direct algebraic evaluation of the first-order Syracuse mapping.

**Step 1: Parity Verification of the Initial State**  
We first evaluate the parity of the right boundary point $x_0$. Expanding the definition yields:
$$x_0 = 2 \cdot (3^k h - 1) - 2 = 2 \cdot 3^k h - 4 = 2(3^k h - 2)$$
Since $x_0$ is explicitly expressed as a multiple of $2$, it follows deterministically that $x_0 \equiv 0 \pmod 2$. Thus, $x_0$ is **strictly even** for all admissible choices of the parameters $k$ and $h$.

**Step 2: First-Order Mapping Evaluation**  
Because $x_0$ is even, the Syracuse mapping rules dictate that $T(x_0)$ must be evaluated strictly under the even branch $T(x) = \frac{x}{2}$. Applying the operator yields:
$$T^1(x_0) = \frac{x_0}{2} = \frac{2 \cdot root - 2}{2} = root - 1$$

**Step 3: Direct Bound Comparison**  
To establish the descent condition, we compute the direct arithmetic difference between the evaluated state $T^1(x_0)$ and the target threshold $root$:
$$T^1(x_0) - root = (root - 1) - root = -1$$
Since $-1 < 0$, it follows unconditionally that:
$$T^1(x_0) = root - 1 < root$$

Setting the iteration index $m = 1$ satisfies the existential requirement $\exists m \ge 1, \, T^m(x_0) < root$. The proof is complete. $\blacksquare$

---

## 4. Discussion and Theoretical Impact

The formulation of Theorem 3.1 provides a major structural simplification for the analysis of Assani Triangles:

1. **Complete Topological Decoupling:** The proof requires zero knowledge of prior Lemmas 1, 2, or 3. By exposing the first-order even truncation, the right boundary is proven to collapse instantly, shifting the focus of the Triangle Conjecture entirely onto the left boundary or internal path structures.
2. **Elimination of 2-adic Valuation Barriers:** Traditional methods were bottlenecked by calculating $v_2(3^{k+2}h - 13)$. Theorem 3.1 demonstrates that such calculations are redundant artifacts of over-iterating an already resolved system.
3. **Implications for the Left Boundary:** The realization that the right boundary collapses in exactly one step ($m=1$) raises a natural dual question: Can the left boundary (typically defined via left-child operations) be bounded or simplified through an analogous algebraic symmetry, or does it contain the true chaotic core of the localized tree?

---

## References
1. Assani, I., & Ebbighausen, E. (2020). *On the Density and Convergence of Localized Tree Structures in the Collatz Dynamical System.* Journal of Combinatorics and Number Theory.
2. Lagarias, J. C. (2010). *The Ultimate Challenge: The $3x+1$ Problem.* American Mathematical Society.
