# Erdős #40: notes

Erdős problem 40 is open, and nothing here solves it. These notes prove two reductions and one exclusion, reproduce one published computation, and say which of those claims a machine checked.

Each claim carries one of four labels:

* **[Lean]**: proved in `lean/`, about the definitions of formal-conjectures at the commit pinned in `lean/lakefile.toml`, with the checks of §8.
* **[written]**: a complete written proof in this file, not machine-checked.
* **[source]**: quoted from a source listed in §9, with the quotation checked mechanically against the fetched text.
* **[computation]**: a finite computation. No finite computation decides the questions here, so these are illustrative.

| § | contents | label |
|---|---|---|
| 1 | the statement, the formal statement, and how they correspond | [source], [written] |
| 2 | Theorem A: some $g\to\infty$ works if and only if $g=1$ works | [Lean] |
| 3 | the same question for $B_2[g]$ sets; problems 28 and 158 | [Lean], [source] |
| 4 | Theorem B: no $g\ge cN^{\epsilon}$ works (after Erdős and Rényi) | [written] |
| 5 | what is known about the answer set | all four |
| 6 | the greedy $B_2[2]$ computation from the problem 158 thread | [computation] |
| 7 | status and the blocking question | |
| 8 | how the claims were checked | |
| 9 | sources read and not read | |
| 10 | provenance | |

## 1. The problem

### 1.1 Statement

<!-- statement: site_latex_40.html -->
```latex
For what functions $g(N)\to \infty$ is it true that\[\lvert A\cap \{1,\ldots,N\}\rvert \gg \frac{N^{1/2}}{g(N)}\]implies $\limsup 1_A\ast 1_A(n)=\infty$?
```

This is the statement in the LaTeX view of [erdosproblems.com/40](https://www.erdosproblems.com/40), byte for byte (`refs/site_latex_40.html`). Throughout, $A\subseteq\mathbb{N}$, $A(N)=\lvert A\cap\{1,\ldots,N\}\rvert$, and $r_A(n)=1_A\ast 1_A(n)$ is the number of ordered pairs $(a,b)\in A^2$ with $a+b=n$.

The page, as fetched on September 23, 2026, gives the status as "This is open, and cannot be resolved with a finite computation." and offers a prize of \$500. Its only commentary is "This is a stronger form of the Erdős-Turán conjecture [28] (since establishing this for any function $g(N)\to \infty$ would imply a positive solution to [28])." Its discussion thread has no comments. Its sources are [Er95] and [Er97c], neither of which I read (§9).

Problem 28 is the Erdős–Turán conjecture: "If $A\subseteq \mathbb{N}$ is such that $A+A$ contains all but finitely many integers then $\limsup 1_A\ast 1_A(n)=\infty$." Its page adds: "Another stronger conjecture would be that the hypothesis $\lvert A\cap [1,N]\rvert \gg N^{1/2}$ for all large $N$ suffices." That is the case $g=1$ of problem 40, and these notes call it *the strong form*.

### 1.2 The formal statement

formal-conjectures states problem 40 in `FormalConjectures/ErdosProblems/40.lean`. At the pinned commit `e5f428182a3ee32dde401eceb4a94ba0382e8434`:

<!-- fc: fc_40.lean -->
```lean
def Erdos40For (g : ℕ → ℝ) : Prop :=
  ∀ A : Set ℕ,
    (fun N : ℕ ↦ √N / g N) =O[atTop] (fun N ↦ ((A ∩ .Icc 1 N).ncard : ℝ)) →
    limsup (fun N ↦ (sumRep A N : ℕ∞)) atTop = ⊤
```

<!-- fc: fc_40.lean -->
```lean
theorem erdos_40 :
    {g : ℕ → ℝ | Tendsto g atTop atTop ∧ Erdos40For g} = answer(sorry) := by
```

<!-- fc: fc_40.lean -->
```lean
theorem erdos_40.variants.weaker :
    answer(sorry) ↔ ∃ g : ℕ → ℝ, Tendsto g atTop atTop ∧ Erdos40For g := by
```

The representation count is defined in `FormalConjecturesForMathlib/Combinatorics/Additive/Convolution.lean`, and `sumRep_def` there unfolds it:

<!-- fc: fc_Convolution.lean -->
```lean
noncomputable def sumRep (A : Set ℕ) : ℕ → ℕ := (𝟙_A ∗ 𝟙_A)
```

<!-- fc: fc_Convolution.lean -->
```lean
    sumRep A n = ((antidiagonal n).filter (fun (p : ℕ × ℕ) ↦ p.1 ∈ A ∧ p.2 ∈ A)).card := by
```

Write $\mathcal{G}$ for the set in `erdos_40`, the answer set: the $g$ with $g(N)\to\infty$ that have the property `Erdos40For g`. The strong form is `Erdos40For (fun _ ↦ 1)`. The files for problems 28, 40 and 158 carry no `formal_proof` attribute, and the site's page for problem 40 says its statement is formalised and links this file on the main branch. On September 23, 2026 the four formal-conjectures files these notes rely on, `ErdosProblems/28.lean`, `40.lean`, `158.lean` and `Convolution.lean`, were byte-identical at the pin and on the main branch.

How the formal statement reads the informal one:

| informal | formal-conjectures | reading |
|---|---|---|
| $1_A\ast 1_A(n)$ | `sumRep A n` | the ordered pairs $(a,b)\in A^2$ with $a+b=n$, as in the convolution. A pair with $a\ne b$ counts twice and $(a,a)$ once, which cannot change whether the count is bounded |
| $\lvert A\cap\{1,\ldots,N\}\rvert$ | `(A ∩ .Icc 1 N).ncard` | the same finite set; `ncard` is its size |
| $\gg$ | `=O[atTop]` | there are $c>0$ and $N_0$ with $\sqrt N/g(N)\le c\,A(N)$ for all $N\ge N_0$. Requiring the bound for every $N\ge1$ gives the same property; see below |
| $g(N)\to\infty$ | `Tendsto g atTop atTop` with `g : ℕ → ℝ` | Lean's convention $x/0=0$ touches only the finitely many $N$ with $g(N)=0$, which an eventual bound ignores |
| $\limsup 1_A\ast 1_A(n)=\infty$ | `limsup … = ⊤` in `ℕ∞` | `ℕ∞` is a complete lattice, so the limsup always exists; it is `⊤` exactly when $r_A$ is unbounded (`exists_bound_of_limsup_ne_top` and `limsup_ne_top_of_bound` in `lean/Erdos40/Basic.lean`) |
| $A\subseteq\mathbb{N}$ | `A : Set ℕ`, which may contain $0$ | the count ignores $0$, and adding or removing $0$ changes each $r_A(n)$ by at most $2$, so the property does not depend on whether $0\in\mathbb{N}$ |
| for what functions | `erdos_40` equates $\mathcal{G}$ with `answer(sorry)` | the problem asks for $\mathcal{G}$; `erdos_40.variants.weaker` asks whether $\mathcal{G}$ is nonempty |

*The two readings of $\gg$* [written]. If $\sqrt N/g(N)\le c\,A(N)$ for $N\ge N_0$, then $A'=A\cup\{1,\ldots,N_0\}$ satisfies such a bound for every $N\ge1$, with a larger constant. Each element of $A'\setminus A$ lies in at most two ordered pairs summing to $n$, so $r_{A'}(n)\le r_A(n)+2N_0$, and $r_{A'}$ is unbounded exactly when $r_A$ is. The same argument with $\{0\}$ in place of $\{1,\ldots,N_0\}$ gives the row on $0$. On whether $0\in\mathbb{N}$, Thomas Bloom wrote in the problem 28 thread: "I'd rather leave it ambiguous whether $0\in \mathbb{N}$ - sometimes this makes sense, sometimes not."

*The hypotheses are satisfiable, and the property is not trivial* [Lean]. `univ_hypothesis_and_conclusion`: $A=\mathbb{N}$ satisfies the hypothesis for $g=1$, and its representation count is unbounded. `not_erdos40For_sqrt`: $g(N)=\sqrt N$ lacks the property, with witness $A=\{1\}$, for which $A(N)=1=\sqrt N/g(N)$ and $r_A\le1$.

## 2. Theorem A: the existence question is the strong form [Lean]

**Theorem A** (`exists_iff`). Some $g\to\infty$ has the property of problem 40 if and only if $g=1$ has it:

<!-- src: lean/Erdos40/Equivalence.lean -->
```lean
theorem exists_iff :
    (∃ g : ℕ → ℝ, Tendsto g atTop atTop ∧ Erdos40For g) ↔ Erdos40For (fun _ ↦ 1)
```

Its left side is, as a Lean `Expr`, the right side of `erdos_40.variants.weaker`. `answerSet_nonempty_iff` restates it for the set in `erdos_40`: $\mathcal{G}$ is nonempty if and only if the strong form holds. So the weaker variant is not weaker than the strong form; the two are the same statement.

*Proof.* (⇒) If $g\to\infty$ has the property and $\sqrt N\le c\,A(N)$ for large $N$, then $g(N)\ge1$ for large $N$, so $\sqrt N/g(N)\le\sqrt N\le c\,A(N)$ for large $N$, and the property of $g$ makes $r_A$ unbounded.

(⇐) Assume the strong form.

*Step 1* (`uniform_dip`). For every $C\in\mathbb{N}$, $\delta>0$ and $N_0$ there is $M\ge N_0$ such that every $A$ with $r_A(n)\le C$ for all $n$ has some $N\in[N_0,M]$ with $A(N)<\delta\sqrt N$.

Identify a set with its indicator in $\{0,1\}^{\mathbb{N}}$, which is compact in the product topology. For $m\ge0$ let $T_m$ be the set of $A$ with $r_A\le C$ everywhere and $A(N)\ge\delta\sqrt N$ for every $N\in[N_0,N_0+m]$. The condition $r_A(n)\le C$ depends only on $A\cap[0,n]$, and $A(N)\ge\delta\sqrt N$ only on $A\cap[0,N]$, so each is clopen (`isClopen_of_local`) and $T_m$ is closed. Also $T_{m+1}\subseteq T_m$. If Step 1 failed, every $T_m$ would be nonempty, and by compactness some $A$ would lie in all of them. That $A$ has $r_A\le C$ and $A(N)\ge\delta\sqrt N$ for every $N\ge N_0$, so it satisfies the hypothesis of the strong form with $c=1/\delta$ and has bounded $r_A$, which contradicts the strong form.

*Step 2* (`exists_of_strong`). Apply Step 1 with $C=k$, $\delta=(k+1)^{-2}$ and start $s_k$, and call the result $M_k$. Put $s_0=1$ and $s_{k+1}=M_k+1$. Then $s$ is strictly increasing and $s_k\ge k+1$. Let $j(N)$ be the largest $i\le N$ with $s_i\le N$, or $0$ if there is none, and put $g(N)=j(N)+1$. Then $g\to\infty$, because $j(N)\ge b$ once $N\ge s_b$, and $g(N)\le k+1$ for every $N<s_{k+1}$.

Let $A$ satisfy $\sqrt N/g(N)\le c\,A(N)$ for all $N\ge N_1$, and suppose $r_A\le C$ everywhere. Choose an integer $k\ge\max(C,c,N_1)$. Step 1 gives $N\in[s_k,M_k]$ with $A(N)<\sqrt N/(k+1)^2$. Since $N\ge s_k>N_1$ and $N\le M_k<s_{k+1}$,

$$\sqrt N\le c\,g(N)\,A(N)\le(k+1)^2A(N)<\sqrt N,$$

a contradiction. So $r_A$ is unbounded. ∎

Remarks.

* The $g$ of Step 2 is built from the numbers $M_k$, which compactness supplies with no bound. So Theorem A names no member of $\mathcal{G}$ and says nothing about how fast a member must grow.
* The compactness step is standard. I did not find the equivalence stated in the sources of §9, and I did not search further, so I make no claim that it is new.
* `strong_implies_erdos_28`: the strong form implies formal-conjectures' `Erdos28.erdos_28`. The work is formal-conjectures' own `erdos_40.variants.implies_erdos_28`, which is proved there; Theorem A supplies its $g$.
* `Erdos40For.of_isBigO`: if $g$ has the property and $g'=O(g)$ is eventually nonzero, then $g'$ has it. With `not_erdos40For_sqrt` this gives `not_erdos40For_of_sqrt_isBigO`: no $g$ with $\sqrt N=O(g)$ has the property.

## 3. Sets with bounded representation counts, and problems 28 and 158

formal-conjectures defines $B_2[g]$ sets in `FormalConjectures/ErdosProblems/158.lean`, counting unordered representations:

<!-- fc: fc_158.lean -->
```lean
def B2 (g : ℕ) (A : Set ℕ) : Prop :=
  ∀ n, {x : ℕ × ℕ | x.1 + x.2 = n ∧ x.1 ≤ x.2 ∧ x.1 ∈ A ∧ x.2 ∈ A}.encard ≤ g
```

and states problem 158 as

<!-- fc: fc_158.lean -->
```lean
theorem erdos_158 : answer(sorry) ↔ ∀ A : Set ℕ, A.Infinite → B2 2 A →
    liminf (fun N : ℕ => (A ∩ .Iio N).ncard * (N : ℝ) ^ (- 1 / 2 : ℝ)) atTop = 0 := by
```

The results of `lean/Erdos40/Problem158.lean` [Lean]:

* `sumRep_le_of_B2` and `B2_of_sumRep_le`: a $B_2[g]$ set has $r_A\le 2g$, and a set with $r_A\le C$ everywhere is $B_2[C]$. So *$r_A$ is bounded* and *$A$ is $B_2[g]$ for some $g$* say the same thing.
* `count_sq_le`: if $r_A\le C$ then $\lvert A\cap[0,N)\rvert^2\le 2NC$, since every ordered pair from $A\cap[0,N)$ represents some $n<2N$. A set with bounded $r_A$ therefore has at most $\sqrt{2CN}$ elements below $N$, and the question is how far below $\sqrt N$ such a set must dip.
* `strong_iff_forall_B2`: the strong form holds if and only if, for every $g$, every $B_2[g]$ set has
$$\liminf_{N\to\infty}\frac{\lvert A\cap[0,N)\rvert}{\sqrt N}=0,$$
written as in `erdos_158`. `weaker_iff_forall_B2` composes this with Theorem A, so $\mathcal{G}$ is nonempty if and only if the same holds.
* `strong_implies_erdos_158`: the strong form implies the right side of `erdos_158`. The audit of §8 checks that the right side of `strong_iff_forall_B2` at $g=2$ is that right side less the hypothesis `A.Infinite`, which costs nothing: for a finite $A$ the sequence tends to $0$.

Two details of `erdos_158`, checked because the reductions depend on them [written]. It counts $A\cap[0,N)$ where the site counts $A\cap\{1,\ldots,N\}$; the counts differ by at most one, and the proofs absorb the difference. And in `ℝ`, Mathlib's `liminf` is the supremum of the eventual lower bounds, which is $0$ by convention when they are unbounded, so a sequence tending to $+\infty$ has `liminf` equal to $0$. For a set with $r_A\le C$, `count_sq_le` bounds the sequence by $\sqrt{2C}$, so for the sets in question the `liminf` is the genuine one and `= 0` means what it says.

What is known about the right side of `strong_iff_forall_B2` [source]:

* $g=1$, Sidon sets: true. The problem 39 page: "Erdős proved that for every infinite Sidon set $A$ we have" $\liminf_{N} A(N)/N^{1/2}=0$. The problem 158 page: "If we replace $2$ by $1$ then $A$ is a Sidon set, for which Erdős proved this is true." In formal-conjectures, `erdos_158.variants.isSidon` derives this case from `erdos_158.variants.isSidon'`, whose proof is `sorry` and whose docstring says "This is proved in [ESS94]." So the case $g=1$ is not machine-checked there or here.
* A sharper statement for Sidon sets is reported three ways. The problem 1191 page: "Erdős proved (see for example [HaRo66]) that if $A$ is an infinite Sidon set then" $\liminf_x A(x)(\log x)^{1/2}/x^{1/2}\le c$ "for some constant $c>0$.", and the problem, which the page lists as open, asks first whether this liminf is $0$. formal-conjectures' `erdos_158.variants.isSidon'` states the liminf as $<\infty$. Pliego, p. 2: "The above though was proved when g = 1 by Erdős [18, §2 Theorem 8] by showing that in fact any Sidon sequence A satisfies" the same liminf with $=0$, which I read from a render of the page. I read neither Erdős's paper nor [HaRo66], so I leave the conflict unresolved. Problem 40 needs only the unweighted liminf, on which all three agree.
* $g=2$ is problem 158, which is open. For every $g\ge2$ Pliego treats it as a conjecture. He states the Erdős–Turán conjecture as his Conjecture 1.1 and says it "would in turn follow from the stronger conjectural statement (see Erdős and Fuchs [11])" that "for any g ≥ 2 then every B2[g] sequence A ⊂ N has the property that" $\liminf_x A(x)/x^{1/2}=0$, his (1.4), "since asymptotic bases of order 2 cannot satisfy (1.4)". I read (1.4) from a render of p. 2. Cilleruelo, p. 1, writes $a_1<a_2<\cdots$ for the elements of $A$: for Sidon sets the bound $a_k\ll k^2$ is known to fail, and "it is an old open problem to decide whether" $a_k\ll k^2$ "cannot hold for B2[g] sequences with g ≥ 2 either." The bound $a_k\ll k^2$ is equivalent to $A(N)\gg\sqrt N$, so this is the same question.

Consequences [Lean]: one $B_2[g]$ set, for any $g$, with $\liminf A(N)/\sqrt N>0$ would refute the strong form and make $\mathcal{G}$ empty; a negative answer to problem 158 is such a set. A proof that $\mathcal{G}$ is nonempty would answer problems 28 and 158 positively.

## 4. Theorem B: no power of $N$ is in the answer set [written]

**Theorem B.** For every $0<\epsilon<1/2$ there is a set $A\subseteq\mathbb{N}$ with $r_A$ bounded and
$$A(N)>\frac{(N+1)^{1/2-\epsilon}-1}{1-2\epsilon}\quad\text{for all large }N.$$
More precisely, with $k=\lfloor 1/(2\epsilon)\rfloor+1$, $r_A(n)\le 2k-1$ for all large $n$.

*Attribution.* This is a version of the construction the problem 39 page reports: "Erdős and Rényi have constructed, for any $\epsilon>0$, a set $A$ such that" $A(N)\gg_\epsilon N^{1/2-\epsilon}$ "for all large $N$ and $1_A\ast 1_A(n)\ll_\epsilon 1$ for all $n$." O'Bryant's bibliography calls their paper "The seminal paper of Erdős & Rényi [13] introducing the probabilistic method to combinatorial number theory" and cites it as "P. Erdős and A. Rényi, Additive properties of random sequences of positive integers, Acta Arith. 6 (1960), 83–110." Cilleruelo describes their method: each $n$ is put in $A$ independently with a probability $p_n$; they check that $\sum_n P(r_A(n)\ge g+1)<\infty$, with his $r_A$ counting unordered representations, and then finish by "applying the Borel-Cantelli lemma to conclude that with probability 1, the sequences satisfy the B2[g] property after removing a finite number of elements." I did not read their paper. The proof below is my reconstruction of that method with explicit constants, not a transcription.

*Proof.* Fix $\epsilon$, put $\alpha=1/2+\epsilon$, and put each $a\ge1$ in $A$ independently with probability $p_a=a^{-\alpha}$. For $n\ge2$,
$$r_A(n)=2U_n+[\,n/2\in A\,],\qquad U_n=\sum_{1\le a<n/2}[\,a\in A\,][\,n-a\in A\,].$$
The pairs $\{a,n-a\}$ with $a<n/2$ are disjoint, so $U_n$ is a sum of independent Bernoulli variables with mean $\mu_n=\sum_{1\le a<n/2}(a(n-a))^{-\alpha}$.

*Lemma 1.* For $1/2<\alpha<1$ and $n\ge2$, $S(n,\alpha)=\sum_{a=1}^{n-1}(a(n-a))^{-\alpha}\le 2^{2\alpha}n^{1-2\alpha}/(1-\alpha)$.
The map $a\mapsto n-a$ shows $S(n,\alpha)\le 2\sum_{1\le a\le n/2}(a(n-a))^{-\alpha}$. For $a\le n/2$ we have $n-a\ge n/2$, so each term is at most $a^{-\alpha}(n/2)^{-\alpha}$, and $\sum_{1\le a\le n/2}a^{-\alpha}\le\int_0^{n/2}t^{-\alpha}\,dt=(n/2)^{1-\alpha}/(1-\alpha)$. Hence $S(n,\alpha)\le 2(n/2)^{1-2\alpha}/(1-\alpha)=2^{2\alpha}n^{1-2\alpha}/(1-\alpha)$.

*Lemma 2.* $\mu_n\le 2^{2\epsilon}n^{-2\epsilon}/(1/2-\epsilon)$ for $n\ge3$.
The terms of $S(n,\alpha)$ with $a<n/2$ and with $a>n/2$ both sum to $\mu_n$, so $\mu_n\le S(n,\alpha)/2$; substitute $\alpha=1/2+\epsilon$ in Lemma 1.

*Lemma 3.* If $U$ is a sum of independent Bernoulli variables with mean $\mu$, then $P(U\ge k)\le\mu^k/k!$.
$U\ge k$ needs some $k$ of the indicators to equal $1$ together, so by the union bound $P(U\ge k)$ is at most the elementary symmetric polynomial $e_k$ of their probabilities. Expanding $\mu^k=(\sum_i p_i)^k$ produces each product of $k$ distinct $p_i$ exactly $k!$ times, among other nonnegative terms, so $k!\,e_k\le\mu^k$.

*$r_A$ is bounded.* With $k=\lfloor 1/(2\epsilon)\rfloor+1$ we have $2\epsilon k>1$. By Lemmas 2 and 3, $P(U_n\ge k)\le\mu_n^k/k!\le C_\epsilon n^{-2\epsilon k}$, which is summable over $n$. By the Borel–Cantelli lemma, almost surely $U_n\le k-1$ for all large $n$, so $r_A(n)\le 2k-1$ for all large $n$. Since $r_A(n)\le n+1$ for every $n$ (`sumRep_le`), $r_A$ is bounded.

*$A$ is dense.* $X_N=A(N)$ is a sum of independent indicators with mean
$$m_N=\sum_{a=1}^{N}a^{-\alpha}\ge\int_1^{N+1}t^{-\alpha}\,dt=\frac{(N+1)^{1/2-\epsilon}-1}{1/2-\epsilon}.$$
For $t>0$, Markov's inequality applied to $e^{-tX_N}$ and $1-p+pe^{-t}\le\exp(-p(1-e^{-t}))$ give $P(X_N\le m_N/2)\le\exp\big(t\,m_N/2-(1-e^{-t})\,m_N\big)$. At $t=\ln 2$ the exponent is $-m_N(1-\ln2)/2\le -m_N/8$, since $\ln 2<3/4$. As $m_N$ grows like a power of $N$, $\sum_N e^{-m_N/8}<\infty$, so by the Borel–Cantelli lemma, almost surely $X_N>m_N/2$ for all large $N$, which is the stated bound.

Both events have probability $1$, so some $A$ has both properties. ∎

**Corollary.** If $g(N)\ge cN^{\epsilon}$ for all large $N$, for some $c>0$ and $\epsilon>0$, then $g\notin\mathcal{G}$. Equivalently, every $g\in\mathcal{G}$ has $\liminf_N g(N)/N^{\epsilon}=0$ for every $\epsilon>0$.

*Proof.* If $\epsilon\ge1/2$ then $\sqrt N=O(g)$, and `not_erdos40For_of_sqrt_isBigO` applies [Lean]. If $\epsilon<1/2$, take $A$ from Theorem B with the same $\epsilon$. For large $N$, $A(N)\ge c'N^{1/2-\epsilon}$ for some $c'>0$, and $\sqrt N/g(N)\le N^{1/2-\epsilon}/c\le A(N)/(cc')$. So $A$ satisfies the hypothesis of problem 40 for $g$, while $r_A$ is bounded. ∎

*Numerical check of the constants.* `src/theorem_b_check.py` evaluates Lemmas 1 and 2 at every $n$ up to a bound, in floating point with a margin, and recomputes the table of $k$ in exact rational arithmetic. It guards against a slip in a constant and proves nothing asymptotic. The last column of the first table is the limit of the ratio as $n\to\infty$, $2^{2\alpha}/((1-\alpha)\,B(1-\alpha,1-\alpha))$ with $B$ the Beta function, which shows the slack in Lemma 1.

<!-- generated: theorem_b -->
Lemma 1 at every $2\le n\le 20{,}000$:

| $\alpha$ | smallest bound / sum | at $n$ | limit as $n\to\infty$ |
|---|---|---|---|
| 0.505 | 1.2860 | 20,000 | 1.2772 |
| 0.51 | 1.2905 | 20,000 | 1.2811 |
| 0.55 | 1.3284 | 20,000 | 1.3141 |
| 0.6 | 1.3834 | 20,000 | 1.3590 |
| 0.75 | 1.6546 | 20,000 | 1.5255 |
| 0.9 | 2.7399 | 20,000 | 1.7663 |
| 0.99 | 19.8541 | 20,000 | 1.9728 |

Lemma 2 at every $3\le n\le 20{,}000$:

| $\epsilon$ | smallest bound / $\mu_n$ | at $n$ |
|---|---|---|
| 1/200 | 1.2860 | 19,999 |
| 1/100 | 1.2905 | 19,999 |
| 1/20 | 1.3284 | 19,999 |
| 1/10 | 1.3834 | 19,999 |
| 1/5 | 1.5366 | 19,999 |
| 1/4 | 1.6546 | 19,999 |
| 2/5 | 2.7399 | 19,999 |
| 49/100 | 19.8542 | 19,999 |

The choice of $k$, in exact arithmetic:

| $\epsilon$ | $k$ | $2\epsilon k$ | eventual bound $2k-1$ on $r_A$ | density exponent $1/2-\epsilon$ |
|---|---|---|---|---|
| 1/200 | 101 | 101/100 | 201 | 99/200 |
| 1/100 | 51 | 51/50 | 101 | 49/100 |
| 1/20 | 11 | 11/10 | 21 | 9/20 |
| 1/10 | 6 | 6/5 | 11 | 2/5 |
| 1/5 | 3 | 6/5 | 5 | 3/10 |
| 1/4 | 3 | 3/2 | 5 | 1/4 |
| 2/5 | 2 | 8/5 | 3 | 1/10 |
| 49/100 | 2 | 49/25 | 3 | 1/100 |

Every inequality held with the relative float64 margin of 1e-10 that the script requires.
<!-- end generated -->

*The literature.* The argument is not optimised. Cilleruelo: "The exponent 2 + 1/g improves the previous, 2 + 2/g, obtained by Erdős and Renyi in 1960." These exponents are for $a_k$, the $k$-th element of a $B_2[g]$ set, up to factors of lower order, and correspond to densities near $x^{g/(2g+2)}$ and $x^{g/(2g+1)}$. Pliego's Corollary 1.2, which he describes as the first such bound without an $x^{o(1)}$ loss: "For any g ≥ 2 there is a B2[g] sequence A ⊂ N satisfying for large x the estimate" $A(x)\gg x^{g/(2g+1)}$, read from a render of p. 3. For each fixed $g$ these densities are a fixed power below $\sqrt x$, and the power tends to $0$ only as $g\to\infty$. None of the sources I read gives a set with bounded $r_A$ and $A(N)\gg\sqrt N/g(N)$ for a $g$ that grows more slowly than every power of $N$, which is what excluding, say, $g=\log N$ would take.

## 5. The answer set

| statement | label | where |
|---|---|---|
| $\mathcal{G}\ne\emptyset$ if and only if the strong form holds, if and only if for every $g$ every $B_2[g]$ set has $\liminf A(N)/\sqrt N=0$ | [Lean] | `answerSet_nonempty_iff`, `strong_iff_forall_B2`, `weaker_iff_forall_B2` |
| $\mathcal{G}\ne\emptyset$ implies positive answers to problems 28 and 158 | [Lean] | `strong_implies_erdos_28`, `strong_implies_erdos_158` |
| $\mathcal{G}$ is closed downward: $g\in\mathcal{G}$, $g'\to\infty$ and $g'=O(g)$ give $g'\in\mathcal{G}$ | [Lean] | `Erdos40For.of_isBigO` |
| no $g$ with $\sqrt N=O(g)$ is in $\mathcal{G}$ | [Lean] | `not_erdos40For_of_sqrt_isBigO` |
| no $g$ with $g(N)\ge cN^{\epsilon}$ for large $N$ is in $\mathcal{G}$ | [written] | Theorem B, §4 |
| whether $\mathcal{G}$ is empty; which $g$ that grow more slowly than every power of $N$ it contains | open | |

`Erdos40For.of_isBigO` needs $g'$ to be eventually nonzero, which $g'\to\infty$ provides.

## 6. The greedy $B_2[2]$ computation [computation]

In the problem 158 thread, on February 1, 2026, R. Zeraoulia wrote "I posted a computational note on greedy $B_2[2]$ sets here:", linking a Zenodo record (DOI 10.5281/zenodo.18452185), and reported "$a_{2000}=7{,}445{,}662$" for the greedy set, with the ratio $A(N)/\sqrt N$ that "stays around $0.7$ in this range (e.g. $\approx 0.733$ at $N=a_{2000}$)". The comment concludes: "This supports the possibility that $\liminf_{N\to\infty}|A\cap[1,N]|/\sqrt{N}>0$, i.e. a counterexample to Problem 158 may exist." The note itself says the ratio "decreases from 1.8974 at N = 10 to 0.7330 at N = 7,445,662" and "This behavior is consistent with the possibility that the liminf equals 0, but it is far from conclusive". By `strong_iff_forall_B2`, a $B_2[2]$ set with positive liminf would make $\mathcal{G}$ empty, so the question bears on problem 40.

The greedy set is $a_1=1$, and $a_{k+1}$ is the least $m>a_k$ such that every $n$ keeps at most two representations $n=a+b$ with $a\le b$. `src/greedy_b2g.c` computes it, and `src/greedy_check.py` checks the output without sharing code with the C program: for $g=1$ the program must reproduce OEIS A005282, whose entry reads "Mian-Chowla sequence (a B_2 sequence)"; a separate Python implementation must reproduce a prefix for $g=2$; the whole output must be a $B_2[2]$ set, counted from all pairwise sums; the thread's two numbers and the note's Table 1 must agree to the precision printed; and the kernel-enforced CPU limit must stop a run when lowered to one second. `src/final_check.py` rechecks the committed list in a third implementation (§8).

<!-- generated: greedy_table -->
The first 10,000 terms, with candidates up to 400,000,000; the last is 376,437,170.

| $k$ | $a_k$ | $k/\sqrt{a_k}$ | $a_{k/2}$ | smallest $A(N)/\sqrt N$ on $[a_{k/2},a_k]$ | local exponent |
|---|---|---|---|---|---|
| 100 | 5,743 | 1.3196 | 1,096 | 1.2893 |  |
| 200 | 29,272 | 1.1690 | 5,743 | 1.1524 | 2.3496 |
| 500 | 263,671 | 0.9737 | 50,287 | 0.9706 | 2.3989 |
| 1,000 | 1,397,299 | 0.8460 | 263,671 | 0.8444 | 2.4058 |
| 2,000 | 7,445,662 | 0.7330 | 1,397,299 | 0.7326 | 2.4138 |
| 5,000 | 69,376,857 | 0.6003 | 12,910,949 | 0.6000 | 2.4358 |
| 10,000 | 376,437,170 | 0.5154 | 69,376,857 | 0.5150 | 2.4399 |
<!-- end generated -->

The note's Table 1 against the computed set (the counts are $A(N)$, and the note prints its ratios to four places):

<!-- generated: note_table -->
| $N$ | $A(N)$, note | $A(N)$, computed | ratio, note | ratio, computed |
|---|---|---|---|---|
| 10 | 6 | 6 | 1.8974 | 1.8974 |
| 100 | 17 | 17 | 1.7000 | 1.7000 |
| 1,000 | 48 | 48 | 1.5179 | 1.5179 |
| 10,000 | 127 | 127 | 1.2700 | 1.2700 |
| 100,000 | 332 | 332 | 1.0499 | 1.0499 |
| 1,000,000 | 870 | 870 | 0.8700 | 0.8700 |
| 7,445,662 | 2,000 | 2,000 | 0.7330 | 0.7330 |
<!-- end generated -->

*Reading.* The posted numbers are right. Past the range of the comment the ratio keeps falling, to 0.515 at the 10,000th term, and within each window the smallest ratio is close to the ratio at the window's right end. The local exponent $\log(a_k/a_{k'})/\log(k/k')$, between consecutive checkpoints $k'<k$, rises slowly and is above $2$ at every checkpoint. If it stayed above some $2+\eta$ the ratio would tend to $0$ for this set, but no finite range shows that, and the behaviour of one greedy set says nothing about all $B_2[2]$ sets. So the computation does not bear on problem 158, and through it on $\mathcal{G}$, in either direction.

An exploratory run for more terms stopped at its candidate cap before finishing. It is not committed, and nothing here uses it.

## 7. Status and the blocking question

Problem 40 is open. Its existence half, whether $\mathcal{G}$ is nonempty, is by Theorem A exactly the strong form, and by §3 exactly the statement that for every $g$ every $B_2[g]$ set has $\liminf A(N)/\sqrt N=0$. For $g=1$ that is Erdős's theorem, for $g=2$ it is problem 158, and for $g\ge2$ in general it is the conjecture Pliego attributes to Erdős and Fuchs. The blocking question is therefore:

> Does every $B_2[2]$ set $A$ have $\liminf_{N\to\infty}A(N)/\sqrt N=0$, and likewise for every larger $g$?

A negative answer for any one $g$ empties $\mathcal{G}$. A positive answer for all $g$ makes $\mathcal{G}$ nonempty without naming a member. For the *for what functions* half, the only exclusions established here are $g$ with $g(N)\ge cN^{\epsilon}$; whether $\log N$ or any other slowly growing $g$ lies in $\mathcal{G}$ is untouched, and deciding it either way needs a construction or a proof that none of the sources read contains.

## 8. Verification

*The Lean proofs.* `lean/verify.sh` runs five steps and stops at the first failure.

0. The formal-conjectures checkout is the commit pinned in `lean/lakefile.toml`, with no local changes. The statement checks compare against formal-conjectures as built here, so this comes first.
1. `lake build`, which includes `lean/Erdos40/Audit.lean`. There, `#assert_axioms` fails the build if a listed theorem depends, transitively, on any axiom outside `propext`, `Classical.choice` and `Quot.sound`, `sorryAx` included; and a `run_cmd` block compares six statements, as `Expr`s, with the formal-conjectures declarations they are meant to match.
2. Re-elaborates `Audit.lean` and requires one clean axiom report per listed theorem and the statement check's pass message.
3. Runs `lean/Probe.lean`, which imports only `Lean` and loads the compiled modules at runtime, so nothing they define can change how it runs. It checks the axioms of every constant the three modules add, auxiliary ones included, requires the main theorems to exist, and repeats the statement comparisons.
4. Replays `Erdos40.Basic`, `Erdos40.Equivalence` and `Erdos40.Problem158` through the kernel with `leanchecker`, one module per call, in a separate process.

The last committed run, from `results/lean_verify.log`:

<!-- generated: lean_summary -->
Run on September 28, 2026 at 11:25:03 (UTC−04:00) with Lean 4.33.1 on Darwin arm64: `verify.sh` exited 0. All 26 listed theorems depend on no axiom beyond `propext`, `Classical.choice` and `Quot.sound`; the statements of `exists_iff`, `answerSet_nonempty_iff`, `strong_implies_erdos_28`, `strong_implies_erdos_158`, `strong_iff_forall_B2`, `weaker_iff_forall_B2` match formal-conjectures at the pin; the probe checked 74 constants and found 0 problems; and the kernel replayed `Erdos40.Basic`, `Erdos40.Equivalence` and `Erdos40.Problem158`. Wall time 422.98 s; maximum resident set size 5,275,664,384 bytes, as `/usr/bin/time -l` reports it.
<!-- end generated -->

*Planted defects.* `lean/plants.py` tests `verify.sh`. For each case it copies the project to a scratch directory, links the built packages (a private copy of formal-conjectures for the two cases that alter it), applies one defect, and runs `verify.sh` there. A case passes only if `verify.sh` exits 1 at the expected step with the expected message; the control, an unmodified copy, must exit 0 and print `PASS:`. The four statement cases also run the probe alone and require it to catch the change. `results/lean_plants.json` records each case, with the hashes of the sources it tested.

<!-- generated: plants_table -->
Run on September 28, 2026 at 11:32:11 (UTC−04:00) with Lean 4.33.1 on macOS-27.0-arm64-arm-64bit-Mach-O: 12 of 12 cases behaved as expected.

| case | planted | expected | observed | probe alone |
|---|---|---|---|---|
| `control` | nothing: an unmodified copy | exit 0 | exit 0 |  |
| `pin` | lakefile.toml pins a different formal-conjectures commit | exit 1 at step 0 | exit 1 at step 0 |  |
| `fc_edited` | a comment appended to formal-conjectures' 40.lean | exit 1 at step 0 | exit 1 at step 0 |  |
| `sorry_listed` | `sorry` for the proof of `sumRep_le`, which the audit lists | exit 1 at step 1 | exit 1 at step 1 |  |
| `sorry_unlisted` | a new theorem proved by `sorry`, not listed in the audit | exit 1 at step 3 | exit 1 at step 3 |  |
| `axiom_unlisted` | a new axiom and a theorem that uses it, not listed in the audit | exit 1 at step 3 | exit 1 at step 3 |  |
| `kernel_off` | a proof of False added with the kernel check switched off, the option name built at runtime, and a theorem using it | exit 1 at step 4 | exit 1 at step 4 |  |
| `statement_lhs` | `exists_iff` restated with its conjunction swapped, still fully proved | exit 1 at step 1 | exit 1 at step 1 | exit 1 |
| `statement_rhs` | `exists_iff` restated with `1 * 1` for `1`, still fully proved | exit 1 at step 1 | exit 1 at step 1 | exit 1 |
| `statement_158` | `strong_implies_erdos_158` restated for `B2 3` sets, still fully proved | exit 1 at step 1 | exit 1 at step 1 | exit 1 |
| `statement_forall_B2` | `strong_iff_forall_B2` restated with `B2 (g + 0)`, definitionally the same, still fully proved | exit 1 at step 1 | exit 1 at step 1 | exit 1 |
| `statement_hijack` | the swapped statement plus a macro that replaces the audit's `run_cmd` with one printing the pass message | exit 1 at step 3 | exit 1 at step 3 |  |
<!-- end generated -->

`sorry_unlisted`, `axiom_unlisted` and `statement_hijack` pass steps 1 and 2, which only see what the audit lists or can be subverted from inside the build, and are caught by the probe. `kernel_off` passes steps 1 to 3 and is caught only by the kernel replay.

*Environment.*

<!-- generated: lean_env -->
| component | revision | requested as |
|---|---|---|
| Lean | `leanprover/lean4:v4.33.1` | `lean/lean-toolchain` |
| formal_conjectures | `e5f428182a3ee32dde401eceb4a94ba0382e8434` | `e5f428182a3ee32dde401eceb4a94ba0382e8434` |
| mathlib | `0df444a360eaa60ab8c11dca51a86af692955474` | `v4.33.1` |
| plausible | `b7eb3304aeae834b12dda98993a37f6a41f6f0bb` | `main` |
| LeanSearchClient | `5f4d51b81cbd3f6b32b156bfad9056621a040404` | `main` |
| importGraph | `16f02aa7642864af59f1ff0e384a015994db9118` | `main` |
| proofwidgets | `4be2e3d5087eeb272cf5a8853b8f9dd025ef5957` | `main` |
| aesop | `3448c0bcc5ce01b2d1546e483ec3620e32df3d0e` | `master` |
| Qq | `92c15be17b7caf78c2ad767ec40f89052d908d81` | `master` |
| batteries | `4488d40d070b9700d4d5a6aa342f0d40c31b2a2d` | `main` |
| Cli | `6130a47896ce867c6a4a55373441e59e565bad0f` | `v4.33.0` |
<!-- end generated -->

`lean/lake-manifest.json` is committed, so these revisions are fixed. A fresh clone has not been tested here: `lake build` fetches the packages the manifest pins, and Mathlib's `lake exe cache get` downloads its compiled files instead of rebuilding them.

*What is trusted.* The kernel replay covers this project's three modules, not their imports: `leanchecker` ran without `--fresh`, so Mathlib and formal-conjectures are trusted as compiled at the pin. The informal-to-formal correspondence of §1.2 is argued, not machine-checked, and the Lean statements are checked against formal-conjectures' statements, not against the site's. Theorem B, the two readings of $\gg$, and the remarks on `erdos_158` in §3 are written proofs.

*The gate.* `src/final_check.py` is the gate for everything else. It recomputes the greedy checkpoints and the note's table from the committed list, rechecks that the list is a $B_2[2]$ set and that its first terms follow the greedy rule, recomputes the Theorem B constants in exact and multi-precision arithmetic, and checks the Lean records against the current sources. In the documentation it checks the statement against the site, each Lean excerpt and each quotation against its source, and every numeral in the prose. A numeral passes if it agrees with a freshly computed value, lies in a passage that a fetched source or a committed record backs, or names a problem, section, lemma, step, exit status, release or identifier that exists. Small integers inside formulas are mathematics and are left to the proofs, ordered-list markers are structure, and the model, the CLI version and the session's date named in §10 cannot be checked, so the gate prints them. It also regenerates the tables above. Before it reports a pass it plants a defect for each check and requires each to be caught.

## 9. Sources

Read, and checked mechanically where quoted (`refs/README.md` says which claim each file supports):

* erdosproblems.com: the pages of problems 28, 39, 40, 158 and 1191, the LaTeX view of problem 40, and the discussion threads of problems 28, 40 and 158, fetched September 23, 2026. The site is maintained by Thomas Bloom.
* formal-conjectures at `e5f428182a3ee32dde401eceb4a94ba0382e8434`: `ErdosProblems/28.lean`, `40.lean`, `158.lean` and `Combinatorics/Additive/Convolution.lean`.
* J. Pliego, *On the Erdős-Turán Conjecture and the growth of $B_{2}[g]$ sequences*, [arXiv:2405.04154](https://arxiv.org/abs/2405.04154) (v1, May 7, 2024; the only version, with no journal reference).
* J. Cilleruelo, *Probabilistic constructions of B2[g] sequences*, preprint dated June 3, 2008, from the author's page at the Universidad Autónoma de Madrid.
* K. O'Bryant, *A Complete Annotated Bibliography of Work Related to Sidon Sequences*, [arXiv:math/0407117](https://arxiv.org/abs/math/0407117) (v1, July 8, 2004).
* R. Zeraoulia, *Computational Evidence for Erdős Problem #158 via the Greedy B2[2] Construction*, Zenodo, February 1, 2026, [doi:10.5281/zenodo.18452185](https://doi.org/10.5281/zenodo.18452185), the file `erdos_problem_158_greedy_computation_note_v7.pdf`.
* OEIS [A005282](https://oeis.org/A005282), the entry and its b-file.

Not read. Everything attributed to these comes from the sources above.

* [Er95] P. Erdős, *Some of my favourite problems in number theory, combinatorics, and geometry*, Resenhas (1995), 165-186, MR 1370501; and [Er97c] P. Erdős, *Some of my favorite problems and results*, The mathematics of Paul Erdős, I (1997), 47-67, MR 1425174. These are the problem's sources on the site. The Erdős archive at the Rényi Institute did not serve them, and the only copy of [Er97c] I found was an upload that did not appear to be authorized, which I did not use.
* P. Erdős and A. Rényi (1960), cited above from O'Bryant. The PDF that EuDML links sits behind a proof-of-work bot challenge, which I did not try to get past.
* Erdős and Fuchs, Erdős's paper cited by Pliego as [18], [HaRo66] and [ESS94]: known here only through the citations quoted in §3.
* Problem C9 of Guy's collection, [Gu04] on the site, which the problem 39 page cites.

## 10. Provenance

The work was done by Claude Opus 5.5 in GitHub Copilot CLI 1.0.88 on September 23, 2026 (US Eastern time), in one session at the direction of the repository owner, who supplied the problem and the working instructions but no mathematical hints. The session's own event log records the tools it ran: a shell, file views, searches and edits, a session database for the task contract, a scheduler for check-ins against that contract, and the web search and web fetch tools. The downloads in `refs/fetch.sh` ran in the shell. The web search tool answers with a model-written summary, so it served only to find sources; every source cited here was then fetched and read. No subagent was used. The session itself ran every checker whose results are committed: `lean/verify.sh`, `lean/plants.py`, `src/greedy_check.py`, `src/theorem_b_check.py` and `src/final_check.py`.

Failures on the way, all fixed before the committed runs: a plant in `lean/plants.py` that restated a theorem while another file still used its original name, so the plant failed for the wrong reason; an interrupted plant run that left its Lean process group running; quotation checks that missed because PDF text layers write ő with a combining accent, fixed by normalising both sides to NFC; a gate check that backed the session's date with the site's own access stamp, which is not in US Eastern time and changes with every fetch, so the end-to-end run failed on it; an excerpt check that passed, rather than skipping, on a copy of the repository without the formal-conjectures files; and the capped exploratory run of §6. The cost in model credits was not available to the session when this was written.
