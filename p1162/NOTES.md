# Erdős problem #1162 — working notes

> **Give an asymptotic formula for the number of subgroups of $S_n$. Is there a statistical theorem on their order?** — a problem of Erdős and Turán, source [Va99, 5.73]

**Status: open. Nothing here solves it.** The tracker records the problem as `OPEN` with the note *"This is open, and cannot be resolved with a finite computation."* That is correct, and no computation below is offered as a proof of an asymptotic statement.

What this directory does contain:

1. a correction to an ambiguity in the problem statement as displayed (§1: the logarithm is base 2);
2. exact counts, far past the published ranges, of the family of subgroups that carries the whole first-order answer (§2–§3);
3. two limiting constants identified in closed form and confirmed to 9–12 decimal places (§3, §4; the first version pointed to §4, §5, but §5 is a table of lower bounds);
4. a second-order term derived from the orbit structure and confirmed by out-of-sample extrapolation (§6);
5. an exactly located lossy step in the one published theorem that addresses the *second* half of Erdős's question, and a conjecture for what replaces it (§7);
6. the evidence that does **not** support the conjecture, stated plainly (§8);
7. a second, larger family, Frattini-closed 2-subgroups with every orbit of size at most 8, counted exactly to $n = 512$, and two short proofs: $\log_2 f(n) \ge n^2/16 + \tfrac78 n\log_2 n - O(n)$, and $e(n)/f(n) \to 0$. The second means the elementary abelian family of items 2–5 is not typical. Conjecture 7.2 is withdrawn in consequence and replaced by Conjecture 9.3, which puts a typical subgroup at $2^{3n/8}$ rather than $2^{n/4}$ (§9);
8. a census of every 2-subgroup of $S_n$ for $n \le 15$, Frattini-closed or not, whose shares run against premise (b) of Conjecture 9.3 and leave (c) untested (§9.8, added November 1, 2026).

`src/final_check.py` re-derives the numbers quoted below and exits non-zero if one is wrong. It compares every decimal number in this file with a fresh computation, to the precision printed, except the timings in §10 and the wrong values that correction notes quote from first versions; `python src/final_check.py --coverage` lists those exceptions. Integers and fractions are compared where a check reads them, which is not everywhere. Not checked: the literature and citation searches, the outputs of the deleted program in §2.1, and the reasoning of §7 and of the proofs in §9, whose arithmetic is checked. The first version of this paragraph said *"Every number quoted below is re-derived"*. That was untrue: the §5 and §8 tables, the §6 ε values and the §3 rows for $n = 64$ and for the ratios were not checked, and the §5 entry for $n = 128$ was wrong (0.5276; the value is 0.5275). Those corrections, and the others made on September 22 and 23, 2026, are marked where they occur.

---

## 1. The statement: which logarithm?

The problem page states Pyber's result as `log f(n) ≍ n²` and Roney-Dougal–Tracey's as `log f(n) = (1/16 + o(1))n²` without fixing a base, and the constant 1/16 is meaningless without one.

It is base 2. From page 1 of [RoTr25], read both from `pdftotext -raw` and from a 150 dpi render of the page:

> "|Sub(S_n)| of subgroups of S_n is at most 2^{ξn²+o(n²)} where ξ = ⅙ log 24 (our logarithms are to the base 2). Pyber conjectured, however, that
> |Sub(S_n)| = 2^{n²/16+o(n²)}."

Independently, the OEIS A005432 comment records Pyber's brackets as `c^{n²(1+o(1))} ≤ a(n) ≤ d^{n²(1+o(1))}` with `c = 2^{1/16}`, `d = 24^{1/6}`, which is the same statement with the base made explicit.

*(The default and `-layout` modes of `pdftotext` both sort text by position and splice superscripts from the displayed math into this sentence — it comes out as "are to the 2 2 base 2)" — so a verbatim search against them fails even though the sentence is on the page. `src/final_check.py` uses `-raw`.)*

**All logarithms in these notes are base 2.**

### What is and is not answered

| Half of #1162 | Status |
| --- | --- |
| Asymptotic formula for $f(n) = \|\mathrm{Sub}(S_n)\|$ | Largely answered. [RoTr25] Theorem 1: $2^{n^2/16+\alpha n\log n} \le f(n) \le 2^{n^2/16+\beta n^{3/2}}$ for absolute constants $\alpha>0,\beta$. The second-order term is not pinned down: $n\log n$ vs $n^{3/2}$. Theorem 9.1 here puts the $n\log n$ coefficient at no less than $7/8$ asymptotically, and Conjecture 9.3 is that $7/8$ is exact. |
| Statistical theorem on their order | Barely touched. [RoTr25] Theorem 6 is the only result and it is one-sided. §7 locates the lossy step in its proof and conjectures larger constants. The first version of this row said §7 *shows* the constants are not optimal; §7 does not show that. §9.6–9.8 give a heuristic answer: $\log_2\|H\| \approx 3n/8$. |

## 2. What was computed exactly

### 2.1 The full subgroup lattice, $n \le 8$

`src/subgroups.c` enumerates every subgroup of $S_n$ by expanding one representative per conjugacy class and closing under conjugation.

| $n$ | subgroups $f(n)$ | classes | agrees with |
| --- | --- | --- | --- |
| 1–7 | 1, 2, 6, 30, 156, 1455, 11300 | 1, 2, 4, 11, 19, 56, 96 | A005432 / A000638 |
| 8 | 151221 | 296 | A005432(8) / A000638(8) |

$n=8$ took 883 s. This is ground truth for the *actual* object of #1162 — the full multiset of subgroup orders — and it is the anchor that the much larger elementary abelian computation is checked against.

**A wrong version of this program is preserved in the history of this directory and is worth recording.** The first attempt pruned with a "covered" set, marking `⟨H,g⟩` covered for each new generator `g`. It produced 1, 2, 6, 30, **144, 1066, 6415** — wrong from $n=5$ on, 43 % low at $n=7$ — because if `g ∈ ⟨H,g'⟩` for an earlier `g'`, then `⟨H,g⟩` may still be a different and *smaller* subgroup, so genuine minimal overgroups were dropped. The rewrite prunes by right cosets instead (`⟨H,g⟩ = ⟨H,hg⟩` for `h ∈ H`, so the whole coset `Hg` may be skipped), which is sound. The error was invisible by inspection and was caught only by comparing against OEIS.

### 2.2 Elementary abelian 2-subgroups, $n \le 512$

Write $a(n,k) = \#\{H \le S_n : H \cong (\mathbb{Z}/2)^k\}$ and $e(n) = \sum_k a(n,k)$. `src/elemab.py` computes these exactly:

* $\#\mathrm{Hom}(G,S_n) = n!\,[x^n]\exp\!\big(\sum_{H\le G} x^{[G:H]}/[G:H]\big)$; for $G=\mathbb{F}_2^k$ the index-$2^j$ subgroups are the codimension-$j$ subspaces, so $a_{2^j} = \binom{k}{j}_2 (2^j-1)!$ and only $j \le \log_2 n$ contribute — each coefficient costs $O(\log n)$ big-integer multiplications;
* $\#\mathrm{Inj}$ from $\#\mathrm{Hom}$ by $q$-binomial inversion at $q=2$, $i_k = \sum_m (-1)^{k-m} 2^{\binom{k-m}{2}} \binom{k}{m}_2 h_m$;
* $a(n,k) = i_k / |\mathrm{GL}_k(2)|$.

For example `rank_counts(8) = [1, 763, 6685, 4440, 350]`, and $e(512)$ is an exact integer of 18 997 bits.

**Why this family.** It is not a proxy chosen for convenience: it already carries the entire first-order answer. $\log e(n) = n^2/16 + \Theta(n\log n)$ (§6), and $\log f(n) = n^2/16 + o(n^2)$, so elementary abelian 2-subgroups account for all of the leading term. §8 records what that does *not* imply, and §9.5 proves the gap is real: $e(n)/f(n) \to 0$, so the family carries the leading term without being typical.

### 2.3 Verification chain for $a(n,k)$

Four independent paths, all agreeing:

| check | range | status |
| --- | --- | --- |
| `src/verify_elemab.py`, brute-force permutation enumeration sharing no code with `elemab.py` | $n \le 10$ | all ranks match ($n=10$: `[1, 9495, 330015, 620325, 158220, 5670]`, 397 s) |
| `src/subgroups.c` lattice, a third implementation | $n \le 8$ | all ranks match |
| $a(n,1) = A000085(n) - 1$ (involutions, independent recurrence) | $n \le 30$ | matches |
| structural assertion $a(n,\lfloor n/2\rfloor+1) = 0$ | $n \le 512$ | holds — a sharp test of a massively cancelling alternating sum |

## 3. The rank distribution: an exact statistical theorem for this family

For each $n \le 512$ the distribution $k \mapsto a(n,k)/e(n)$ is known exactly, so the following are *computations*, not estimates:

| $n$ | mode | mean$/n$ | sd | window holding $\ge 1-10^{-6}$ |
| --- | --- | --- | --- | --- |
| 64 | 16 | 0.249997 | 0.849353 | $\pm 4$ |
| 128 | 32 | 0.250000 | 0.849306 | $\pm 4$ |
| 256 | 64 | 0.250000 | 0.849306 | $\pm 4$ |
| 512 | 128 | 0.250000 | 0.849306 | $\pm 4$ |

The mode is exactly $n/4$ for every $n \equiv 0 \pmod 4$ tested from 32 up. **The width does not grow with $n$**: the standard deviation converges to an absolute constant and the $1-10^{-6}$ window stays at $\pm 4$ while $n$ grows by a factor of 8.

The limit law is identified in closed form. For $n \equiv 0 \pmod 4$,

$$\frac{a(n,\ n/4+j)}{a(n,\ n/4)} \longrightarrow 2^{-j^2},$$

and this holds to better than $10^{-6}$ already at $n=128$:

| $j$ | $-3$ | $-2$ | $-1$ | $0$ | $+1$ | $+2$ | $+3$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $n=512$ | 0.001953 | 0.062500 | 0.500000 | 1 | 0.500000 | 0.062500 | 0.001953 |
| $2^{-j^2}$ | 0.001953 | 0.062500 | 0.500000 | 1 | 0.500000 | 0.062500 | 0.001953 |

The standard deviation of that law is

$$\sqrt{\frac{\sum_{j\in\mathbb{Z}} j^2 2^{-j^2}}{\sum_{j\in\mathbb{Z}} 2^{-j^2}}} = 0.8493059610\ldots$$

against a measured 0.849305961022908 — and the agreement is far better than the stored artifact suggests. `data/elemab_stats.json` keeps `sd` as a double, and the double is much less precise than its format: `ratio` in `src/asymptotics.py` divides two integers through their floating-point base-2 logarithms, which are in the thousands for these $n$, and the variance then cancels about four more digits. The cautions below give the resulting errors. Recomputing the variance from the exact integer rank counts at 60 decimal places gives

| $n$ | 64 | 128 | 256 | 512 |
| --- | --- | --- | --- | --- |
| $\|\mathrm{sd}(n) - \text{closed form}\|$ | 4.7e−5 | 8.0e−9 | 1.0e−17 | 2.7e−36 |
| decimal places | 4 | 7 | 16 | **35** |

The first version printed 1.1e−17 for $n = 256$, from rounding the stored 1.05e−17 a second time; the value is $1.047\times10^{-17}$. It also said that the double *"truncates at about $10^{-16}$ relative"*, which is the precision of the format, not of the stored values.

At $n=512$ the exact standard deviation agrees with $0.849305961022908280743697585694$ in all 30 digits shown. The error is roughly squaring at each doubling of $n$ — the $\log_{10}$ errors are $-4.3, -8.1, -17.0, -35.6$, ratios $1.87, 2.10, 2.09$ — i.e. it decays like $2^{-cn}$ with $c \approx 0.23$, not polynomially. §9.3 argues, but does not prove, that the error is of order $n^2 2^{-n/4}$. Under that rate the apparent $c$ between $n$ and $2n$ would be $1/4 - 2/n$, which is 0.219, 0.234 and 0.242 for $n = 64, 128, 256$; the measured values are 0.196, 0.230 and 0.241. That is consistent with $c \to 1/4$ but does not establish it. The first version printed the last of the three ratios as 2.10; it is 2.0947.

Two cautions. The deviations read off the stored doubles instead (2.8e−10 at $n=256$, 1.2e−9 at $n=512$, and non-monotone) are **float noise in the artifact, not convergence error**; I reported them as convergence before checking, and they are wrong. And this is convergence of an exact finite quantity to a closed form, checked at four values of $n$ — it is extremely strong numerical evidence, but it is not a proof. The limit law itself is now proved, as Proposition 9.5(d) in §9.10; the rate of convergence is not.

**Where it comes from.** $a(n, n/4+j)$ behaves like $2^{nk/2-k^2}$ near $k=n/4$ (§6), and putting $k = n/4+j$ gives $n^2/16 - j^2$: the quadratic in $k$ has its vertex at $n/4$, and the $2^{-j^2}$ profile is exactly the same mechanism that makes $\binom{r}{r/2+j}_2 \approx 2^{r^2/4-j^2}$. §9.3 makes this literal: numerically the rank law is that of the dimension of a uniformly random subspace of $\mathbb{F}_2^{\lfloor n/2\rfloor}$, to a total variation distance of order $n^2 2^{-n/4}$.

So, **for this family**, the answer to the second half of #1162 is sharp: $\log_2 |H| = n/4 + O(1)$, with an explicit limiting distribution. The family is not typical (§9.5), and in the larger family of §9.6 the order concentrates near $3n/8$ instead, with a spread of order $n^{1/4}$.

## 4. A single maximal elementary abelian is not enough

The classical Pyber lower bound takes one maximal elementary abelian $E \cong \mathbb{F}_2^{\lfloor n/2\rfloor}$ and counts its subspaces, giving the Galois number $G_{\lfloor n/2\rfloor}$. Exactly:

$$\log_2 G_r - \frac{r^2}{4} \longrightarrow \log_2\!\frac{\vartheta}{\eta} = 2.8820499654\ldots \ (r \text{ even}),
\qquad \vartheta=\sum_{j\in\mathbb{Z}}2^{-j^2},\quad \eta=\prod_{i\ge1}(1-2^{-i}).$$

For odd $r$ the limit is $\log_2(\vartheta'/\eta) = 2.8820461863\ldots$, with $\vartheta' = \sum_{j\in\mathbb{Z}}2^{-(j+1/2)^2} = 2.1289312505\ldots$ in place of $\vartheta = 2.1289368272\ldots$; here $\eta = 0.2887880951\ldots$ and $\vartheta/\eta = 7.372$. Measured: $2.882049965$ at $r = 128, 256, 512$ (even) and $2.882046186$ at $r = 127, 255, 511$ (odd), stable to 9 decimals; the deviation from the closed form falls like $2^{-r/2}$ ($5.2\times10^{-5}$ at $r=32$, $8.0\times10^{-10}$ at $r=64$). The first version gave only the even limit and quoted the odd measurement beside it as if it converged to the same constant; the two limits differ by $3.8\times10^{-6}$.

$\eta^{-1} = 3.46275\ldots$ is exactly the constant $\zeta_2$ of [RoTr25] Lemma 2.4 ("$P$ has at most $\zeta_p p^{k(\ell-k)}$ subgroups of order $p^k$ … $\zeta_p < 4$").

**Consequence.** One maximal $E$ contributes $n^2/16 + O(1)$ — a *bounded* excess. So the entire $\alpha n\log n$ of [RoTr25] Theorem 1's lower bound comes from the multiplicity of elementary abelian subgroups, not from any single one.

## 5. An explicit lower bound for the constant $\alpha$

[RoTr25] Theorem 1 gives $f(n) \ge 2^{n^2/16+\alpha n\log n}$ with $\alpha>0$ unspecified. Since $e(n) \le f(n)$ and $e(n)$ is exact, every computed $n$ certifies a value, and the two families of §9 certify larger ones: $B(n)$ is the Frattini-closed family of §9.1 and $F_8(n)$ its members with every orbit of size 8 carrying a conjugate of 8T22 (§9.4).

| $n$ | 16 | 32 | 64 | 128 | 256 | 512 |
| --- | --- | --- | --- | --- | --- | --- |
| $(\log_2 e(n)-n^2/16)/(n\log_2 n)$ | 0.4487 | 0.4790 | 0.5039 | 0.5275 | 0.5487 | 0.5671 |
| $(\log_2 F_8(n)-n^2/16)/(n\log_2 n)$ | 0.4075 | 0.5076 | 0.5596 | 0.5998 | 0.6321 | 0.6581 |
| $(\log_2 B(n)-n^2/16)/(n\log_2 n)$ | 0.5463 | 0.5711 | 0.5972 | 0.6232 | 0.6467 | 0.6673 |

The $n = 128$ entry of the first row read 0.5276 in the first version; the value is 0.527548. Each entry is a statement about one $n$, not a proof that $\alpha \ge 0.667$ for all $n$, and no such proof is possible: Theorem 1 asks for an $\alpha$ valid for all $n > 1$, and $n = 2$ alone forces $\alpha \le 3/8$ (§9.4). What is proved is asymptotic. The second row tends to $7/8$, so the $\liminf$ of the ratio for $f(n)$ is at least $7/8$ (Theorem 9.1). All three rows are still rising at $n = 512$.

## 6. The second-order term is $\tfrac34 n\log_2 n$

**Derivation (orbit structure).** Let $\mathbb{F}_2^k$ act on $[n]$ with every orbit of size $m=2^j$. Exactly $N_j(n,k) = \frac{n!}{t!}\big(\binom{k}{j}_2/m\big)^t$ with $t=n/m$, and dividing by $|\mathrm{GL}_k(2)|$,

$$\log_2 \frac{N_j(n,k)}{|\mathrm{GL}_k(2)|} = (1-2^{-j})\,n\log_2 n \;+\; \frac{jkn}{2^j} - k^2 \;+\; O_j(n).$$

Maximising the $n^2$ part gives $k = jn/2^{j+1}$ with value $j^2n^2/2^{2j+2}$:

| $j$ | orbit size | optimal $k/n$ | $n^2$ coefficient | $n\log_2 n$ coefficient |
| --- | --- | --- | --- | --- |
| 1 | 2 | 1/4 | **1/16** | 1/2 |
| 2 | 4 | 1/4 | **1/16** | **3/4** |
| 3 | 8 | 3/16 | 9/256 | 7/8 |
| 4 | 16 | 1/8 | 1/64 | 15/16 |

Sizes 2 and 4 **tie** at the leading constant $1/16$ — the constant of Pyber and of Theorem 1 — and size 4 wins the tie-break at order $n\log n$. Mixtures do not help: with $c_j$ the fraction of points in orbits of size $2^j$, the $n^2$ coefficient is $\kappa\sum_j c_j j2^{-j} - \kappa^2$ and $\sum_j c_j j 2^{-j} \le 1/2$ with equality iff $c_j$ is supported on $\{1,2\}$; on that face the $n\log n$ coefficient is $1/2 + c_2/4$, maximised at $c_2=1$. Hence

$$\log_2 e(n) = \frac{n^2}{16} + \frac34 n\log_2 n + O(n).$$

§9.5 turns this display into a theorem: Proposition 9.2 proves the upper bound, and the construction in the proof of Theorem 9.1, run with 4-blocks, proves the lower bound. The derivation above remains the explanation of the constant.

**Check 1 — the closed form.** `pure_type_log2` was validated against 15 brute-forced hom counts (every $(n,k,j)$ with $j \ge 1$, $2^j \mid n$, and either $n \in \{4,6\}$, $1 \le k \le 3$ or $n = 8$, $1 \le k \le 2$), exact match. The first version said "all $(n,k,j)$ with $n\le 8$, $k\le3$", which is more than was brute-forced. Its optimal rank equals the predicted $j/2^{j+1}$ *exactly* in every case tested ($j=1,2,3,4$; $n=64,128,256$), and the measured $\mathrm{excess}/(n\log_2 n)$ approaches $1-2^{-j}$ with an $O(n)$ correction whose size is confirmed constant after multiplying by $\log_2 n$: the values of $(1 - 2^{-j} - \mathrm{excess}/(n\log_2 n))\log_2 n$ at $n = 64, 128, 256$ are 0.686, 0.703, 0.712 for $j=1$ and 1.685, 1.706, 1.717 for $j=2$. The first version printed 0.687 and 0.713 for $j = 1$; the values are 0.68558 and 0.71240.

**Check 2 — the orbit profile.** Exactly, $\mathbb{E}[\#\text{orbits of size }m] = \binom{n}{m}a_m h_k(n-m)/h_k(n)$. At $k=n/4$ the expected fraction of points in orbits of size 4 is

| $n$ | 64 | 128 | 256 | 512 |
| --- | --- | --- | --- | --- |
| size 4 | 0.735008 | 0.804914 | 0.857902 | 0.897348 |
| size 2 | 0.264883 | 0.195086 | 0.142098 | 0.102652 |
| size $\ge 8$ | $1.0\times10^{-4}$ | $4.0\times10^{-9}$ | $2.2\times10^{-18}$ | $2.6\times10^{-37}$ |

rising toward 1 with the residual entirely on size 2, exactly as the mixture argument predicts. (`src/final_check.py` verifies $\sum_m m\,\mathbb{E}[\#\text{orbits}_m]=n$ as an arithmetic self-check to within $2\times10^{-9}$; the largest error is $1.4\times10^{-9}$, at $n = 512$, from rounding in the float logarithms that `ratio` uses. The first version said $10^{-9}$, which the check did not enforce: it allowed $10^{-6}n$.)

The first version printed $1.1\times10^{-4}$ for size $\ge 8$ at $n = 64$; the value is $1.049\times10^{-4}$.

**Check 3 — out-of-sample extrapolation.** A direct fit is not enough: $A$ converges slowly, and a basis chosen by in-sample fit proves nothing when seven candidate shapes are available. So $A_{\text{local}}(n)$ was measured by second differences (which annihilate the $Bn+C$ terms exactly), fitted on $n\in[96,296]$ only, and scored on $n\in[400,496]$, which no fit ever saw:

| basis | $A_\infty$ | train max-res | **held-out max-res** |
| --- | --- | --- | --- |
| $n^{-1/2}$ | **0.74962** | 3.3e−5 | **6.7e−5** |
| $1/\log_2^2 n$ | 0.75598 | 1.7e−4 | 2.9e−4 |
| $1/\log_2 n$ | 0.78532 | 1.1e−4 | 5.0e−4 |
| $n^{-1/4}$ | 0.77239 | 3.3e−4 | 1.1e−3 |
| $\log_2 n/n$ | 0.74102 | 3.9e−4 | 1.1e−3 |
| $n^{-1/2}\log_2 n$ | 0.76416 | 4.5e−4 | 1.2e−3 |
| $n^{-1}$ | 0.73823 | 7.3e−4 | 1.9e−3 |

Scanning a free exponent $p$ in $A_\infty + cn^{-p}$ the same way selects $p = 0.48$ and $A_\infty = \mathbf{0.75057}$, within $5.7\times10^{-4}$ of the predicted $3/4$, held-out max residual $1.6\times10^{-5}$.

The first version of the table printed the held-out residuals of $n^{-1/2}$ and $n^{-1/2}\log_2 n$ as 6.8e−5 and 1.3e−3; the values are $6.749\times10^{-5}$ and $1.248\times10^{-3}$. The ranking is unchanged.

**A rejected hypothesis, kept.** The first version of the model predicted $\varepsilon \sim c/\log_2 n$ for the fraction $\varepsilon$ of points in 2-orbits, hence $A_{\text{local}} = 3/4 - c/\log_2 n$. That is **false**: measured $\varepsilon\log_2 n$ falls 1.589 → 0.924 over $n=64..512$ while $\varepsilon\sqrt{n}$ only drifts 2.119 → 2.323, and the $1/\log_2 n$ basis extrapolates to 0.785, missing the target by 0.035 despite a good in-sample fit. The first version called $\varepsilon$ the fraction of points *off* 4-orbits. The numbers are the 2-orbit fraction; the off-4 fraction also counts fixed points and orbits of size $\ge 8$, and gives 1.590 and 2.120 at $n = 64$ and the same values as above at 512. The exponent $1/2$ is selected by held-out error. The first version said *"I have no derivation of why it is $1/2$"*; §9.2 now gives one, from the saddle point of the set-partition count $M(n)$, which predicts $\varepsilon\sqrt n \to \sqrt6 = 2.449$ and the coefficient $-\sqrt6/8$.

## 7. The gap in [RoTr25] Theorem 6, located exactly

Theorem 6 is the only published result addressing the second half of #1162:

> "Let µ be a real number in [0, 1/16) and let ν be a real number in [0, ½ − √3⁄4). Then a random subgroup of S_n has an elementary abelian 2-section of order at least 2^{µn}, and a Sylow 2-subgroup of order at least 2^{νn}."

Both constants come from the same place. From the proof (p. 35):

> "use Lemma 2.4(i) to see that the number of R-subgroups of a Sylow 2-subgroup of S_n is at most 4·2^{νn(n−νn)} … Since ν(1 − ν) < 1/16, the result follows"

and for the section, "Since µ < 1/16 the result follows from the lower bound in Theorem 1". Solving $\nu(1-\nu) = 1/16$ gives $\nu = \tfrac12 - \tfrac{\sqrt3}{4} = 0.06698729810778\ldots$, which is their printed constant exactly — so this is the binding step, not a by-product. The first version said this reproduced the constant "to 25 digits"; the paper prints the closed form $\tfrac12 - \tfrac{\sqrt3}{4}$, so the agreement is algebraic and a digit count adds nothing.

**The lossy step.** Lemma 2.4(i) is applied to a Sylow 2-subgroup of $S_n$, of order $2^{\ell}$ with $\ell = n - s_2(n) \approx n$. The bound $\zeta_2 2^{k(\ell-k)}$ is attained only when the ambient group is elementary abelian of rank $\ell$. But the maximal elementary abelian 2-subgroups of $S_n$ have rank $\lfloor n/2\rfloor$, not $n$. Replacing $\ell \approx n$ by $n/2$ turns the threshold equation into

$$\nu\left(\tfrac12 - \nu\right) = \tfrac{1}{16},$$

whose discriminant is **exactly zero**, with a double root at $\nu = \tfrac14$. The same substitution sends µ's threshold $\mu < 1/16$ to $\mu(\tfrac12-\mu) < \tfrac1{16}$, again a double root at 1/4.

That the barrier is a *double* root rather than an interval endpoint is the signature of a tight extremal problem, and $1/4$ is precisely the rank at which §3 finds the elementary abelian family concentrated.

**Conjecture 7.1.** Theorem 6 holds for all $\mu < 1/4$ and all $\nu < 1/4$, and 1/4 is optimal for both.

Improvement factors: $4.000$ for µ, $3.732$ for ν.

**Conjecture 7.2 (the statistical theorem Erdős asked for).** For a uniformly random $H \le S_n$ with $n \equiv 0 \pmod 4$,
$$\log_2|H| = \left(\tfrac14 + o(1)\right)n \quad\text{in probability.}$$

Posed along $n \equiv 0 \pmod 4$ deliberately: [RoTr25] Theorem 4 disproves Kantor's conjecture only for $n \equiv 3 \pmod 4$ ("Let n be congruent to 3 modulo 4 … the probability that a random subgroup of $S_n$ is nilpotent is bounded away from 1"), so $n \equiv 0 \pmod 4$ is the residue class where a random subgroup may still be nilpotent, hence a 2-group by Theorem 5 ("A random nilpotent subgroup of $S_n$ is a 2-group").

The $\ge$ direction of 7.2 is Conjecture 7.1. The $\le$ direction is not addressed by anything in [RoTr25].

**Status after §9 (September 22, 2026).** Conjecture 7.2 is withdrawn, not disproved. Its only support was the elementary abelian family, and §9.5 proves that family is a vanishing fraction of all subgroups. A larger family, counted exactly in §9, concentrates at $\log_2|H| \approx 3n/8$ instead (§9.6), and Conjecture 9.3 replaces 7.2. Of Conjecture 7.1, the claim that 1/4 is optimal for ν is withdrawn: under Conjecture 9.3 a random subgroup is, for even $n$, a 2-group of order about $2^{3n/8}$, so ν could be taken up to 3/8. The double root is no evidence for 1/4 either. The same substitution applied to that family gives a double root at 3/8, and §9.9 shows a double root appears for any family whose count reaches $2^{n^2/16}$, so the paragraph above that calls it "the signature of a tight extremal problem" is true but says nothing about which order is typical. Nothing here bears on µ.

## 8. Evidence against, and what is not established

**The elementary abelian family may not be typical.** This is the load-bearing gap. $\log e(n)$ and $\log f(n)$ agree to $o(n^2)$, but Theorem 1's upper bound permits $f(n)/e(n)$ as large as $2^{\beta n^{3/2}}$, so elementary abelian subgroups could be a vanishing fraction of all subgroups forever. Everything in §3 is then a theorem about a measure-zero family, and Conjecture 7.2 would not follow. **Nothing here rules that out.** *Added September 22, 2026:* it is now proved to happen. §9.5 shows $e(n)/f(n) \le 2^{-\frac18 n\log_2 n + O(n)}$, so §3 describes a vanishing family, and Conjecture 7.2 is withdrawn (§7, last paragraph).

**The small-$n$ data points the other way.** From the exact lattices:

| $n$ | 4 | 5 | 6 | 7 | 8 |
| --- | --- | --- | --- | --- | --- |
| $e(n)/f(n)$ | 0.4667 | 0.2949 | 0.1863 | 0.1165 | **0.0809** |
| median $\|H\|$ | 4 | 6 | 8 | 12 | 16 |
| $\mathbb{E}[\log_2\|H\|]/n$ | 0.474 | 0.491 | 0.507 | 0.500 | **0.4999** |

The elementary abelian fraction is *falling*, and the mean of $\log_2|H|$ sits at $n/2$, not $n/4$. This is not evidence for Conjecture 7.2; if anything it runs against it. The honest reading is that $n \le 8$ is nowhere near the asymptotic regime — the $n^2/16$ term only overtakes competing structures for $n$ far beyond anything enumerable — and that small-$n$ lattice data simply cannot discriminate here. It is recorded because omitting it would be misleading.

The larger family of §9 is closer to these numbers, though not close. At $n = 8$ the Frattini-closed family is 40.9% of all subgroups, against 8.1% for the elementary abelian family, and its mean $\mathbb{E}[\log_2|H|]/n$ is 0.4388, against 0.2947 for the elementary abelian family and 0.4999 for all subgroups. The same caveat applies: $n = 8$ is not asymptotic, and the §9 prediction, $3/8$, is approached from above (§9.6).

**Not proved.** The limits in §3 and §4, and the rate of approach in §6, are established as exact computations at each tested $n$ and as heuristic derivations, not as theorems. The headline of §6, $\log_2 e(n) = n^2/16 + \tfrac34 n\log_2 n + O(n)$, was in the same state until §9.5 proved it. No saddle-point analysis of $a(n,k)$ was carried out; §9.2 carries one out, heuristically, for the simpler count $M(n)$. Theorem 9.1 and Propositions 9.2, 9.4 and 9.5 are the only asymptotic statements proved here. Conjectures 7.1, 7.2 and 9.3 are conjectures, and 7.2 and the ν-optimality half of 7.1 are withdrawn (§7).

**Novelty not established.** A search of arXiv (metadata, not full text), the Semantic Scholar citation list of 2503.05416, and a web search found no statement of the $n/4$ constant. That citation list grew from zero to three between the first check and the final one, which is a useful reminder that this kind of negative result has a short shelf life. The three, re-fetched at the end, are: *Proportion of Simple Subgroups in Finite Groups and Their Applications* (arXiv:2606.17488), which studies $\mathcal{V}(G)=\mathrm{Simp}(G)/|L(G)|$ — a different statistic on the same lattice, and not an asymptotic for the order of a random subgroup; *Homological Nielsen realization for the manifolds $\#_n\mathbb{CP}^2$* (arXiv:2605.27537), unrelated; and *Groups having 12 cyclic subgroups* (arXiv:2210.11788), which predates 2503.05416 and is presumably a Semantic Scholar linkage artifact. None overlaps §7. Absence of search results is weak evidence regardless; the arXiv `all:` field does not index full text. Re-run `refs/fetch.sh` before relying on any of this.

**Novelty of §9 not established either.** Theorem 9.1 and Propositions 9.2 and 9.4 were checked against [RoTr25] alone. It does not state them, and the strings "7/8", "0.875" and "extraspecial" do not occur in its text. It does isolate the same group. By its Proposition 5.3(i) a transitive 2-group of degree 8 has excess 2 exactly when $d(G) = 4$, which by §9.1 means 8T22, and its upper-bound proof treats that case separately (Theorem 5.7(i), "The group G has excess 2 and degree 8", and Theorem 5.5(ii)(b), "If no Gi, for i ≤ r, has excess 2 and degree 8"). Its lower bound, Proposition 7.3, counts subdirect products of a group "generated by m disjoint p-cycles" instead. So Theorem 9.1 is the lower-bound side of a case [RoTr25] already singled out, and may well be known to its authors. No other literature was searched for §9.

**Open questions raised here.** (i) Why is the approach exponent in §6 equal to $1/2$? Answered heuristically in §9.2 by the saddle point of $M(n)$, which also predicts the coefficient; not proved. (ii) Is $e(n)/f(n) \to 0$, and if so how fast? Answered: yes, and $e(n)/f(n) \le 2^{-\frac18 n\log_2 n + O(n)}$ (§9.5); the true rate is not known. (iii) Can the $\vartheta$-law of §3 be proved by a saddle-point argument? Still open; §9.3 reduces it to (vi). *Added October 19, 2026:* proved, by a direct count rather than a saddle point, as Proposition 9.5(c) and (d) in §9.10. (iv) Is Conjecture 9.3 true? (v) Are the 2-subgroups that are not Frattini-closed negligible? At $n = 8$ they are 3150 of 65 054, or 4.8%. *Added November 1, 2026:* the share rises to 16.5% at $n = 14$ (§9.8), the opposite of what premise (b) needs. (vi) Prove $e(n) = G_{\lfloor n/2\rfloor}M(n)(1+o(1))$ for even $n$, and the analogue for odd $n$ (§9.2). *Added October 19, 2026:* proved as Proposition 9.5(a) in §9.10, with an explicit rate far weaker than the data show.

## 9. A larger family: Frattini-closed 2-subgroups (added September 22, 2026)

This section mixes three grades of claim, and each is marked where it occurs. **Proved:** the facts in §9.1, Theorem 9.1, and Propositions 9.2, 9.4 and 9.5 (§9.10, added October 19, 2026). The proofs are short and given in full; `src/final_check.py` checks their arithmetic, not their logic. **Computed exactly:** the family count $B(n)$ and its order statistics for every $n \le 512$, and, for $n \le 15$, the number of all 2-subgroups of $S_n$ by order and by Frattini-closure (§9.8, added November 1, 2026). **Heuristic:** every limit and rate in §9.2, §9.3, §9.6 and §9.7, and Conjecture 9.3. The computations are `src/census8.py` (brute force on permutations for $n \le 8$, 33 checks) and `src/block_families.py` (62 gates), and for the census of §9.8 `src/two_subgroups.py` (122 checks) with its independent recount `src/verify_two_subgroups.py`. All pass, and all were tested by planting defects (§10).

### 9.1 The family

Let $H \le S_n$ be a 2-group with orbits $O_1, \dots, O_t$, and let $P_i$ be the transitive 2-group that $H$ induces on $O_i$. Call $H$ *Frattini-closed* (FC) if it contains $N = \prod_i \Phi(P_i)$. Elementary abelian subgroups are FC, since every $\Phi(P_i)$ is trivial. $\langle(1234)(5678)\rangle$ is not: it contains $(13)(24)(57)(68)$ but not $(13)(24)$. For FC $H$, $V = H/N$ is a subspace of $\prod_i P_i/\Phi(P_i) = \mathbb{F}_2^D$, $D = \sum_i d(P_i)$, that projects onto every factor. Conversely, the preimage in $\prod_i P_i$ of such a $V$ is FC with the same orbits and the same $P_i$, because a subgroup of $P_i$ that maps onto $P_i/\Phi(P_i)$ is $P_i$. So the FC subgroups with given orbit data correspond one-to-one with these $V$, and $|H| = 2^{\dim V}\prod_i|\Phi(P_i)|$.

Möbius inversion on the subspace lattice counts the $V$ that project onto every factor as $L\big(\prod_i w_{d(P_i)}\big)$. Here $w_d(z) = \prod_{i<d}(z - 2^i)$ is the block weight, by the $q$-binomial theorem, and $L$ is the linear map $z^D \mapsto G_D$. The exponential formula then counts $B(n)$, the FC subgroups of $S_n$ whose orbits all have size at most 8:

$$B(n) = n!\,[x^n]\; L\Big(\exp \sum_{b \in \{1,2,4,8\}} W_b(z)\,\frac{x^b}{b!}\Big), \qquad W_b = \sum_P w_{d(P)},$$

where $P$ runs over the transitive 2-subgroups of $S_b$. The same recurrence, with the first two moments of $\log_2|\Phi|$ carried along, gives the exact mean and variance of $\log_2|H|$ over the family. Distinct orbit data give distinct subgroups, so every $B(n)$ is a rigorous lower bound for $f(n)$.

`src/census8.py` finds every 2-subgroup of $S_8$ by brute force on permutations: 65 054 of them, of which 61 904 are FC. The 3150 that are not have orders 4 (630 of them) and 8 (2520). For $n \le 7$ every 2-subgroup is FC (1, 2, 4, 20, 76, 631 and 3417 of them for $n = 1, \dots, 7$), and the counts by order match the lattices of §2.1 for every $n \le 8$. The transitive 2-subgroups of $S_b$ number 1, 1, 7 and 12 945 for $b = 1, 2, 4, 8$. The 26 conjugacy classes in degree 8 have orders 8 (5 classes), 16 (6), 32 (8), 64 (6) and 128 (1), the order profile of 8T1–8T5, 8T6–8T11, 8T15–8T22, 8T26–8T31 and 8T35 in [GN]. Every one has $d(P) \le b/2$; this is read off the census, not taken from the literature. Gate 1 of `block_families.py` recomputes $B(n)$, its order histogram and the first two moments of $\log_2|H|$ from the census types by the formula above, and all three equal the brute-force values for every $n \le 8$.

Only types with $d(P) = b/2$ can give $D = n/2$. There are 1, 4 and 105 of them in degrees 2, 4 and 8: $C_2$; the regular $V_4$ and the three $D_4$; and the 105 conjugates of one group of order 32, 8T22. Their exponential generating function is $\exp(x^2/2 + x^4/6 + x^8/384)$, and its coefficients $B^*(n)$ play the part for $B$ that $M(n)$ plays for $e$ in §9.2.

8T22 is the extraspecial group $2^{1+4}_+ = D_4 \circ D_4$, which is how [GN] lists it (SmallGroup(32,49)). The proof uses only the census invariants. It has order 32 and $d = 4$, so $|\Phi| = 2$. A normal subgroup of order 2 is central, and $P/\Phi$ is elementary abelian, so $P' \le \Phi \le Z(P)$. $P$ is not abelian, because an abelian transitive group is regular and would have order 8; so $P' = \Phi$. The census gives $|Z(P)| = 2$, so $Z(P) = P' = \Phi(P)$ has order 2 and $P$ is extraspecial. Of the two extraspecial groups of order 32, $2^{1+4}_+$ has $2^4 + 2^2 - 1 = 19$ involutions and $2^{1+4}_-$ has $2^4 - 2^2 - 1 = 11$. The census counts 19.

### 9.2 $e(n)$ is close to $G_{n/2}M(n)$ (numerical)

Let $M(n) = n!\,[x^n]\exp(x^2/2 + x^4/24)$, the number of partitions of $[n]$ into blocks of sizes 2 and 4. Such a partition, with a transposition on each 2-block and the regular $V_4$ on each 4-block, has $D = n/2$, which suggests $e(n) \approx G_{n/2}M(n)$ for even $n$. For odd $n$ the model is $G_{(n-1)/2}\,n\,M(n-1)$, with one fixed point. The exact residuals:

| $n$ | 8 | 16 | 32 | 64 | 128 | 256 | 512 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $\log_2 e(n) - \log_2 G_{n/2}M(n)$ | −0.938 | −0.566 | −0.0412 | 5.2e−4 | 7.9e−8 | 1.0e−16 | 2.7e−35 |

At $n = 9, 17, 33$ the odd model gives −1.378, −0.765 and −0.0657. The table is computed from the exact integers. Over the 49 values $n \ge 128$ stored in `data/elemab_stats.json`, whose doubles cannot resolve residuals below about $10^{-12}$, the largest residual in absolute value is $7.9\times10^{-8}$, at $n = 128$. This is numerical evidence for $e(n) = G_{\lfloor n/2\rfloor}M(n)(1 + o(1))$, open question (vi), not a proof. Proposition 9.5 in §9.10 now proves the statement, with a rate of $2^{-n/64}$ against the $2^{-n/4}$ seen here.

The saddle point of $M$ explains the exponent $1/2$ of §6. At the saddle, $x^2 + x^4/6 = n$, the expected number of points in 2-blocks is $u = x^2 = \sqrt{6n+9} - 3$, so the 2-orbit fraction is $\varepsilon = u/n \approx \sqrt{6/n}$. That gives 0.26288, 0.19433, 0.14182 and 0.10255 at $n = 64, 128, 256, 512$, against the exact 0.264883, 0.195086, 0.142098 and 0.102652 of §6 Check 2. So $\varepsilon\sqrt n \to \sqrt6 = 2.449$, slowly: the exact value is 2.323 at $n = 512$. Carrying $\log_2 M$ one term further gives $A(n) = \tfrac34 - \tfrac{\sqrt6}{8}n^{-1/2} + \tfrac{9\sqrt6}{32}n^{-3/2} + O(n^{-2})$ for the local coefficient of §6 Check 3, with no $n^{-1}$ term. Three checks agree. The §6 fit found $-0.298$ for the coefficient of $n^{-1/2}$, against $-\sqrt6/8 = -0.306$. $A_{\text{local}}$ computed from $G_{n/2}M(n)$ equals the stored values computed from $e(n)$ to $1.8\times10^{-7}$ for $n = 128, \dots, 504$. And $\big(A_{\text{local}}(n) - \tfrac34 + \tfrac{\sqrt6}{8}n^{-1/2}\big)n^{3/2}$, computed from exact $M$ at $n = 512, 1024, \dots, 32768$, runs 0.62799, 0.64695, 0.65982, 0.66864, 0.67473, 0.67897, 0.68192, with Richardson limit 0.68905 against $9\sqrt6/32 = 0.68892$. The saddle-point step itself is not justified here.

### 9.3 The rank law is that of a random subspace (numerical)

Let $m = \lfloor n/2\rfloor$. The total variation distance between the rank law $a(n,\cdot)/e(n)$ and the law of the dimension of a uniformly random subspace of $\mathbb{F}_2^m$, $\binom{m}{k}_2/G_m$, is:

| $n$ | 16 | 32 | 64 | 128 | 256 |
| --- | --- | --- | --- | --- | --- |
| TV | 9.97e−2 | 7.12e−3 | 8.40e−5 | 1.28e−8 | 1.67e−17 |
| $\log_2(\mathrm{TV}\cdot 2^{n/4}/n^2)$ | −7.33 | −9.13 | −9.54 | −8.22 | −7.73 |

The rate is argued, not proved. The largest correction to $e \approx G_{n/2}M$ at the saddle is a single regular $(\mathbb{Z}/2)^3$ orbit of size 8. It has $8!/(8\cdot168) = 30$ embeddings, so its weight is $30x^8/8! = x^8/1344$ with $x^8 = u^4 \approx 36n^2$, and it lowers $D$ by 1, which costs $G_{m-1}/G_m \approx 2^{-n/4}$. That predicts a TV of order $n^2 2^{-n/4}$. The last row stays between −9.54 and −7.33 while the TV falls by more than 15 orders of magnitude, but it is not constant, and the argument does not predict it ($\log_2(36/1344) = -5.2$). The same row for the sd error of §3 is −10.37, −8.91, −8.41 and −8.15 at $n = 64, 128, 256, 512$. If the TV bound were proved, the $\vartheta$-law of §3 would follow from $\binom{m}{m/2+j}_2/\binom{m}{m/2}_2 \to 2^{-j^2}$, so open question (iii) reduces to (vi). Proposition 9.5(c) in §9.10 proves that the TV tends to $0$, at the weaker rate $2^{-n/64}$, and (d) draws the $\vartheta$-law from it.

### 9.4 Theorem 9.1: the $n\log n$ coefficient of $f$ is at least $7/8$

**Theorem 9.1.** $\log_2 f(n) \ge \frac{n^2}{16} + \frac78 n\log_2 n - O(n)$. The subgroups counted are 2-groups, so the same bound holds for the number of 2-subgroups and for the number of nilpotent subgroups of $S_n$.

*Proof.* Let $n = 8m$. Partition $[n]$ into $m$ blocks of 8, put a conjugate of 8T22 on each block, and let $H$ range over the FC subgroups with these orbit data. There are $\frac{(8m)!}{m!\,8!^m}105^m$ such orbit data, and they give distinct subgroups, since $H$ determines its orbits and the groups it induces on them. For each, the number of $H$ is the number $L(w_4^m)$ of subspaces of $\mathbb{F}_2^{4m}$ that project onto each of the $m$ factors $\mathbb{F}_2^4$. A subspace that misses factor $i$ lies in the preimage of one of the 15 hyperplanes of that factor, which is a hyperplane of $\mathbb{F}_2^{4m}$, so $L(w_4^m) \ge G_{4m} - 15m\,G_{4m-1}$. Summing the $q$-Pascal rule $\binom{r}{k}_2 = \binom{r-1}{k-1}_2 + 2^k\binom{r-1}{k}_2$ over $k$ gives $G_r = G_{r-1}(1 + \mathbb{E}\,2^K)$, where $K$ is the dimension of a uniformly random subspace of $\mathbb{F}_2^{r-1}$. By the symmetry $\binom{r-1}{k}_2 = \binom{r-1}{r-1-k}_2$, $\mathbb{E}K = (r-1)/2$, so Jensen's inequality gives $G_r \ge 2^{(r-1)/2}G_{r-1}$. Hence $L(w_4^m) \ge G_{4m}\big(1 - 15m\,2^{-(4m-1)/2}\big)$, and the bracket is positive for $m \ge 3$. Finally $G_{4m} \ge \binom{4m}{2m}_2 \ge 2^{4m^2} = 2^{n^2/16}$, and Stirling's formula gives $\log_2\frac{(8m)!\,105^m}{m!\,8!^m} = \frac78 n\log_2 n + cn + O(\log n)$ with $c = \frac38 - \frac78\log_2 e - \frac18\log_2 384 = -1.96048$. For $n = 8m + r$ with $0 < r < 8$, fix $r$ points; $(n-r)^2/16 \ge n^2/16 - 7n/8$. ∎

The bracket is 0.00563, 0.669 and 0.896 at $m = 3, 4, 5$, against true values of $L(w_4^m)/G_{4m}$ of 0.469, 0.765 and 0.918. The count in the proof is exact, not only a bound: `block_families.py` computes $F_8(n) = \frac{n!}{m!\,8!^m}105^m L(w_4^m)$ for every $n = 8m \le 512$. Subtracting $\frac{n^2}{16} + \frac78 n\log_2 n + cn + \frac32 + \log_2\frac{\vartheta}{\eta} - \frac{7\log_2 e}{12n}$ from $\log_2 F_8(n)$ leaves $-3.2\times10^{-3}$, $8.8\times10^{-7}$, $1.2\times10^{-7}$ and $1.5\times10^{-8}$ at $n = 64, 128, 256, 512$. From $n = 256$ on this equals the next Stirling term, $\frac{\log_2 e}{360m^3}(1 - 8^{-3})$, to three digits; at $n = 64$ it is dominated by $\log_2(L(w_4^m)/G_{4m})$. $F_8(n)$ first exceeds $e(n)$ at $n = 24$, and $\log_2(F_8(n)/e(n))$ is 4.6, 21.4, 64.8, 170.8 and 419.6 at $n = 32, \dots, 512$. The approach to $7/8$ is slow: $(\log_2 F_8(n) - n^2/16)/(n\log_2 n)$ is 0.658 at $n = 512$, because the linear term costs $1.96/\log_2 n = 0.218$ there.

Consequences for [RoTr25]:

* (a) Theorem 1 is stated "for all integers n > 1", and $n = 2$ alone forces $\alpha \le 3/8$ there, so Theorem 9.1 does not let one put $\alpha = 7/8$ in it. What it improves is the constant for large $n$: $\liminf_n (\log_2 f(n) - n^2/16)/(n\log_2 n) \ge 7/8$. The paper takes $\alpha$ from its Proposition 7.3 ("The lower bound is Proposition 7.3"), which sets $\alpha_p = \min\{0.08,\ (1-(2p-1)^2/(4p^2))/(2p\log(2p)),\ \alpha_R,\ (1-\alpha_0)(1-1/p)\}$ with $\alpha_0 = 448/453$, read from a render of page 23. At $p = 2$ the last entry is $5/906$, so the explicit constant is at most $5/906 < 0.0056$; the second entry is $7/128$, since the paper's logarithms are to base 2. That construction counts subgroups "generated by m disjoint p-cycles", whose natural count is $n^2/16 + \frac12 n\log_2 n$ (the $j = 1$ row of §6).
* (b) The upper-bound constants of its Theorems 2 and 3 must satisfy $\beta_2 \ge 7/8$ and $\gamma(2) \ge 7/8$.
* (c) Any $\gamma$ that answers its first question, $|\mathrm{Sub}(S_n)| \le 2^{n^2/16 + \gamma n\log n}$, is at least $7/8$.

### 9.5 Proposition 9.2: elementary abelian subgroups are a vanishing fraction

**Proposition 9.2.** $\log_2 e(n) = \frac{n^2}{16} + \frac34 n\log_2 n + O(n)$. With Theorem 9.1, $e(n)/f(n) \le 2^{-\frac18 n\log_2 n + O(n)} \to 0$.

*Proof of the upper bound.* An elementary abelian $E \le S_n$ acts regularly on each of its orbits, so they have sizes $2^{j_i}$, and $E$ is a subspace of $\mathbb{F}_2^D$ with $D = \sum_i j_i$. Hence $e(n) \le \sum G_D$, summed over orbit data. Let $n_0$ be the number of fixed points and $n_+$ the number of points in orbits of size at least 8. A point in an orbit of size $2^j$ contributes $j/2^j$ to $D$, which is $1/2$ for $j = 1, 2$ and at most $3/8$ for $j \ge 3$. So $D = n/2 - \delta$ with $\delta \ge n_0/2 + n_+/8$, and since $\delta \le n/2$, $D^2/4 \le n^2/16 - n\delta/8$. From $\binom{D}{k}_2 \le 2^{k(D-k)}/\eta$ and $\sum_k 2^{k(D-k)} \le \vartheta\,2^{D^2/4}$, $G_D \le (\vartheta/\eta)\,2^{D^2/4} < 8\cdot2^{D^2/4}$. To count the orbit data, label each point as fixed, in a block of size 2 or 4, or in a block of size at least 8 ($3^n$ ways). Partition the second kind into 2- and 4-blocks, each with its unique regular group ($M(n')$ ways, $n' = n - n_0 - n_+$). List the third kind in order and cut the list into consecutive blocks, each identified with $\mathbb{F}_2^j$ by position (at most $n_+!\,2^{n_+} \le n^{2n_+}$ ways). So

$$e(n) \le 8\cdot3^n\,2^{n^2/16}\max_{n_0,\,n_+} M(n')\,n^{2n_+}\,2^{-n n_0/16 - n n_+/64}.$$

For $n \ge 2^{11}$, $n^2 \le 2^{n/64}$, so $e(n) \le 8\cdot3^n\,\bar M(n)\,2^{n^2/16}$ with $\bar M(n) = \max_{n' \le n} M(n')$. Cauchy's bound at $x = n'^{1/4}$ gives $M(n') \le n'!\,e^{x^2/2 + x^4/24}x^{-n'}$, so $\log_2\bar M(n) \le \frac34 n\log_2 n + O(n)$. ∎

*Proof of the lower bound.* Run the construction of Theorem 9.1 with 4-blocks: partition $[4t]$ into $t$ blocks of 4, each carrying the regular $V_4$. The $(4t)!/(t!\,24^t)$ partitions give distinct subgroups, and each carries $L(w_2^t) \ge G_{2t}\big(1 - 3t\,2^{-(2t-1)/2}\big)$ elementary abelian subgroups, with 3 hyperplanes per factor. The bracket is positive for $t \ge 5$: it is −0.061, 0.337 and 0.602 at $t = 4, 5, 6$, against true ratios $L(w_2^t)/G_{2t}$ of 0.459, 0.601 and 0.729. Then $G_{2t} \ge 2^{t^2} = 2^{n^2/16}$ and $\log_2\frac{(4t)!}{t!\,24^t} = \frac34 n\log_2 n + O(n)$. For $n \not\equiv 0 \pmod 4$, fix up to 3 points as in Theorem 9.1. ∎

So the headline of §6 is a theorem, with the $O(n)$ term unquantified. The true rate at which $e(n)/f(n) \to 0$ is not known; the proofs bound it by $2^{-\frac18 n\log_2 n + O(n)}$. For comparison, the exact $\log_2(B(n)/e(n))$ is 6.2, 14.7, 35.8, 85.7, 200.8 and 462.1 at $n = 16, 32, \dots, 512$.

### 9.6 Order statistics in the family (exact for $n \le 512$; limits heuristic)

Exact for each $n$, with $f_b$ the expected fraction of points in orbits of size $b$:

| $n$ | $\mathbb{E}\log_2\|H\|/n$ | sd | sd$/n^{1/4}$ | $f_1$ | $f_2$ | $f_4$ | $f_8$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | 0.4388 | 1.050 | 0.625 | 0.055 | 0.282 | 0.454 | 0.209 |
| 16 | 0.4162 | 1.293 | 0.646 | 0.013 | 0.208 | 0.473 | 0.306 |
| 32 | 0.3918 | 1.413 | 0.594 | 7.1e−4 | 0.152 | 0.491 | 0.357 |
| 64 | 0.3885 | 1.584 | 0.560 | 1.8e−6 | 0.099 | 0.413 | 0.489 |
| 128 | 0.3875 | 1.863 | 0.554 | 1.7e−11 | 0.062 | 0.323 | 0.615 |
| 256 | 0.3856 | 2.203 | 0.551 | 2.4e−21 | 0.038 | 0.245 | 0.717 |
| 512 | 0.3835 | 2.610 | 0.549 | 8.1e−41 | 0.023 | 0.182 | 0.795 |

The mean falls toward $3/8$, the value when every point lies in an 8T22 orbit: $\frac14$ from $\dim V = D/2$, plus $\frac18\log_2|\Phi|$ per point. The sd grows like $n^{1/4}$. Two exact quantities tie the table to the dense types of §9.1. First, $\mathbb{E}\log_2|H|/n$ is close to $\frac14 + \frac3{16}f_4 + \frac18 f_8$, which is its value when every orbit type is dense and $\dim V$ averages $D/2$ (the four dense 4-types have $|\Phi| = 1$ once and $|\Phi| = 2$ three times). The difference, computed exactly, is 0.078 at $n = 8$, 0.0052 at $n = 32$, $5.2\times10^{-10}$ at $n = 128$, $1.4\times10^{-19}$ at $n = 256$ and $8.2\times10^{-39}$ at $n = 512$. Second, $\log_2 B(n) - \log_2 G_{n/2} - \log_2 B^*(n)$ is $5.1\times10^{-8}$, $2.9\times10^{-17}$ and $3.7\times10^{-36}$ at $n = 128, 256, 512$.

The limits come from the saddle point of $B^*$ and are heuristic. There $y = x^2$ solves $y + \frac23y^2 + \frac{y^4}{48} = n$, and the expected numbers of points in orbits of size 2, 4 and 8 are $y$, $\frac23y^2$ and $\frac{y^4}{48}$. At $n = 512$, $y = 11.8225$ predicts $f_2, f_4, f_8 = 0.0231, 0.1820, 0.7949$, against the exact 0.0231, 0.1823, 0.7945. Substituting into the identity gives $\mathbb{E}\log_2|H| - \frac{3n}{8} \approx \frac{y^2}{24} - \frac{y}{8}$: 0.845, 1.586, 2.701 and 4.346 at $n = 64, 128, 256, 512$, against the exact 0.864, 1.600, 2.713 and 4.355. For the spread, Gaussian conditioning at the saddle gives $\mathrm{sd}^2 \approx \frac{y^2}{24} + \frac{y}{32} - \frac16 + 0.7213$, where the last term is the variance of the §3 law. That gives an sd of 1.541, 1.835, 2.184 and 2.598, against the exact 1.584, 1.863, 2.203 and 2.610, a shortfall that falls from 0.042 to 0.012. As $n \to \infty$ the heuristic gives $\log_2|H| = \frac{3n}{8} + \big(\frac{1}{\sqrt{12}} + o(1)\big)\sqrt n + O_p(n^{1/4})$ inside the family, with sd $\sim 48^{1/4}24^{-1/2}\,n^{1/4} = 0.5373\,n^{1/4}$. At $n = 512$ the mean is still $0.0085\,n$ above $3n/8$.

### 9.7 Odd $n$, and the other two questions of [RoTr25]

For odd $n$ a 2-group has a fixed point, which costs half a dimension, and a second family takes over. Let $S_3(n)$ be the set of $H \le S_n$ for which, for some 3-set $T$, the group $A_3$ on $T$ is a normal subgroup of $H$ and $H/A_3$ is FC with orbits of size at most 8. Here $H/A_3$ acts on $T$ through $S_3/A_3 = C_2$, as on a block of size 2. Taking the quotient by $A_3$ is a bijection onto such FC groups. The $T$-block enters with weight $z$ rather than $w_1(z)$, because $H$ may act on $T$ as $A_3$ or as $S_3$. So $|S_3(n)| = \binom n3 L(z\,F_{n-3})$, where $F_{n-3}$ is the polynomial whose image under $L$ is $B(n-3)$. $T$ is unique, since $H/A_3(T)$ is a 2-group, and $S_3(n)$ is disjoint from $B(n)$, whose members are 2-groups. `block_families.py` checks $|B(n)| + |S_3(n)| \le f(n)$ against A005432 for $n \le 18$. Along odd $n$, $S_3$ overtakes $B$:

| $n$ | 9 | 17 | 33 | 65 | 129 | 257 | 511 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $\log_2(\|S_3(n)\|/B(n))$ | −0.872 | −0.648 | −0.288 | 0.072 | 0.396 | 0.696 | 0.979 |
| $\log_2(y/6)$ | −1.124 | −0.687 | −0.291 | 0.066 | 0.391 | 0.693 | 0.978 |

The second row is the saddle-point estimate $S_3(n)/B(n) \approx \binom n3 B^*(n-3)/\big(nB^*(n-1)\big) \approx y/6$, with $y$ from §9.6 at the same $n$. Both families carry the Galois factor $G_{(n-1)/2}$, and $B^*(n-3)/B^*(n-1) \approx y/n^2$. Along even $n$, $S_3$ needs a fixed point and loses a factor of about $2^{n/4}$: $\log_2(|S_3(n)|/B(n))$ is −1.896, −2.911, −5.762, −13.027, −28.379, −59.776 and −123.205 at $n = 8, 16, \dots, 512$, which is −0.948, −0.728, −0.720, −0.814, −0.887, −0.934 and −0.963 per $n/4$.

All but a fraction $2^{-\Theta(n)}$ of $S_3(n)$ act on $T$ as $S_3$, so they are not nilpotent, while $B(n)$ consists of 2-groups. Suppose $B \cup S_3$ carries almost all subgroups (premises (a)–(c) of §9.8, unproved). Then the probability that a random subgroup is nilpotent tends to 1 along even $n$, and to 0 along odd $n$, like $6/y \asymp n^{-1/4}$. The paper's second question asks "is it possible that the probability that a random subgroup of Sn is nilpotent tends to 0"; on this heuristic the answer is no, because it fails along even $n$. That is consistent with its Theorem 4, which bounds the probability away from 1 for $n \equiv 3 \pmod 4$. The third question asks whether a random subgroup has no orbits longer than some absolute constant $C$. The family gives $C = 8$ only by definition, so this is premise (c) restated, not evidence for it. The support for (c) is dimension. By [RoTr25] Proposition 5.3(i), "If G is excessive then G has excess i ∈ {1,2} and n ≤ 32", so a transitive 2-group of degree $b \ge 16$ has $d \le b/4 + 2 \le 3b/8$, against $b/2$ for blocks of 8T22. An orbit of size $b \ge 16$ therefore lowers $D$ by at least $b/8$. Whether the extra orbit data compensate is not settled here.

### 9.8 Conjecture 9.3, and what it rests on

**Conjecture 9.3.** (i) $\log_2 f(n) = \frac{n^2}{16} + \frac78 n\log_2 n + O(n)$. (ii) For a uniformly random $H \le S_n$, $\log_2|H| = (\frac38 + o(1))\,n$ in probability.

It rests on three premises, none proved. (a) Subgroups that are not 2-groups, outside $S_3(n)$, are negligible. (b) 2-subgroups that are not FC are negligible; at $n = 8$ they are 3150 of 65 054, or 4.8%. By $n = 14$ the share is 16.5% (see the census below). (c) FC subgroups with an orbit longer than 8 are negligible. Under (a)–(c), part (i) follows from Proposition 9.4, and part (ii) follows if the concentration of §9.6 persists. For odd $n$ it is carried by $S_3$, whose members have order $3\cdot|H/A_3|$, and $\log_2|H/A_3|$ concentrates in the same way.

**Proposition 9.4.** $\log_2 B(n) = \frac{n^2}{16} + \frac78 n\log_2 n + O(n)$. Hence Conjecture 9.3(i) is equivalent to $f(n) \le 2^{O(n)}B(n)$.

*Proof.* $B(n) \ge F_8(n)$ gives the lower bound, by Theorem 9.1. For the upper bound, every transitive 2-group of degree $b \le 8$ has $d \le b/2$ (census), so $D \le n/2$, and each orbit datum carries at most $G_D \le 8\cdot2^{n^2/16}$ subgroups, by the bound in the proof of Proposition 9.2. The orbit data are counted by $T(n) = n!\,[x^n]\exp(x + x^2/2 + 7x^4/24 + 12945x^8/8!)$, and Cauchy's bound at $x = n^{1/8}$ gives $\log_2 T(n) \le \frac78 n\log_2 n + O(n)$. ∎

The evidence at small $n$ is weak. Shares of $f(n)$ (A005432), with $\alpha_f(n) = (\log_2 f(n) - n^2/16)/(n\log_2 n)$:

| $n$ | $e/f$ | $B/f$ | $\|S_3\|/f$ | outside $B \cup S_3$ | $\alpha_f$ |
| --- | --- | --- | --- | --- | --- |
| 6 | 0.1863 | 0.4337 | 0.1512 | 0.4151 | 0.5324 |
| 7 | 0.1165 | 0.3024 | 0.2261 | 0.4715 | 0.5293 |
| 8 | 0.0809 | 0.4094 | 0.1100 | 0.4807 | 0.5503 |
| 9 | 0.0479 | 0.2831 | 0.1547 | 0.5621 | 0.5479 |
| 10 | 0.0380 | 0.3674 | 0.0725 | 0.5601 | 0.5590 |
| 11 | 0.0237 | 0.2538 | 0.1655 | 0.5807 | 0.5526 |
| 12 | 0.0181 | 0.3618 | 0.0672 | 0.5710 | 0.5649 |
| 13 | 0.0115 | 0.2549 | 0.1515 | 0.5936 | 0.5569 |
| 14 | 0.0101 | 0.3779 | 0.0583 | 0.5638 | 0.5649 |
| 15 | 0.0062 | 0.2476 | 0.1694 | 0.5830 | 0.5575 |
| 16 | 0.0054 | 0.4121 | 0.0548 | 0.5331 | 0.5663 |
| 17 | 0.0034 | 0.2756 | 0.1759 | 0.5485 | 0.5576 |
| 18 | 0.0031 | 0.4544 | 0.0486 | 0.4970 | 0.5658 |

The share outside $B \cup S_3$ is between 41.5% and 59.4% for $6 \le n \le 18$, and 49.7% at $n = 18$. Along even $n$, $B/f$ rises from 0.362 at $n = 12$ to 0.454 at $n = 18$; along odd $n$ it stays between 0.248 and 0.302. The formula that (a)–(c) would give, $G_{n/2}B^*(n)$ for even $n$ and $G_{(n-1)/2}\big(nB^*(n-1) + \binom n3 B^*(n-3)\big)$ for odd $n$, is recorded but not conjectured. $f(n)$ divided by it is 1.21, 1.40, 1.82, 2.10, 2.24, 2.27 and 2.19 at $n = 6, 8, \dots, 18$, and 0.75, 0.69, 0.92, 1.22, 1.45, 1.74, 1.90 and 1.95 at $n = 3, 5, \dots, 17$; both sequences are still drifting. $\alpha_f(n)$ lies between 0.557 and 0.566 for $12 \le n \le 18$. With the linear term $cn$, $c = -1.96$, of §9.4, that does not conflict with $7/8$, and nothing at $n \le 18$ can separate $3/4$ from $7/8$.

**Added November 1, 2026: every 2-subgroup of $S_n$ for $n \le 15$.** Until now the evidence for (b) was the single value at $n = 8$. `src/two_subgroups.py` counts all 2-subgroups of $S_n$, FC or not, for $n \le 15$. Fixed points split off: if $t(k)$ is the number of fixed-point-free 2-subgroups of $S_k$, then $S_n$ has $\sum_k \binom nk t(k)$ of them. A fixed-point-free 2-subgroup whose orbits $O_1, \dots, O_t$ have sizes $b_i \in \{2, 4, 8\}$ lies in a conjugate of $R = \prod_i R_{b_i}$, a Sylow 2-subgroup of $\prod_i \mathrm{Sym}(O_i)$. So the subgroups with these orbits are counted by listing the subgroups of $R$ that are transitive on every $O_i$ and weighting each by $\prod_i n_2(b_i)/s_{b_i}(H|_{O_i})$, where $n_2(b) = b!/2^{b-1}$ is the number of Sylow 2-subgroups of $S_b$ and $s_b(P)$ the number of them that contain $P$. For transitive $P$, $s_b(P)$ is the number of maximal chains of subgroups from a point stabiliser $P_x$ up to $P$: a Sylow 2-subgroup of $S_b$ is the stabiliser of a nest of block systems with blocks of sizes $2, 4, \dots, b/2$, and the block systems of $P$ correspond to the subgroups between $P_x$ and $P$. The listing descends through maximal subgroups and stops at the first one that is intransitive on some $O_i$. The largest $R$, of order 2048, has 69 104 subgroups transitive on every orbit, for orbit type $(8, 4, 2)$; the longest list has 123 864 subgroups, for $(4, 4, 4, 2)$. The program passes its 122 checks. `src/verify_two_subgroups.py`, which shares no code with it, lists the same subgroups on explicit permutations, finds each $s$ by testing every Sylow 2-subgroup of each $\mathrm{Sym}(O_i)$, and for $n \le 12$ also conjugates the list over a transversal and counts the distinct subgroups, compared exactly as sorted lists of their elements, with no weights at all. The two agree at every $n \le 15$, by order and by FC status. The FC counts equal $B(n)$, computed from the formula of §9.1, at every $n \le 15$, the elementary abelian counts equal $a(n, k)$, and for $n \le 8$ every count equals `census8.py`. An OEIS search for 1, 2, 4, 20, 76, 631, 3417, 65054 (the counts for $n \le 8$) returned no sequence on November 1, 2026.

| $n$ | 2-subgroups | not FC, of 2-subgroups | not FC, of $f$ | 2-subgroups, of $f$ | $\|S_3\|/f$ | neither, of $f$ |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 65 054 | 0.0484 | 0.0208 | 0.4302 | 0.1100 | 0.4598 |
| 9 | 508 150 | 0.0558 | 0.0167 | 0.2998 | 0.1547 | 0.5454 |
| 10 | 11 666 461 | 0.0680 | 0.0268 | 0.3942 | 0.0725 | 0.5333 |
| 11 | 110 251 571 | 0.0698 | 0.0190 | 0.2728 | 0.1655 | 0.5617 |
| 12 | 4 498 662 367 | 0.1480 | 0.0628 | 0.4246 | 0.0672 | 0.5082 |
| 13 | 52 935 128 039 | 0.1562 | 0.0472 | 0.3021 | 0.1515 | 0.5464 |
| 14 | 2 557 694 647 600 | 0.1649 | 0.0746 | 0.4525 | 0.0583 | 0.4891 |
| 15 | 34 736 519 933 806 | 0.1656 | 0.0491 | 0.2968 | 0.1694 | 0.5339 |

The last column is the share of $f(n)$ that is neither a 2-group nor in $S_3(n)$, which is what premise (a) asks to be negligible. Read even and odd $n$ separately, because the shares alternate with parity rather than move together: 2-subgroups are between 39.4% and 45.3% of $f(n)$ at even $n$ and between 27.3% and 30.2% at odd $n$, while the non-FC share of 2-subgroups is slightly higher at each odd $n$ than at the even $n$ before it.

* (c) is not tested. The orbits of a 2-group have 2-power length, so the first orbit longer than 8 has 16 points, and for $n \le 15$ the share that (c) asks to be negligible is zero by necessity, not by evidence. An orbit of 16 points needs the Sylow 2-subgroup of $S_{16}$, of order $2^{15}$, and the census was not run that far.
* (b) runs the wrong way. The share of 2-subgroups that are not FC is 4.8%, 6.8%, 14.8% and 16.5% at $n = 8, 10, 12, 14$, and 5.6%, 7.0%, 15.6% and 16.6% at $n = 9, 11, 13, 15$. As a share of $f(n)$ it is 2.1%, 2.7%, 6.3% and 7.5% along even $n$. A 2-subgroup with at most one orbit of length at least 4 is FC, because squares are trivial on orbits of length 2: for that orbit $O$, $(h|_O)^2 = h^2 \in H$ for every $h \in H$. So every non-FC subgroup has two or more orbits of length 4 or 8, and the census confirms that the other orbit types contribute none. Those orbit types take over the count: the 2-subgroups with two or more orbits of length at least 4 are 24.7%, 43.2%, 71.8% and 86.9% of all 2-subgroups at $n = 8, 10, 12, 14$, and among them the share that is not FC is 19.6%, 15.8%, 20.6% and 19.0%. So the rise comes from these orbit types spreading. The non-FC share of all 2-subgroups is the product $q_n r_n$, where $q_n$ is the share with two or more orbits of length at least 4 and $r_n$ the non-FC share among them. It would tend to about a fifth, not to zero, only if both $q_n \to 1$ and $r_n$ stays near one fifth; up to $n = 14$ the census shows $q_n$ rising and $r_n$ roughly steady, and it proves neither. Up to $n = 15$ the share has risen at every step along both parities.
* (a) carries the largest share. The subgroups that are neither 2-groups nor in $S_3(n)$ are between 46.0% and 56.2% of $f(n)$ for $8 \le n \le 15$. Along even $n$ they are 46.0%, 53.3%, 50.8% and 48.9%, falling since $n = 10$, slowly and from about half. The column "outside $B \cup S_3$" of the table above is this share plus the non-FC 2-subgroups: at $n = 14$, 56.4% is 48.9% plus 7.5%.

So for $8 \le n \le 15$ the data support none of the three premises: (c) cannot be tested, the share that (b) needs to vanish grows, and the share that (a) needs to vanish stays near half. The 4.8% quoted for (b) is the smallest non-zero value in the range. This does not disprove Conjecture 9.3, which is asymptotic and remains a conjecture, but premise (b) now has small-$n$ evidence against it rather than for it.

### 9.9 The double root of §7, redone

Apply the substitution of §7 to the ambient group $\prod 2^{1+4}_+$ over $n/8$ blocks of 8, which has order $2^{5n/8}$ and a Frattini subgroup of order $2^{n/8}$. Its FC subgroups of order $2^{\nu n}$ number about $2^{n^2(\nu - 1/8)(5/8 - \nu)}$, and $(\nu - \tfrac18)(\tfrac58 - \nu) = \tfrac1{16}$ is $\nu^2 - \tfrac34\nu + \tfrac9{64} = 0$: discriminant exactly zero, double root $\nu = 3/8$. This is automatic. Suppose a family has $2^{n^2 g(\nu) + o(n^2)}$ members of order $2^{\nu n}$, with $g$ concave. Then $\max g \le 1/16$, because $f(n) = 2^{n^2/16 + o(n^2)}$. If the family reaches $2^{n^2/16}$, then $g(\nu) = 1/16$ holds only at the maximiser, which is a double root when $g$ is smooth there. So the double root of §7 says that the elementary abelian family reaches the leading term and concentrates at $n/4$, which §3 already shows. It says nothing about which order is typical among all subgroups.

### 9.10 Proposition 9.5: $e(n)$ is asymptotically $G_{\lfloor n/2\rfloor}$ times the number of 2-, 4-block data, and the rank law is a theorem (added October 19, 2026)

**Proved, not machine-checked.** This answers open questions (vi) and (iii) of §8 with explicit but very weak rates: the proofs lose a factor of about $2^{n/64}$ against the $2^{n/4}$ the data show, and their thresholds lie beyond the computed range, so nothing here is checkable against the computations of §9.2 and §9.3, whose residuals are far smaller than the bounds proved. Write $m = \lfloor n/2\rfloor$, let $M(n)$ be as in §9.2, and put $M'(n) = M(n)$ for even $n$ and $M'(n) = n\,M(n-1)$ for odd $n$.

**Proposition 9.5.** (a) $e(n) = G_m\,M'(n)\,(1 + \varepsilon_n)$ with $|\varepsilon_n| \le 24n\,2^{-(n-2)/64}$ for even $n \ge 1000$ and $|\varepsilon_n| \le 640n^2\,2^{-(n-3)/32}$ for odd $n \ge 1000$. (b) For $n \ge 1000$ and $k \ge n/8$, $a(n,k) = \binom{m}{k}_2 M'(n)\,(1 + \varepsilon_{n,k})$ with $|\varepsilon_{n,k}| \le 640 n^2 2^{-n/64}$. (c) For $n \ge 2000$ the total variation distance between the rank law $a(n,\cdot)/e(n)$ and the law $\binom{m}{\cdot}_2/G_m$ of the dimension of a uniformly random subspace of $\mathbb{F}_2^m$ is at most $2000n^2 2^{-n/64}$; in particular it tends to $0$. (d) For $n \equiv 0 \pmod 4$ and fixed $j$, $a(n, n/4+j)/a(n, n/4) \to 2^{-j^2}$, the law of §3.

*Setup.* An elementary abelian $E \le S_n$ acts regularly on each orbit, so its orbits have sizes $2^{j_i}$ and the group it induces on an orbit of size $2^j$ is a regular $\mathbb{F}_2^j$. Call the orbits together with these induced groups the *orbit data* of $E$. By the correspondence of §9.1 (every $\Phi(P_i)$ is trivial here), the elementary abelian subgroups with given orbit data are exactly the subspaces of $\prod_i P_i \cong \mathbb{F}_2^D$, $D = \sum_i j_i$, that project onto every factor. Three counting facts, used throughout:

1. *Hyperplanes.* A subspace $V \le \mathbb{F}_2^D$ fails to project onto the factor $\mathbb{F}_2^{j_i}$ exactly when its image lies in one of the $2^{j_i}-1$ hyperplanes of that factor, that is, when $V$ lies in one of $2^{j_i}-1$ hyperplanes of $\mathbb{F}_2^D$. A hyperplane contains $G_{D-1}$ subspaces and $\binom{D-1}{k}_2$ of dimension $k$, and $\binom{D-1}{k}_2 / \binom{D}{k}_2 = (2^{D-k}-1)/(2^D-1) \le 2^{-k}$, while $G_{D-1}/G_D \le 2^{-(D-1)/2}$ by the Jensen step in the proof of Theorem 9.1. Iterating, $G_{m-\delta}/G_m \le 2^{-\delta(2m-\delta-1)/4} \le 2^{-\delta(m-1)/4}$ and $\binom{m-\delta}{k}_2 \le 2^{-\delta k}\binom{m}{k}_2$ for $0 \le \delta \le m$; both are decreasing in $\delta$.
2. *Orbit data are few.* Order the $n$ points and cut the order into consecutive blocks whose sizes are powers of 2; identify the $i$-th point of a block of size $2^j$ with the binary expansion of $i$ in $\mathbb{F}_2^j$ and let the block carry the translations. Every regular $\mathbb{F}_2^j$ on a block arises from some ordering of the block, so every orbit datum arises from some ordering and cut. Hence there are at most $T(n) = n!\,2^{n-1} \le n^{2n}$ orbit data in all, and at most $t!\,2^t$ on a given set of $t$ points.
3. *Dimension loss.* A point in an orbit of size $2^j$ contributes $j/2^j$ to $D$: $1/2$ for $j = 1, 2$, at most $3/8$ for $j \ge 3$, and $0$ for a fixed point. With $n_0$ fixed points and $n_+$ points in orbits of size at least 8, $D \le n/2 - n_0/2 - n_+/8$.

Call the $s = n_0 + n_+$ points that are fixed or in orbits of size at least 8 *special*. Choosing the special points, which of them are fixed, the orbit data on the non-fixed special points (fact 2) and a partition of the rest into 2- and 4-blocks gives at most
$$\sum_{n_0+n_+=s}\binom{n}{s}\binom{s}{n_0}\,n_+!\,2^{n_+}\,M(n-s) \;\le\; (s+1)\,(4n)^s\,M(n-s)$$
orbit data with $s$ special points, since $\binom{n}{s} \le n^s/s!$, $\binom{s}{n_0} \le 2^s$, $n_+! \le s!$ and $2^{n_+} \le 2^s$. Here $M(n-s) \le M(n')$ whenever $n - s \le n'$ and both are even, by extending a partition with 2-blocks. Two sums recur: for $0 < q \le 1/2$, $\sum_{s\ge1}(s+1)q^s = q(2-q)/(1-q)^2 \le 6q$ and $\sum_{r\ge2}(r+2)q^r = q^2\bigl(q/(1-q)^2 + 4/(1-q)\bigr) \le 10q^2$.

*Proof of (a), even $n$.* The data with no special point have every orbit of size 2 or 4, so $D = m = n/2$, and they number $M(n)$: a partition of $[n]$ into 2-blocks and 4-blocks, each 2-block carrying its transposition and each 4-block the unique regular $V_4$ of $S_4$. By fact 1 such a datum carries between $G_m - (t_1 + 3t_2)G_{m-1}$ and $G_m$ subgroups, where $t_1, t_2$ are the numbers of 2- and 4-blocks; $t_1 + 3t_2 \le \tfrac32(t_1 + 2t_2) = 3n/4$, so it carries at least $G_m(1 - \tfrac34 n\,2^{-(m-1)/2})$. A datum with $s \ge 1$ special points has $D \le m - s/8$ by fact 3, hence carries at most $G_m 2^{-(s/8)(m-1)/4} = G_m 2^{-s(n-2)/64}$ subgroups by fact 1, and $n - s$ is even, so $M(n-s) \le M(n)$. With $q = 4n\,2^{-(n-2)/64}$, which is below $1/2$ for $n \ge 1000$,
$$M(n)G_m\bigl(1 - \tfrac34 n\,2^{-(n-2)/4}\bigr) \;\le\; e(n) \;\le\; M(n)G_m\Bigl(1 + \sum_{s\ge1}(s+1)q^s\Bigr) \;\le\; M(n)G_m(1 + 6q),$$
and both error terms are at most $24n\,2^{-(n-2)/64}$.

*Proof of (a), odd $n$.* Every orbit has even size or is a fixed point, so $n_0 \ge 1$ and $s$ is odd, and by fact 3, $D \le (n-1)/2 - (n_0-1)/2 - n_+/8 \le m - (s-1)/8$. The data with $s = 1$ have $D = m$ and number $n\,M(n-1) = M'(n)$, and each carries at least $G_m(1 - \tfrac34 n\,2^{-(n-3)/4})$ subgroups, as before. Data with $s \ge 3$ number at most $(s+1)(4n)^s M(n-1)$ and carry at most $G_m 2^{-(s-1)(n-3)/64}$ each. With $q = 4n\,2^{-(n-3)/64} \le 1/2$ their total is at most $M(n-1)G_m\sum_{s\ge3}(s+1)(4n)^s 2^{-(s-1)(n-3)/64} = M'(n)G_m\cdot 4\sum_{r\ge2}(r+2)q^r \le 40q^2\,M'(n)G_m$, and $40q^2 = 640n^2 2^{-(n-3)/32}$. ∎

*Proof of (b).* Repeat the argument for subspaces of dimension $k$. A datum without special points carries between $\binom{m}{k}_2(1 - \tfrac34 n\,2^{-k})$ and $\binom{m}{k}_2$ of them, by fact 1, and $\tfrac34 n 2^{-k} \le n2^{-n/8}$. For even $n$ a datum with $s \ge 1$ special points has $D \le m - \lceil s/8\rceil$, so carries at most $2^{-ks/8}\binom{m}{k}_2$ by fact 1; with $q' = 4n\,2^{-k/8} \le 4n\,2^{-n/64} \le 1/2$ the special data contribute at most $6q'\binom{m}{k}_2 M(n) \le 24n\,2^{-n/64}\binom{m}{k}_2M(n)$. For odd $n$, $D \le m - \lceil (s-1)/8\rceil$ and the special data with $s \ge 3$ contribute at most $40q'^2\binom{m}{k}_2 M'(n) \le 640n^2 2^{-n/32}\binom{m}{k}_2 M'(n)$. ∎

*Proof of (c).* Let $p_k = \binom{m}{k}_2/G_m$ and $r_k = a(n,k)/e(n)$, so $\sum_k p_k = \sum_k r_k = 1$. For $k \ge n/8$, (a) and (b) give $|r_k - p_k| \le p_k(|\varepsilon_{n,k}| + |\varepsilon_n|)/(1 - |\varepsilon_n|) \le 1500n^2 2^{-n/64}p_k$. For $k < n/8$, $m/2 - k > n/8 - 1/4$, so $(m/2-k)^2 \ge n^2/64 - n/8$; from $\binom{m}{k}_2 \le 2^{k(m-k)}/\eta$ and $G_m \ge \binom{m}{\lfloor m/2\rfloor}_2 \ge 2^{(m^2-1)/4}$, $p_k \le 2^{1/4-(m/2-k)^2}/\eta \le 2^{1/4 - n^2/64 + n/8}/\eta$; and since every datum has $D \le m$, fact 2 gives $a(n,k) \le T(n)\binom{m}{k}_2$, so $r_k \le n^{2n}p_k/(1-|\varepsilon_n|) \le 2n^{2n}p_k$. The tail $\sum_{k<n/8}(p_k + r_k)$ is therefore at most $(n/8)\cdot 3n^{2n}\cdot 2^{1/4-n^2/64+n/8}/\eta$, whose logarithm is below $2n\log_2 n + n/8 + 2\log_2 n - n^2/64 + 2 < -18000$ for $n = 2000$ and decreases from there. Adding the bulk and the tail gives the bound stated. ∎

*Proof of (d).* By (b) with $m = n/2$ even and $k = m/2 + j$, the ratio is $\binom{m}{m/2+j}_2/\binom{m}{m/2}_2$ times $1 + o(1)$. Writing $\binom{m}{k}_2 = 2^{k(m-k)}\prod_{i=1}^{k}\frac{1 - 2^{-(m-k+i)}}{1 - 2^{-i}}$, the numerator product lies between $1 - k2^{-(m-k)}$ and $1$, so for fixed $j$ both $\binom{m}{m/2 \pm j}_2$ and $\binom{m}{m/2}_2$ equal $2^{k(m-k)}(\eta^{-1} + o(1))$, and $k(m-k) = m^2/4 - j^2$. ∎

**What this does and does not say.** (a) proves the model of §9.2, and (c) and (d) prove the limit law of §3 and the convergence of §9.3, so open questions (iii) and (vi) are closed as stated. The rates are not the truth: the proofs count orbit data by $(4n)^s$, where the data show that one extra 8-orbit, weight about $n^2$ against a loss of $2^{-n/4}$, is the whole correction, and the thresholds lie beyond the computed range $n \le 512$. So the bounds proved here are not testable against `data/elemab_stats.json`; the theorem and the measurements are consistent, and that is all. The saddle-point step of §9.2, which gives the $-\sqrt6/8\,n^{-1/2}$ coefficient of §6, is still not justified. Nothing here bears on Conjecture 9.3, whose premises concern subgroups outside the elementary abelian family.

## 10. Reproduction

```sh
bash scripts/reproduce.sh              # everything, in order (~35 min, plus about
                                       # a quarter of an hour for the census of §9.8)
FAST=1 bash scripts/reproduce.sh       # ~1 min before the census; skips the n=8 lattice,
                                       # the n=10 brute force, everything above n=128,
                                       # and the census recount above n=12 (not re-timed)
python src/final_check.py              # check the artifacts and this file; exit 0 = pass
python src/final_check.py --coverage   # also list the numbers it does not compare
```

`final_check.py` runs 429 checks. The first part checks the artifacts in `data/` against OEIS, brute force and closed forms. The stored rows of `data/elemab_stats.json` are recomputed from fresh rank counts for $n \le 256$ and $n = 512$, 88 of the 119; recomputing the other 31 would take about three minutes, so they are checked against the §3 closed form for the sd and against $G_{n/2}M(n)$ from §9.2 for $\log_2 e(n)$, which catches an error larger than $10^{-7}$ in a stored sd or $10^{-9}$ in a stored $\log_2 e(n)$. The second part reruns `census8.py` and `block_families.py` into a temporary directory, requires their output to equal the stored JSON, and then compares the numbers printed in this file with fresh values. It finds each number by the words just before it, so a sentence edited without its check fails instead of passing unread.

The census of §9.8, added November 1, 2026, is checked the same way. `final_check.py` reruns `two_subgroups.py` to $n = 12$ and requires its output to equal the stored JSON there. The stored rows for $n = 13$, 14 and 15 are checked instead by identities that do not use the enumeration: each row must be the binomial transform of the stored fixed-point-free rows, each fixed-point-free total the sum over its orbit types, the number of subgroups of each order $2^j$ odd (Frobenius) and of the largest order $n!/2^{v_2(n!)}$, the FC counts those of the formula of §9.1, and the elementary abelian counts $a(n,k)$. The stored output of `verify_two_subgroups.py` must equal the census at every $n \le 15$; it is read, not rerun, because the recount takes a quarter to half an hour (timings below).

A clean pass of an unchallenged checker would prove nothing, and `final_check` distinguishes *skip* from *pass*. Under `FAST=1` the data stop at $n = 128$: every comparison that needs a larger $n$ is reported as skipped and counted in the summary line, and the check count above and the coverage check run only when nothing was skipped.

The first version of this section gave the count as 118 and listed five plants, one of them *"one stored standard deviation perturbed by $10^{-4}$"*, caught by *"the theta-law checks"*. Those checks read the stored rows for $n = 128$, 256 and 512 only. With the row for $n = 64$ perturbed instead, the checker as it stood on September 23, 2026 failed only its check count, which was stale. The stored-row recomputation was added for that reason.

Every plant below was run on September 23, 2026, in a copy of the repository, against the checker as committed; each made `final_check.py` exit with status 1. The first five are the five of the first version.

| planted defect | caught by |
| --- | --- |
| `data/elemab_stats.json`: the sd for $n = 64$ raised by $10^{-4}$ | the stored-row recomputation |
| `data/lattice_n8.json`: `total_subgroups` raised by 1 | the A005432 comparison and the §2.1 table |
| `src/elemab.py`: exponent $\binom{d}{2} \to \binom{d}{2}+1$ in the $q$-binomial inversion | 45 checks |
| `data/sd_precision.json`: `decimal_places` for $n = 512$ lowered from 35 to 12 | the high-precision sd floor |
| `data/sd_precision.json`: the last three of the 30 digits for $n = 512$, `694` → `999` | the 30-digit string comparison |
| §9.6: $f_8$ at $n = 512$ raised by one in its last digit | the comparison with a fresh value, and coverage |
| §9.4: $\alpha$ of $F_8$ at $n = 512$ raised by one in its last digit | the comparison with a fresh value, and coverage |
| §6: the Check 3 residual of the $n^{-1/2}$ basis raised by one in its last digit | the comparison with a fresh value, and coverage |
| `data/census8.json`: the involution count of 8T22 changed from 19 to 11 | the fresh rerun of `census8.py` |
| §9.9: $\tfrac9{64} \to \tfrac9{32}$ in the quadratic | the exact-string check |
| `data/block_families.json`: the sd of $\log_2\lvert H\rvert$ at $n = 512$ raised by $10^{-2}$ | the fresh rerun of `block_families.py` |
| a colon changed to a semicolon in a sentence that anchors a check | the anchor search, and coverage |
| a new sentence containing an uncompared decimal | coverage |
| `data/two_subgroups.json`: a count by order at $n = 14$ raised by 2, which keeps it odd | the binomial transform of the fixed-point-free rows, and the comparison with the recount |
| `data/two_subgroups.json`: the fixed-point-free total at $k = 14$ raised by 1 | the binomial transform, and the sum over orbit types |
| `data/two_subgroups.json`: an FC count by order at $n = 11$ raised by 1 | the fresh rerun to $n = 12$, the binomial transform, the formula of §9.1, and the comparison with the recount |
| `results/verify_two_subgroups.json`: the total at $n = 15$ lowered by 1 | the comparison with the recount |
| §9.8: the non-FC share at $n = 14$ in the census table raised by one in its last digit | the comparison with a fresh value, and coverage |

The last five rows were added on November 1, 2026 and run the same way, against the checker with the census checks in it. The two §9 programs were tested the same way, by a substitution in the source or by replacing one function at run time. All ten plants made the program exit with status 1. The seven plants of the census programs of §9.8, run on November 1, 2026 in a copy of the repository, did the same; `verify_two_subgroups.py` was run with `--n-max 10` for them. Its first version keyed the explicit conjugation on a 128-bit hash of each subgroup, so its agreement was a fingerprint comparison with a small chance of collision; it now keys on the exact sorted list of elements, and the transversal plant was rerun against that version.

| program | planted defect | failed |
| --- | --- | --- |
| `census8.py` | the class of order 128 dropped at degree 8 | 1 check: the class counts by order, against [GN] |
| `census8.py` | the Frattini rank $d$ capped at 3 | 1 check: 8T22 as the one class with $d = 4$ |
| `census8.py` | some intransitive rows counted as transitive | 4 checks: the class lists for degrees 2, 4 and 8, and 8T22 |
| `block_families.py` | $M(n)$ without the 4-blocks, once in the source and once by replacing the function | 3 gates each time |
| `block_families.py` | the second difference of $\log_2 M$ taken naively in floating point | 1 gate: the Richardson limit |
| `block_families.py` | the rank law of $n-1$ in place of that of $n$ | 2 gates: the total variation distance and its rate |
| `block_families.py` | $w_d$ with roots $2^{i+1}$ instead of $2^i$ | 25 gates |
| `block_families.py` | $L_{\mathrm{sq}}$ returning $\sum \dim V$ instead of $\sum (\dim V)^2$ | 7 gates |
| `block_families.py` | $f(12)$ lowered to $4\times10^9$ | 1 gate: $\lvert B(12)\rvert + \lvert S_3(12)\rvert \le f(12)$ |
| `two_subgroups.py` | the chains for $s_b(P)$ counted from the trivial subgroup instead of from $P_x$ | 43 checks, the first being $\sum_P n_2(b)/s_b(P)$ against `census8.py` |
| `two_subgroups.py` | the FC test skipping orbits of length 4 | 6 checks: the FC counts against `census8.py` and against the formula of §9.1 |
| `two_subgroups.py` | the factor $m_b!$ left out of the number of set partitions of an orbit type | 43 checks |
| `verify_two_subgroups.py` | the test $\Phi(H\vert_O) \le H$ skipped on orbits of length 4 | 3 checks: the comparison with `two_subgroups.json` at $n = 8$, 9 and 10 |
| `verify_two_subgroups.py` | one transversal element left out of the explicit conjugation | 13 checks |
| `verify_two_subgroups.py` | one Sylow subgroup of each block left out of the test for $s$ | the assertion that every weighted count is an integer |
| `verify_two_subgroups.py` | a lossy key in the explicit conjugation: the two smallest elements instead of the whole sorted list | 8 checks |

One plant that looked natural was not caught: the FC test of `two_subgroups.py` skipping orbits of length 8. It need not be. For $n \le 15$ an orbit type with an 8-orbit has at most one other orbit of length at least 4, and when exactly two orbits $O, O'$ have length at least 4, $h^2 = (h|_O)^2(h|_{O'})^2$ for every $h \in H$, so $(h|_O)^2 \in H$ if and only if $(h|_{O'})^2 \in H$; the test on the other orbit decides.

Timings on an idle 11-core arm64 machine, CPython 3.14.7. Ranges are over two runs; single figures are one run and should be read as indicative. The rows dated November 1, 2026 were measured on the same machine under heavy competing load, so they are not comparable with the others.

| step | time |
| --- | --- |
| `subgroups 8` | 883 s |
| `subgroups 1..7` | 8.8–11.6 s |
| `asymptotics.py --n-max 512` | 218 s |
| `verify_elemab.py --n-max 9` | 43.9–50.8 s |
| `verify_elemab.py --n-max 10` | 397 s and 913 s on two runs — it is very sensitive to competing load |
| `sd_precision.py` (to $n=512$, 60 dp) | 115 s |
| `census8.py` | 5.6–8.7 s |
| `block_families.py --n-max 512` | 4.4–6.1 s |
| `final_check.py` | 65.7–67.0 s over three runs |
| `two_subgroups.py` (to $n = 15$, November 1, 2026) | 43.7 s and 62.1 s, at a load average near 24 |
| `verify_two_subgroups.py --n-max 15 --dedup-max 12` (November 1, 2026) | 851 s at the same load with exact keys, 802 s of CPU and a peak of 446 MB; the explicit conjugation for orbit type $(8, 4)$ alone took 625 s. The two earlier runs with hashed keys took 1391 s and 1741 s |

$n=9$ was attempted (`subgroups 9`) and deliberately abandoned. After 11 minutes it had not finished the first class, whose cost is a full $9! = 362880$-element coset scan; scaling the $n=8$ run by the $9\times$ larger permutation group and the $1.87\times$ larger class count puts the whole run at four hours at best. It would also have re-verified a value already listed in A005432 rather than extending anything, and the enumerator is already validated against OEIS at all eight of $n=1,\ldots,8$. `final_check.py` carries the expected $n=9$ values (1694723 subgroups, 554 classes), so the check will run automatically if anyone completes it.

## References

* **[RoTr25]** C. M. Roney-Dougal and G. Tracey, *Subgroups of symmetric groups: enumeration and asymptotic properties*, arXiv:2503.05416v1 [math.GR], submitted March 7, 2025. <https://arxiv.org/abs/2503.05416> (v1 only as of the last fetch, no journal reference; 3 Semantic Scholar citations, none bearing on §7 — see §8.)
* **[Py93]** L. Pyber, *Enumerating finite groups of given order*, Ann. of Math. 137 (1993), 203–220. [RoTr25] reference 16.
* **[GN]** T. Dokchitser, *GroupNames*, "Transitive groups of degree up to 15", <https://people.maths.bris.ac.uk/~matyd/GroupNames/T15.html>, read September 22, 2026. Used for the orders of 8T1–8T35 and for the entry "Extraspecial group; = D4○D4", ES+(2,2), 32,49, at 8T22 (§9.1).
* **[Va99]** Various, *Some of Paul's favorite problems*, booklet for the conference "Paul Erdős and his mathematics", Budapest, July 1999 — the source of #1162 (5.73). Not consulted directly; not available to me.
* OEIS [A005432](https://oeis.org/A005432) (subgroups of $S_n$), [A000638](https://oeis.org/A000638) (conjugacy classes of subgroups), [A000085](https://oeis.org/A000085) (involutions), [A006116](https://oeis.org/A006116) (Galois numbers).
