# Erdős problem #200 — notes

**Question.** Let $L(N)$ be the length of the longest arithmetic progression of primes contained in $\{1,\dots,N\}$. Is $L(N) = o(\log N)$?

**Status: open.** From a local mirror of [erdosproblems.com/200](https://www.erdosproblems.com/200) (scraped 2026-09-22), the status is `OPEN`, the statement is *"Does the longest arithmetic progression of primes in $\{1,\ldots,N\}$ have length $o(\log N)$?"*, the recorded progress is *"It follows from the prime number theorem that such a progression has length $\leq(1+o(1))\log N$"*, the sources are `["ErGr79", "ErGr80"]`, and the linked sequence is `A005115`. The site itself is JS-rendered behind Cloudflare and has no JSON API, so the mirror — not a scrape — is the source used here. Re-checked 2026-09-22 against the upstream dataset [`teorth/erdosproblems`](https://github.com/teorth/erdosproblems) `data/problems.yaml`, which also records `state: "open"`; the live page carries no comments and no proof claims. The DeepMind Lean formalisation `FormalConjectures/ErdosProblems/200.lean` states the question and the known variant $L(N)\le(1+o(1))\log N$; both are `sorry`, so nothing is proved there. [`lean/`](lean/) and [`lean-pnt/`](lean-pnt/) discharge the second (§9) — unconditionally, with the literature constant $1$.

**This document does not solve the problem.** It records an explicit rigorous upper bound, an independent exhaustive re-derivation of $L(N)$ as far as compute allows, a validated heuristic predicting that the answer is *yes*, and a quantitative statement of why the standard method cannot reach it.

---

## 0. Summary of claims and their status

| # | Claim | Status |
|---|---|---|
| 1 | Structure theorem for APs of primes (Prop. 1) | proved below; checked numerically against all 26 records; **machine-checked in Lean** (`structure_of_primeAP`) |
| 2 | $\theta(L(N)) \le \log N + 1/(L(N)-1)$, hence $L(N)\le(1+o(1))\log N$ (Prop. 2) | proved below; not new — this is the standard PNT bound made explicit; **machine-checked in Lean** (`theta_lt_sharp`) |
| 3 | $L(N)$ exactly, for all $N \le 5.73\times10^{11}$ | **independently recomputed here** (exhaustive) |
| 4 | $L(N)$ exactly, for all $N < a(27)$ — certainly for all $N < a(26)=1.28\times10^{16}$, since $a(27)$ is known only to satisfy $a(27)\le 6.96\times10^{17}$ | cited (OEIS A005115, Luhn), verified for consistency, not recomputed |
| 5 | Hardy–Littlewood model of the AP count | validated against exact counts, geometric mean ratio 0.9939 over 30 tests |
| 6 | $L(N)\sim 2\log N/\log\log N$, so the answer is **yes** | heuristic only — *not* a proof |
| 7 | Selberg's sieve (per difference or joint), the large sieve and Gallagher's larger sieve cannot yield $L(N)<c\log N$ for any fixed $c<1$ | proved for the larger sieve at every $N$, and for the others as $N\to\infty$ at any fixed level $N^A$ (§6); computed at the 19 records, where the main term is at least $e^{1.41}$, and up to $N=10^{10^6}$, where the level needed exceeds $N^{2.33}$. **Corrected September 22, 2026**: the earlier "margin $\ge10^{3.2}$" measured a formula, not a sieve (§6.2) |
| 8 | Best known lower bound: $L(N)\gg\log_2^{(7)}N$ | cited (Green–Tao), verified verbatim from the source |
| 9 | $L(N)\le(1+o(1))\log N$, **unconditionally** | **machine-checked in Lean 4** (§9) — DeepMind's `erdos_200.variants.upper`, which they leave `sorry`. The PNT is imported from `PrimeNumberTheoremAnd`; `lean/` alone yields the weaker constant $1/\log2$ |

Nothing here improves the known upper bound. Claim 7 explains why that should be expected from sieve methods, and is the main negative content of these notes.

---

## 1. Structure of an AP of primes

Write $W_k=\prod_{p\le k}p$ for the primorial, so $\log W_k=\theta(k)$.

> **Proposition 1.** Let $a,a+d,\dots,a+(k-1)d$ be primes with $k\ge3$ and $d\ge1$. Then exactly one of the following holds.
>
> **(A)** $a>k$ and $W_k \mid d$.
> **(B)** $k$ is prime, $a=k$, and $(W_k/k)\mid d$.

*Proof.* First, a general observation. Let $p\le k$ be prime with $p\nmid d$. Then $d$ is invertible mod $p$, so as $i$ runs through $0,1,\dots,p-1$ — a subrange of $0,\dots,k-1$ since $p\le k$ — the terms $a+id$ run through every residue class mod $p$. Hence some term is divisible by $p$; being prime, that term *equals* $p$. As $a\le a+id$, this forces $a\le p\le k$.

Contrapositive: if $a>k$ then $p\mid d$ for every prime $p\le k$, and since these primes are distinct, $W_k\mid d$. That is (A).

Now suppose $a\le k$. Then $a$ is a prime $q\le k$. If $q\mid d$ then $a+d=q+d$ is a multiple of $q$ strictly larger than $q$, contradicting primality; so $q\nmid d$. The term of index $q$, namely $a+qd=q(1+d)$, is a multiple of $q$ strictly larger than $q$, so it cannot be a term of the progression: the index $q$ must exceed $k-1$, i.e. $q\ge k$. With $q\le k$ this gives $q=k$, so $a=k$ and $k$ is prime. Finally, for any prime $p\le k$ with $p\ne k$: if $p\nmid d$ then by the observation $a\le p<k=a$, absurd; so $p\mid d$, and therefore $(W_k/k)\mid d$. That is (B).

The two cases are mutually exclusive since $a>k$ in (A) and $a=k$ in (B). $\blacksquare$

Case B is not vacuous: $3,5,7$ and $5,11,17,23,29$ and $7,157,\dots,907$ are exactly the records for $k=3,5,7$. It is, however, rare — no record for $k\ge8$ is of Case B, and none of the 5,503,986 progressions enumerated in §3 for $k\ge8$ is either. Forgetting Case B is the natural bug here, and it changes the answer: see §3.2.

## 2. An explicit upper bound

> **Proposition 2.** If $\{1,\dots,N\}$ contains a $k$-term AP of primes, $k\ge3$, then
> $$N \;\ge\; \begin{cases}
> \min\bigl(p_k^{+}+(k-1)W_k,\;\; k+(k-1)W_k/k\bigr), & k \text{ prime},\\[2pt]
> p_k^{+}+(k-1)W_k, & k \text{ composite},
> \end{cases}$$
> where $p_k^{+}$ is the least prime $>k$. In either case $$W_k \;<\; \frac{k}{k-1}\,N, \qquad\text{equivalently}\qquad \theta(k) \;<\; \log N+\frac{1}{k-1}.$$

*Proof.* By Prop. 1 the progression is of Case A or Case B. In Case A, $d\ge W_k$ and $a\ge p_k^{+}$, so the last term is $a+(k-1)d\ge p_k^{+}+(k-1)W_k$. In Case B, $d\ge W_k/k$ and $a=k$, so the last term is $\ge k+(k-1)W_k/k>(k-1)W_k/k$. The Case B bound is the smaller of the two, giving $(k-1)W_k/k<N$; take logarithms and use $\log\frac{k}{k-1}=-\log(1-\tfrac1k)\le\frac{1}{k-1}$. $\blacksquare$

> **Corollary.** $\theta(L(N))<\log N+\tfrac{1}{L(N)-1}$, and since $\theta(x)\sim x$ (prime number theorem), $L(N)\le(1+o(1))\log N$. With an effective lower bound of the shape $\theta(x)>x(1-1/\log x)$ this becomes $$L(N)\;\le\;\log N+\bigl(1+o(1)\bigr)\frac{\log N}{\log\log N}.$$

This is the bound quoted on the problem page; Prop. 2 just makes the $o(1)$ explicit and identifies the binding constraint as Case B. The operational form used in `analysis.py` is sharper still: $L(N)\le\max\{k: \text{(Prop. 2 bound)}\le N\}$, computed directly from the primorial. That form is exact and needs no estimate for $\theta$.

*A note on the effective $\theta$ estimate.* The inequality $\theta(x)>x(1-1/\log x)$ for $x\ge41$ is attributed to Rosser–Schoenfeld (1962). I verified it numerically here for $41\le x<5\times10^6$ (348,500 evaluation points, one just below each prime where the inequality is tightest; zero violations, minimum margin 0.400), but I **could not retrieve the paper** — Project Euclid serves a JS wall and the AMS mirror returns 403 — so the equation number is not confirmed and is deliberately not cited. Nothing in §§1–3 depends on it; it enters only the asymptotic restatement in the Corollary.

**Sanity check.** At $N=a(26)=1.2783\times10^{16}$, Prop. 2 gives $L\le43$ while the truth is $L=26$. The bound's ratio to $\log N$ decreases to $1.0058$ at $N=10^{10^4}$, confirming that it really does realise the $(1+o(1))\log N$ ceiling and not something weaker.

## 3. Exhaustive computation

### 3.1 Method

Prop. 1 turns the search into something cache-friendly. For Case A every term lies in one residue class $r \bmod W_k$ with $\gcd(r,W_k)=1$; for each of the $\varphi(W_k)$ classes we sieve the arithmetic progression $r, r+W_k, r+2W_k,\dots$ into a bitset and look for $k$ positions in AP by intersecting shifted copies of the bitset. Case B is a separate, much smaller loop. `src/apsearch.c` implements this with pthreads.

### 3.2 Validation of the search

Prop. 1 is an *assumption* from the search's point of view, so it is checked against a program that does not use it. `src/ref_bruteforce.py` loops over every difference $d$ and every first term $a$ with no structural shortcut, shares no code with `apsearch.c`, and is written in numpy rather than C.

* **Agreement:** exact, on all 11 $(k,B)$ pairs tested — both the total number of progressions and the minimal last term.
* **Negative control:** recompiling with Case B deliberately disabled makes the two disagree, as it must — $k=3$: 186267 vs 186647 progressions, minimum 17 vs 7; $k=5$: 13789 vs 13835, minimum 127 vs 29; $k=7$: 2508 vs 2515, minimum 1307 vs 907. So the cross-check has teeth; it is not passing vacuously.

### 3.3 Results

Exhaustive enumeration reproduced $a(k)$ — the least possible last term of a $k$-AP of primes, OEIS [A005115](https://oeis.org/A005115) — together with its first term and common difference, for $k=3,\dots,20$.

The $k=20$ run went to $B=a(20)=572945039351$ and found **exactly one** 20-term progression, namely $214861583621+18846497670i$, $0\le i\le19$. A 21-term AP with last term $\le B$ would contain two distinct 20-term APs with last term $\le B$ (its terms $1..20$ and its terms $2..21$), so a single 20-term AP rules out any 21-term AP below $B$. The same argument at $k=19$, whose run to $B=1.5\times10^{11}$ also returned exactly one progression, independently rules out a 20-term AP below $1.5\times10^{11}$. Hence:

> $L(N)$ has been **independently re-derived here for every $N\le5.73\times10^{11}$**.

Beyond that range, $L(N)=\max\{k:a(k)\le N\}$ is determined by the tabulated records $a(21),\dots,a(26)$, which are cited rather than recomputed. $L(N)$ is unknown for $N\ge a(27)$; the least known 27-term AP ends at $696112717486210091$, and its minimality is not proved.

### 3.4 The record table

All 26 known values of $a(k)$ were checked four independent ways: the terms form an AP, `sympy.isprime`, an independently written deterministic Miller–Rabin, distinctness, and agreement with Norman Luhn's records page (a source independent of OEIS). 26 records, 0 failures. Minimality of $a(k)$ for $k=23,\dots,26$ is itself a computational result of Perrenet and Petukhov, cited from Luhn's page, not reproved here.

The relevant column is $k/\log a(k)$ — the largest value of $L(N)/\log N$ attained when $L$ jumps to $k$:

| $k$ | $a(k)$ | $k/\log a(k)$ |
|---|---|---|
| 3 | 7 | 1.5417 |
| 10 | 2089 | 1.3081 |
| 13 | 725663 | 0.9633 |
| 17 | 4827507229 | 0.7624 |
| 20 | 572945039351 | 0.7387 |
| 23 | 449924511422857 | 0.6817 |
| 26 | 12783396861134173 | 0.7011 |

The ratio has fallen from $1.54$ to about $0.70$ and is drifting down slowly and non-monotonically. §5 explains why this is exactly what the conjecture predicts and why it is not evidence either way.

## 4. The Hardy–Littlewood model

For a fixed difference $d$, the expected number of $k$-APs of primes with first term near $t$ is $\mathfrak S(d)\prod_{i=0}^{k-1}1/\log(t+id)$, where the singular series is
$$\mathfrak S(d)=\prod_p\frac{1-\omega_d(p)/p}{(1-1/p)^k},\qquad
\omega_d(p)=\begin{cases}1,&p\mid d\\ \min(k,p),&p\nmid d.\end{cases}$$
For $W_k\mid d$ this factors as $\mathfrak S(W_km)=A_kT_k\,g(m)$ with $A_k=\prod_{p\le k}(1-1/p)^{1-k}$, $T_k=\prod_{p>k}(1-k/p)(1-1/p)^{-k}$, and $g(m)=\prod_{p\mid m,\,p>k}\frac{1-1/p}{1-k/p}$. Summing over $d=W_km$ gives
$$C_k(N)\;\approx\;\frac{A_kT_k\,\mathbb E[g]\,N^2}{2(k-1)W_k(\log N)^k}.$$
The shape $N^2/(2(k-1)\log^kN)$ times a singular series is the same normalisation Tao and Teräväinen state in Example 1.7 of [arXiv:2107.02158](https://arxiv.org/abs/2107.02158), which is an independent check on the constant.

**Validation.** Predicted counts were compared with *exhaustively measured* counts at 30 $(k,N)$ points, $k=8,\dots,16$ and $N$ from $9.8\times10^5$ to $8\times10^9$, with measured counts ranging from 22 to 4,513,153:

> geometric mean predicted/observed = **0.9939**, geometric s.d. **1.036**, full range 0.900 – 1.093.

The model was fitted to nothing: every constant is the Hardy–Littlewood one. Reproduced by `src/compare_counts.py`, which reads the histograms emitted by the exhaustive runs and is re-checked by `src/final_check.py`.

**Calibration.** $a(k)$ is a *minimum*, so the right test is not whether $a(k)/\text{threshold}\approx1$ — it will not be — but whether $U_k=\exp(-C_k(a(k)))$ is uniform, as a Poisson model implies. Over the 21 records whose extremal AP is of Case A (the model covers Case A only; $k=3,5,7$ are excluded):

> Kolmogorov–Smirnov $D=0.128$, $p=0.86$; mean $U=0.494$ against an expected $0.500$.

A regression of $\log C_k(a(k))$ on $k$ gives slope $+0.098\pm0.042$ ($t=2.31$) over all 21 points, which looks like a mild systematic drift. It is not: restricting to $k\ge11$ gives $+0.007\pm0.058$ ($t=0.11$), $k\ge14$ gives $-0.090\pm0.058$, $k\ge17$ gives $-0.070\pm0.082$. The apparent trend is driven entirely by the small-$k$ records, where $a(k)$ is 23, 157, 1669 — far too small for an asymptotic prime-density model to apply.

**Cross-check against Granville.** Granville (*Prime Number Patterns*, Amer. Math. Monthly **115** (2008) 279–296, §2.1, eq. (2.1); [DOI:10.1080/00029890.2008.11920529](https://doi.org/10.1080/00029890.2008.11920529)) conjectures that the largest prime of the smallest $k$-AP is about $(e^{1-\gamma}k/2)^{k/2}$, and states that dividing $a(k)$ by it leaves a quotient "in $(2/5,2)$ for each $n$, $15\le n\le21$". My implementation reproduces that interval exactly for all seven values — 1.9915, 0.4025, 1.6709, 0.9773, 0.7749, 0.8356, 1.4021 — which confirms the formula is implemented as he wrote it. Extending to the records published since: 1.197, 2.206, 0.855, 0.569, **0.173** for $k=22,\dots,26$. The tight interval does not persist; $a(26)$ in particular is well below both his formula and the Hardy–Littlewood threshold. Under the Poisson model that is a $U=0.81$ event, unremarkable on its own, but it is an honest caveat against reading too much into the fit.

*(Eq. (2.1) and the tower below were read from a 400 dpi render of the page, not from the PDF text layer, which mangles both.)*

## 5. What the model predicts, and why no computation can confirm it

Setting $C_k(N)=1$ and using Mertens' $\prod_{p\le k}(1-1/p)^{-1}\sim e^{\gamma}\log k$:
$$2\log N \;\approx\; k\bigl(\log\log N-\log\log k+1-\gamma\bigr),$$
so
$$\boxed{\;L(N)\;\sim\;\frac{2\log N}{\log\log N}\;}$$
more precisely $L(N)\approx 2\log N/(\log\log N-\log\log\log N+1-\gamma)$. This is $o(\log N)$: **the conjectured answer to Erdős's question is yes.** It agrees with Granville's (2.1), which is the same statement solved for $a(k)$ — though the inversion is mine: Granville prints a conjecture about the *endpoint* $a(k)$, hedged with "about", not a ratio asymptotic for $L(N)$. No published source states $L(N)\sim2\log N/\log\log N$ in that form, so treat the constant 2 as re-derived here and corroborated by him, not quoted from him.

The decay is, however, only $2/\log\log N$, and that is the honest punchline:

| $N$ | conjectured $L(N)/\log N$ |
|---|---|
| $10^{16}$ | 0.678 |
| $10^{30}$ | 0.584 |
| $10^{100}$ | 0.456 |
| $10^{1000}$ | 0.316 |
| $10^{10^6}$ | 0.159 |
| $10^{10^9}$ | 0.105 |

The measured value at $N=a(26)\approx10^{16.1}$ is 0.701, against 0.678 predicted. The ratio is still above $2/3$ at $N=10^{30}$ and does not fall below $1/2$ until about $N=10^{100}$. **No computation will ever make the $o(1)$ visible.** Data sitting near 0.7 is not evidence against the conjecture; it is precisely what the conjecture predicts. Anyone arguing from the numerical table in either direction is over-reading it.

## 6. Why sieve methods give nothing

This is the substantive negative result. It is computed by `src/sieve_limits.py`, which writes `results/sieve_limits.json`, and re-derived by section G of `src/final_check.py`, which shares no code with it.

> **Result.** Selberg's $\Lambda^2$ upper-bound sieve, applied to each common difference separately or to all candidate progressions jointly, and the large sieve cannot prove $L(N)<c\log N$ for any fixed $c<1$ at any fixed level $N^A$ once $N$ is large: their main term alone tends to infinity (§6.1). Gallagher's larger sieve, and congruences modulo primes up to about $\log N$, give no more than the primorial bound of Prop. 2 (§6.3). For $c>1$ nothing is left to prove, since Prop. 2 already gives $L(N)<c\log N$ for large $N$. The shortfall is not a constant factor: at $k=c\log N$ a proof must save a factor $e^{(2/c-1)k}$ over the number of candidates, and these sieves save $e^{o(k)}$.

*Correction, September 22, 2026.* The previous version of this section argued from $B(k,N)=2^kk!\,C_k(N)$, the classical Selberg bound for a $k$-tuple summed over admissible $d$, and concluded from $B\ge10^{3.2}$ everywhere that sieve methods are vacuous. The conclusion survives, but that argument did not establish it; §6.2 says why. `src/barrier.py` still computes $B$, its docstring now says what $B$ is and is not, and `results/barrier.json` is byte-identical to the version first committed.

### 6.1 The main term

By Prop. 1 a $k$-AP of primes in Case A has $a>k$ and $d=mW_k$. Call the pairs $(a,m)$ with $a>k$, $m\ge1$ and $a+(k-1)mW_k\le N$ the *candidates*. For $d=mW_k$ the first term takes $X_d=N-k-(k-1)d$ values, so there are
$$T_0=\sum_{m=1}^{M}X_{mW_k}=M(N-k)-(k-1)W_k\,\frac{M(M+1)}{2},\qquad M=\Bigl\lfloor\frac{N-k-1}{(k-1)W_k}\Bigr\rfloor,$$
candidates, computed here in exact integers. $T_0=0$ once $(k-1)W_k\ge N-k$, which is Prop. 2's territory. A sieve proof of $L(N)<k$ has to show that fewer than one candidate has all $k$ terms prime.

**Per difference.** Fix $d=mW_k$ and sift the $X_d<N$ possible first terms. Write the local densities as $h(p)=g(p)/(1-g(p))$: $h(p)=1/(p-1)$ for $p\le k$, since every term is then congruent to $a$; $h(p)=k/(p-k)$ for $p>k$ with $p\nmid d$; and $h(p)=1/(p-1)\le k/(p-k)$ for $p>k$ with $p\mid d$. Selberg's bound with weights $\lambda_q$ supported on $q\le\xi$ is $X_dQ(\lambda)+R$, where $Q(\lambda)\ge1/G_d(\xi)$ with $G_d(\xi)=\sum_{q\le\xi}\mu^2(q)\prod_{p\mid q}h(p)$. The method bounds $R$ by a nonnegative quantity, so its output is at least $X_d/G_d(\xi)$. Replacing each $h(p)$ by $h^*(p)$, which is $1/(p-1)$ for $p\le k$ and $k/(p-k)$ above, gives a single $G^*(\xi)\ge G_d(\xi)$ for every $d$, so summed over $d$ the main term is at least $T_0/G^*(\xi)$. Progressions with a term that is itself a prime $\le\xi$ are removed by the sieve although they may be genuine; adding them back only raises the bound. The large sieve's bound $(X_d+\xi^2)/G_d(\xi)$ is at least the same $X_d/G_d(\xi)$. (With $\lambda=\mu$ and unlimited support the $\Lambda^2$ bound is exact, so all of this concerns main term plus remainder *estimate*. Such estimates are useful only up to the natural level, which for one $d$ is about $X_d<N$; $\xi=\sqrt N$ covers it.)

**Jointly.** Sifting all the pairs $(a,m)$ at once, rather than one $d$ at a time, gives $h(p)=1/(p-1)$ for $p\le k$ and, for $p>k$,
$$h(p)=\frac{1+(p-1)k}{(p-1)(p+1-k)}\;<\;\frac{k}{p-k},$$
since after cross-multiplying the two sides differ by $p(1-k)<0$. So $G^*$ dominates the joint sieve too, and its main term is again at least $T_0/G^*(\xi)$. Its natural level is below $T_0$: a class of pairs modulo $q$ holds $T_0/q^2+O(N/q+1)$ of them, so the error bound matches the main term once $q$ is near $T_0/N$. The records below nonetheless evaluate it at level $\max(T_0,N)$. A larger level can only lower this bound for the main term, so a value above 1 at a generous level settles every smaller one.

**Closed form.** Split $q=q_1q_2$ with $q_1\mid W_k$ and every prime factor of $q_2$ above $k$. The sum over $q_1$ is $\prod_{p\le k}\bigl(1+\frac{1}{p-1}\bigr)=W_k/\varphi(W_k)$. A squarefree $q_2\le\xi$ has at most $J=\lfloor\log\xi/\log p_k^+\rfloor$ prime factors, and the $j$-th elementary symmetric function of the $h^*(p)$, $k<p\le\xi$, is at most $S^j/j!$. Hence
$$G^*(\xi)\;\le\;\frac{W_k}{\varphi(W_k)}\,E_J(S),\qquad E_J(S)=\sum_{j=0}^{J}\frac{S^j}{j!},\qquad S=\sum_{k<p\le\xi}\frac{k}{p-k},$$
and the main term is at least $T/E_J(S)$, where $T=T_0\,\varphi(W_k)/W_k$. This needs only $S$, which is summed exactly below $2\times10^8$ and bounded beyond it by Dusart's inequality $\sum_{p\le x}1/p\le\log\log x+0.2614972\ldots+\frac{1}{10\log^2x}+\frac{4}{15\log^3x}$ for $x\ge10372$, the constant being Meissel–Mertens' (Theorem 6.10 of [arXiv:1002.0442](https://arxiv.org/abs/1002.0442), read from the rendered page; published as [DOI:10.1007/s11139-016-9839-4](https://doi.org/10.1007/s11139-016-9839-4), whose theorem numbering I have not checked). `sieve_limits.py` tests the inequality at all 11,077,665 primes in $[10372,2\times10^8]$ and finds no violation; that covers every real $x$ in the range, since the left side only jumps at primes.

**Asymptotics.** $S=O(k\log k)$: the primes in $(k,2k]$ contribute at most $kH_k$, because the $p-k$ are distinct integers in $[1,k]$, and each $p>2k$ contributes less than $2k/p$. Put $k=c\log N$ and $\xi=N^{A/2}$. Then $J\sim(A/2c)\,k/\log k$ and, since $S^j/j!\le(eS/j)^j$,
$$\log E_J(S)\;\le\;\log(J+1)+J\log\frac{eS}{J}\;=\;O\Bigl(\frac{k\log\log k}{\log k}\Bigr)\;=\;o(k).$$
Meanwhile $\log T=(2/c-1)k-o(k)$ provided $M\to\infty$, which holds for $c<1$. (For $c>1$, $W_k>N$ eventually and $T_0=0$.) So for fixed $c<1$ and fixed $A$ the main term is at least $e^{(2/c-1)k-o(k)}$, which tends to infinity: none of these sieves proves $L(N)<c\log N$ for large $N$.

**At the records.** At $N=a(k)-1$ for $k=8,\dots,26$ no $k$-AP of primes exists, so a sieve proof would have something true to prove. $G^*(\sqrt N)$ is computed exactly, by a depth-first sum over squarefree $q\le\sqrt N$ that agrees with brute force at $\xi=20000$ for $k=5,12,26$. Logarithms are natural except in the last column, which compares the old $B$.

| $k$ | $M$ | $\log T_0$ | $\log G^*(\sqrt N)$ | $\log(T_0/G^*)$ | $\log(T/E_J)$ at level $N$ | $\log(T/E_J)$ at level $\max(T_0,N)$ | $\log_{10}(B/T_0)$ |
|---|---|---|---|---|---|---|---|
| 8 | 1 | 5.247 | 2.978 | 2.269 | 1.640 | 1.640 | 4.80 |
| 9 | 1 | 5.242 | 3.317 | 1.924 | 1.269 | 1.269 | 5.61 |
| 10 | 1 | 5.236 | 3.817 | 1.419 | 0.783 | 0.783 | 6.37 |
| 11 | 10 | 14.014 | 4.826 | 9.188 | 7.296 | 7.160 | 4.99 |
| 12 | 10 | 14.024 | 5.521 | 8.503 | 6.550 | 6.449 | 5.66 |
| 13 | 2 | 12.822 | 4.760 | 8.062 | 6.028 | 6.028 | 7.16 |
| 14 | 94 | 21.266 | 6.186 | 15.080 | 11.694 | 11.384 | 6.15 |
| 15 | 412 | 24.298 | 7.059 | 17.239 | 14.020 | 11.497 | 6.30 |
| 16 | 441 | 24.502 | 7.904 | 16.598 | 13.264 | 10.514 | 6.86 |
| 17 | 591 | 27.985 | 7.853 | 20.132 | 17.285 | 14.673 | 7.42 |
| 18 | 1960 | 30.444 | 8.983 | 21.461 | 16.303 | 13.564 | 7.64 |
| 19 | 478 | 30.624 | 8.016 | 22.608 | 17.708 | 17.312 | 8.89 |
| 20 | 3108 | 34.423 | 8.820 | 25.603 | 20.923 | 18.257 | 9.02 |
| 21 | 32316 | 39.157 | 9.870 | 29.287 | 24.943 | 19.849 | 8.97 |
| 22 | 175472 | 42.589 | 11.127 | 31.463 | 27.285 | 21.728 | 9.06 |
| 23 | 91670 | 44.473 | 10.175 | 34.298 | 28.057 | 25.229 | 10.13 |
| 24 | 237293 | 46.420 | 10.749 | 35.671 | 29.561 | 26.642 | 10.57 |
| 25 | 1078257 | 49.490 | 11.488 | 38.001 | 32.122 | 26.849 | 10.82 |
| 26 | 2292031 | 51.039 | 12.143 | 38.896 | 33.140 | 27.686 | 11.27 |

All three main-term columns are positive in every row: no sieve of these kinds proves even the known fact $L(a(k)-1)<k$. The margin is smallest at $k=8,9,10$, where $M=1$ and there are only 188 to 190 candidates, and least at $k=10$: a main term of at least $e^{1.419}\approx4.1$ with the exact $G^*$, and $e^{0.78}\approx2.2$ from the closed form. From $k=11$ it grows, not monotonically, to $e^{38.9}$ at $k=26$. The joint level $\max(T_0,N)$ reaches $N^{1.376}$ at $k=26$; at $k=8,9,10,13$, $T_0<N$ and the two closed-form columns coincide.

**On a grid.** For $k=\lfloor c\log N\rfloor$ with $c\in\{0.3,0.5,0.7,0.9,1.0\}$ and $\log_{10}N\in\{10,16,30,100,10^3,10^4,10^5,10^6\}$, the closed form leaves a main term above 1 at level $N$ and at level $N^2$ in all 80 cases, by at least $e^{3.05}$ (at $\log_{10}N=10$, $c=1$, level $N^2$). The table gives the exponent $A$ of the level $N^A$ at which the closed form first allows a main term below 1. Below that level no Selberg or large-sieve proof exists; above it, the closed form no longer rules one out.

| $c$ | $\log_{10}N=10$ | $16$ | $30$ | $100$ | $10^3$ | $10^4$ | $10^5$ | $10^6$ |
|---|---|---|---|---|---|---|---|---|
| 0.3 | 39.3 | 18.3 | 17.6 | 9.11 | 9.77 | 11.7 | 13.4 | 14.7 |
| 0.5 | 4.90 | 4.48 | 5.33 | 5.92 | 6.83 | 8.20 | 9.36 | 10.6 |
| 0.7 | 3.20 | 3.84 | 3.91 | 4.29 | 5.15 | 6.25 | 7.23 | 8.17 |
| 0.9 | 2.72 | 2.94 | 3.04 | 3.39 | 4.09 | 4.83 | 5.56 | 6.22 |
| 1.0 | 2.34 | 2.55 | 2.68 | 2.98 | 3.59 | 4.14 | 4.87 | 5.53 |

Every entry exceeds 2, the least being $A=2.3398$, whereas the per-difference sieve's natural level is below $N$ and the joint sieve's is below $T_0<N^2$. From $\log_{10}N=100$ on, every row grows. To leading order $A$ grows like $(2-c)\log k/\log\log k$, but the lower-order terms are large: at $\log_{10}N=10^6$ that expression gives 5.46 against the computed 5.53 at $c=1$, and 8.80 against 14.7 at $c=0.3$.

Per term of the progression, a proof must save $\log T/k\approx2/c-1$ nats (5.668, 3.001, 1.858, 1.223 and 1.001 at $\log_{10}N=10^6$ for the five values of $c$). The sieve supplies $\log E_J/k$: between 0.23 and 1.62 at level $N$, and between 0.42 and 2.79 at level $N^2$, over the whole grid, falling as $N$ grows from $\log_{10}N=16$ on. For scale, at $\log_{10}N=16$ and $c=0.7$ ($k=25$) the asymptotic Hardy–Littlewood count that `barrier.py` uses predicts a saving of 1.959 per term against the 1.951 needed, so it expects $e^{25(1.951-1.959)}\approx0.8$ such progressions below $10^{16}$. That asymptotic form runs below the exact heuristic of §4 by a factor of about 1.9 (geometric mean 0.539 over 12 test points, `barrier.py` section 1), and one such progression does exist, since $a(25)=5.77\times10^{15}$. The sieve supplies 0.62 per term there.

Up to $\log_{10}N=10^4$, which is 30 of the 40 $(c,N)$ points, $T_0$ is an exact integer. Beyond that $\log T_0$ is taken as $2\log N-\log(2(k-1))-\theta(k)$, which agrees with the exact value to within $5\times10^{-11}$, that is to rounding error, at the 14 points where both are available and $M>e^{40}$. The grid is computed in floating point, and beyond $\xi=2\times10^8$ it depends on Dusart's inequality.

### 6.2 What was wrong with the previous argument

The previous version computed $B(k,N)=2^kk!\,C_k(N)$ and read $B\ge1$ as "the sieve is vacuous". The conclusion survives only because §6.1 proves it differently; the argument itself had six faults.

1. **The tuple bound is for fixed $k$.** $\bigl(2^kk!+o(1)\bigr)\mathfrak S(H)\,x/(\log x)^k$ holds as $x\to\infty$ with $k$ fixed, and the $o(1)$ depends on $k$. At $k=c\log N$ it is not justified as written, so even $B<1$ would not have been a proof.
2. **$B$ is not a possible sieve output.** Every Selberg main term is at most $T_0$, since $G\ge1$. Of the 106 points in `results/barrier.json`, the 82 with $\log_{10}N\le10^6$ can be compared; beyond that, $k$ exceeds the prime table's $2\times10^8$. $B>T_0$ at 69 of them, including all 57 with $c\ge0.7$ or $\log_{10}N\ge300$, and at all 19 records, by $10^{4.80}$ to $10^{11.27}$. There $B$ exceeds anything a sieve main term can be, so $B\ge1$ says nothing about the sieve. At the other 13 points, all with $c\le0.5$ and $\log_{10}N\le100$, the comparison is inconclusive. Over the 82, $\log_{10}(B/T_0)$ runs from $-2.67$ to $2.92\times10^6$.
3. **The *a fortiori* step was backwards.** Granting the tuple bound for every $d$ inflates $B$, which makes $B>1$ easier to reach and therefore less informative, not more. Showing that a sieve is vacuous needs a *lower* bound for its output, which is what $T/E_J$ is.
4. **The margins measured the formula.** $\min\log_{10}B=3.21$, and $\log_{10}B=2.8\times10^9$ at $\log_{10}N=10^9$, describe $B$, not any sieve. The margin that matters is the main term's: at level $N$ the lower bound of §6.1 exceeds 1 at all 82 comparable points, least $e^{2.40}$ at $\log_{10}N=2$, $c=0.9$. The $(\log k)^k$ growth of the old "overshoot" is $B$'s own excess over $T_0$: at the 6 comparable points with $k\ge10^5$, $\log(B/T_0)$ matches $\log\bigl[(2e^{\gamma-1}c\log k)^k\,e^{-k/(2\log k)}\sqrt{2\pi k}/(e^\gamma\log k)\bigr]$ to within 0.021 %, where $2e^{\gamma-1}\approx1.31$.
5. **"It does not even recover the primorial bound" was an artefact.** $T_0=0$ once $(k-1)W_k\ge N-k$, which is Prop. 2's Case A condition with $k+1$ in place of $p_k^+$, so any sieve over the candidates recovers the primorial bound automatically. Only the smooth formula for $C_k$ misses it.
6. **The diagnoses built on $B$ do not survive.** The claims that Elliott–Halberstam saves at most $2^k$ and that the gap is the parity problem were both read off $B$. Elliott–Halberstam does not enter, because the objects sifted are integers, not primes, and their counts in residue classes are elementary; there is no level of distribution for it to raise. The obstruction §6.1 measures is dimension against level: in dimension $k\asymp\log N$ a sieve at level $N^A$ sees at most $J\approx A\log N/(2\log k)$ prime factors above $k$ and saves $e^{o(k)}$, where $e^{(2/c-1)k}$ is needed.

### 6.3 The larger sieve, and congruences at small primes

Gallagher's larger sieve (*A larger sieve*, Acta Arith. **18** (1971) 77–81, [DOI:10.4064/aa-18-1-77-81](https://doi.org/10.4064/aa-18-1-77-81)) bounds a set of integers in $[1,N]$ that occupies at most $\nu(p)$ residue classes modulo each prime $p\le z$ by
$$\frac{\sum_{p\le z}\log p-\log N}{\sum_{p\le z}\log p/\nu(p)-\log N},$$
when the denominator is positive. Apply it to the $k$ terms of a Case A progression. Then $\nu(p)=1$ for $p\le k$, since $p\mid d$, and $\nu(p)\le k$ for $k<p\le z$; a proof must cover $d=W_k$, where $\nu(p)=k$ for all of them. With those values the numerator minus $k$ times the denominator is $(k-1)(\log N-\theta(k))$, so the bound is below $k$ exactly when $\theta(k)>\log N$, whatever $z$ is. That is the primorial condition $W_k>N$. Prop. 2 needs only $p_k^++(k-1)W_k>N$, so the larger sieve gives slightly less than Prop. 2. `sieve_limits.py` confirms the equivalence in all 15,210 cases with $3\le k\le80$, $N=10^2,10^3,\dots,10^{40}$ and $z\in\{k,2k,k^2,10^4,10^7\}$.

More generally, congruences modulo primes up to about $\log N$ cannot exclude a progression of any length $k\le\log N$. Put $W_z=\prod_{p\le z}p$. For every $k$, the progression $p_z^++iW_z$, $0\le i<k$, has no term divisible by a prime $\le z$, since each term is congruent to $p_z^+\not\equiv0$ modulo every such prime, and it lies in $[1,N]$ whenever $p_z^++(k-1)W_z\le N$. For $k\le\log N$ that allows every $z\ge3$ with $\theta(z)\le\log N-\log\log N$. So whatever excludes long progressions of primes must use primes above about $\log N$, and §6.1 bounds what a sieve can extract from those.

### 6.4 What is and is not covered

**Covered:** Selberg's $\Lambda^2$ sieve, per difference and jointly, with any nonnegative remainder estimate, at every level up to $N$ at the records (for the joint sieve, up to $\max(T_0,N)$) and up to at least $N^{2.33}$ on the grid; the large sieve over the same range of levels; Gallagher's larger sieve with any $z$; and congruences modulo primes up to about $\log N$. The asymptotic statement is for fixed $c<1$ and fixed $A$; it is not uniform as $c\to1$ or $A\to\infty$.

**Not covered:** combinatorial and $\beta$-sieves in dimension $\asymp\log N$, whose main terms are not of Selberg's form; sieves that require the first term to be prime and sift the other $k-1$ terms, where Bombieri–Vinogradov or Elliott–Halberstam would enter; weighted sieves; and every method that is not a sieve. None of this bears on whether $L(N)=o(\log N)$ is *true*, which §5 predicts; it concerns only which methods cannot show it.

## 7. Lower bounds: what is actually known

For completeness, in the other direction. The only *explicit* threshold in the literature is the Green–Tao expository note *A bound for progressions of length k in the primes*, which extracts from their theorem that the first $k$-AP of primes is bounded by an exponential tower. Verbatim from the abstract:

> "the argument shows that those primes are bounded by an exponential tower in $k$ of height seven"

and the final display is $O(NW)\le 2\uparrow2\uparrow2\uparrow2\uparrow2\uparrow2\uparrow2\uparrow O(k)$ — seven 2's. (Check it with `pdftotext refs/green_tao_envelope.pdf -`: the last line containing $\uparrow$ is `O(NW ) ⩽ 2 ↑ 2 ↑ 2 ↑ 2 ↑ 2 ↑ 2 ↑ 2 ↑ O(k).`, seven arrows. Do not count arrows document-wide; the note uses $\uparrow$ for ordinary exponentiation throughout, 185 times.) Green and Tao add the remark that "$100k$ will certainly suffice" for the implied constant. Inverting gives
$$L(N)\;\gg\;\log_2^{(7)}N,$$
a seven-fold iterated logarithm. That is the best explicit lower bound I could source.

Two corrections to things commonly repeated:

* **Granville renders the tower with eight 2's**, $2^{2^{2^{2^{2^{2^{2^{2^{100k}}}}}}}}$, one level taller than Green–Tao's own stated height of seven. Both readings were confirmed from the rendered pages. Granville's "8 twos" is what he wrote; it is not what Green–Tao claim.
* **Leng–Sah–Sawhney ([2402.17994](https://arxiv.org/abs/2402.17994), [2402.17995](https://arxiv.org/abs/2402.17995)) do not state any bound for the primes.** They improve the inverse theorem to quasipolynomial and Szemerédi to $r_k(N)\ll N\exp(-(\log\log N)^{c_k})$, both for *dense sets*. The word "prime" occurs in them only as a technical modulus and in historical remarks.

The real modern advance is Tao–Teräväinen ([2107.02158](https://arxiv.org/abs/2107.02158), JEMS 2023), which makes the $k$-AP-in-primes *count* effective with error $O((\log\log N)^{-c})$. In principle that yields a vastly better $N_0(k)$ than the tower, but they do not write one down and the $k$-dependence of the constants is left implicit, so no explicit $L(N)\gg(\log\log N)^{c}$ statement exists in the literature that I could find. I am not asserting one.

Teräväinen–Wang, *On the Green–Tao theorem for sparse sets* ([2603.09281](https://arxiv.org/abs/2603.09281), 10 Mar 2026), is the most recent quantitative work in this direction: if $\mathcal{A}$ has relative density $\delta$ in the primes up to $N$ and contains no nontrivial $k$-AP with $k\ge4$, then $\delta\ll\exp(-(\log\log\log N)^{c_k})$. Taken at $\delta=1$ — $\mathcal{A}$ the whole set of primes — this forces a $k$-AP once $(\log\log\log N)^{c_k}$ exceeds the implied constant, which is a *triple*-logarithmic threshold rather than a seven-fold tower. **It still does not produce an explicit lower bound for $L(N)$**, because inverting it requires knowing how $c_k$ and the implied constant degrade in $k$, and the abstract states neither. Recorded here as the current frontier, not as a bound. *I read the abstract only; I have not worked through the paper.*

**Upper bound.** No constant $c<1$ with $L(N)\le(c+o(1))\log N$ is known. That is the whole problem. A literature search (2026-09-22) found no improvement on the constant $1$ — not unconditionally, not under GRH, not under Elliott–Halberstam, and not under uniform Hardy–Littlewood. That is a negative search result, not a proof that none exists.

## 8. What would count as progress

1. Any $c<1$ with $L(N)\le c\log N$ for large $N$. By §6 this cannot come from Selberg's sieve, the large sieve or the larger sieve, and it is worth being precise about *why*, because the obvious escape route is already closed. One might hope to exploit the **rigidity** of an AP — all $k$ terms lie in one progression, rather than forming an arbitrary admissible $k$-tuple. But §6.1 already exploits that fully: by Prop. 1 only $d=mW_k$ are sifted, and the joint sieve over $(a,m)$ uses the shape of the progression at every prime. Rigidity is therefore already spent, and what remains is quantitative. At $k=c\log N$ a proof must save a factor $e^{(2/c-1)k}$ over the candidate count, while a sieve at level $N^A$ saves $e^{o(k)}$; the closed form of §6.1 does not allow a main term below 1 until $A$ reaches 2.3398 even in the most favourable case tested, and the required $A$ grows like $(2-c)\log k/\log\log k$, far beyond any level at which remainders are controlled. I know of no method that supplies the missing $e^{(2/c-1)k}$, and the neighbouring problems with patterns of that length (§8.1) are all open.
2. An unconditional bound $L(N)=o(\log N)$ with any decay rate, however slow — the conjecture says the truth is $2/\log\log N$, so even $\log N/\log\log\log N$ would be a first.
3. In the other direction, extracting an explicit $N_0(k)$ from Tao–Teräväinen would replace the seven-fold iterated logarithm with something enormously better, and is apparently just bookkeeping through Manners' constants — but it is bookkeeping nobody has published.

An accessible-*sounding* computational contribution would be $a(27)$: the least 27-term AP of primes is known only as an upper bound, $\le 696112717486210091$. But the cost model above rules it out here. For $k=27$, $W_{27}=223092870$ and $\varphi(W_{27})=36495360$, so $M=B/(26W_{27})\approx1.20\times10^{8}$ and the enumerator needs $\varphi(W)M^2/(2\cdot26\cdot64)\approx1.6\times10^{20}$ word-operations — about $3.2\times10^{12}$ core-seconds, or $10^{5}$ core-years at the $4.9\times10^{7}$/core-second measured in §3. That is four to five orders of magnitude beyond a single machine *for this enumerator*.

It does not follow that $a(27)$ is out of reach, and an earlier draft of these notes wrongly said the $a(23)$–$a(26)$ minimality proofs were distributed efforts. They were not. Joris Perrenet's [`arithmetic-progression`](https://github.com/jorisperrenet/arithmetic-progression/blob/6507fcb038a91a3711e486365c60b3aa37662778/README.md) records $a(25)$ "proven minimal-end" in **~2 days** and $a(26)$ in **~3 days**, both on a **single RTX 4090**, searching the restricted shape $(a+bn)\cdot23\#+c$ rather than all integers. The relevant gap is therefore not core-count but algorithm and hardware: a generic CPU enumerator over all residues is the wrong instrument by some five orders of magnitude. Perrenet attaches his own disclaimer — "Sadly I lost track of which ones were used to produce the records above" — so the runtimes are reported, not reproduced here. Verified independently: his AP-26 witness $(15626261+1666981n)\cdot23\#+59138353$ has last term $12783396861134173$, exactly OEIS $a(26)$.

### 8.1 Neighbouring problems

The local tracker in `../tracker-and-solver`, a separate project that is not part of this repository, was read on September 22, 2026 for problems close to #200: ones that ask for many values of a pattern to be prime at once, and ones about sifting or primorials at the scale $\log n$. `results/sieve_limits.json` records the tracker build that was read. The selection is mine and is not exhaustive. Statuses are the site's, as mirrored by the tracker, and agree with [`teorth/erdosproblems`](https://github.com/teorth/erdosproblems); the pattern column is my annotation of how many values must be prime at once.

| # | status | pattern | question | how it was settled, or what is known |
|---|---|---|---|---|
| [219](https://www.erdosproblems.com/219) | PROVED (LEAN) | $k$ fixed | arbitrarily long APs of primes? | yes: Green–Tao |
| [1187](https://www.erdosproblems.com/1187) | SOLVED | $k$ fixed | monochromatic prime $k$-APs in every finite colouring? | yes, from Green–Tao; its second question, with prime common difference, is trivially no |
| [457](https://www.erdosproblems.com/457) | PROVED (LEAN) | — | for some $\epsilon>0$, do all primes $p\le(2+\epsilon)\log n$ divide $\prod_{i\le\log n}(n+i)$ for infinitely many $n$? | yes: a CRT construction by GPT-5.2 Pro, prompted by Barreto, with an elaboration by Tao |
| [1202](https://www.erdosproblems.com/1202) | SOLVED | — | is there a $k$ such that removing any $(p_i-1)/2$ classes modulo each of any $k$ primes $p_i<n^{1-\epsilon}$ leaves at most $\epsilon n$ integers $m\le n$? | no: a construction by Price and GPT-5.4 Pro; the site notes that the large sieve gives yes when $p_k<n^{1/2}$ |
| [1140](https://www.erdosproblems.com/1140) | DISPROVED | $\sim\sqrt{n/2}$ | infinitely many $n$ with $n-2x^2$ prime for all $2x^2<n$? | no: Epure–Gica, with Mollin–Williams; all such $n$ are known, up to at most one exception |
| [1141](https://www.erdosproblems.com/1141) | DISPROVED (LEAN) | $\sim\sqrt n$ | infinitely many $n$ with $n-k^2$ prime for all $k^2<n$ coprime to $n$? | no: an internal OpenAI model (the site cites [APSSV26b]), deduced from a result of Pollack |
| [1142](https://www.erdosproblems.com/1142) | OPEN | $\log_2n$ | infinitely many $n$, or any $n>105$, with $n-2^k$ prime for all $1<2^k<n$? | none up to $2^{44}$ (Mientka–Weitzenkamp); Vaughan bounds their number up to $N$ by $N^{1-c\log\log\log N/\log\log N}$ |
| [236](https://www.erdosproblems.com/236) | OPEN | $\log_2n$ | is $\#\{k\ge0:n-2^k\text{ prime}\}=o(\log n)$? | Erdős: $\gg\log\log n$ for infinitely many $n$ |
| [852](https://www.erdosproblems.com/852) | OPEN | $\sim\log x$ | is $h(x)=o(\log x)$, and is $h(x)>(\log x)^c$, where $h(x)$ is the longest run of distinct consecutive prime gaps below $x$? | Brun's sieve gives $h(x)\to\infty$ |
| [1181](https://www.erdosproblems.com/1181) | OPEN | — | is $q(n,\log n)<(1-c)\log^2n$ for all large $n$, where $q(n,k)$ is the least prime not dividing $\prod_{i\le k}(n+i)$? | the primorial gives $(1+o(1))\log^2n$ |
| [141](https://www.erdosproblems.com/141) | OPEN | $k$ fixed | $k$ consecutive primes in AP, for every $k$? | verified for $k\le10$; infinitely many is open even for $k=3$ |
| [200](https://www.erdosproblems.com/200) | OPEN | $\sim\log N$ | this problem | §§2–7 |

Twelve problems that I picked myself are a small sample, so what follows is an observation, not evidence. The solved ones fix $k$ (Green–Tao and what follows from it), or are settled by a construction (#457, #1202), or involve about $\sqrt n$ values at once, where finiteness arguments have room (#1140, #1141). Every problem here whose pattern has length about $\log n$ (#1142, #236, #852 and #200) is open. None of the solved techniques transfers: a construction could only show that $L(N)$ is *large*, which the Hardy–Littlewood model of §4 predicts is false, and Green–Tao gives only the lower bound of §7.

The tracker also labels how each solved problem was settled. Over the 16 solved problems tagged `primes` the labels are sieve method 9, analytic estimate 6, explicit construction 6, elementary number theory 5, reduction to known theorem 3, and one each of covering/density argument, Fourier analysis and greedy algorithm. A problem can carry several labels, and they were written by the tracker's models (claude-sonnet-5 and gemini-3.8-flash), not by the site. Sieve methods are the commonest label, and §6 is the test of whether they transfer.

## 9. Machine-checked formalisation

Propositions 1 and 2 above are informal proofs, and the bound they imply is the one Erdős and Graham already credit to the prime number theorem. What was missing was a machine check. DeepMind's `FormalConjectures/ErdosProblems/200.lean` states the problem and the known variant, and both declarations are `sorry` — verified directly in the revision pinned here (`e5f4281`), at lines 44 and 52 of that file.

[`lean/`](lean/) closes the second one. The development is 673 lines of Lean 4 against Mathlib `0df444a3` (toolchain `v4.33.1`), and depends on `formal-conjectures` itself, so the statement being proved is *their* statement rather than a paraphrase of it.

| theorem | content |
|---|---|
| `structure_of_primeAP` | Proposition 1 |
| `mul_primorial_lt` | $(k-1)W_k<kN$ |
| `theta_lt_sharp` | $\theta(k)<\log N+1/(k-1)$ — Proposition 2 |
| `upper_unconditional` | $L(n)\le(1/\log2+o(1))\log n$, **no hypotheses** |
| `erdos_200_variants_upper` | `PNT →` $L(n)\le(1+o(1))\log n$ |

**The unconditional constant there is $1/\log 2 = 1.4427$, not $1$.** Mathlib does not contain the prime number theorem; the sharpest Chebyshev-type lower bound it has is `Chebyshev.theta_ge`, $\theta(n)\ge n\log2-\log(n+1)-2\sqrt n\log n$, whose leading constant is $\log 2=0.693$. So the unconditional result *in that project* is 44 % weaker than the literature. Recovering the literature constant needs $\theta(x)/x\to1$, which appears in the development as a `def PNT : Prop`, an *assumed* hypothesis, proved nowhere in `lean/`. Both results come from one transfer lemma, `eventually_longest_le`: any $\theta(k)\ge(c-o(1))k$ gives $L(n)\le(1/c+o(1))\log n$.

### 9.1 Discharging the hypothesis

That hypothesis is no longer needed. Mathlib still has no PNT, but [`PrimeNumberTheoremAnd`](https://github.com/AlexKontorovich/PrimeNumberTheoremAnd) (PNT+) proves `chebyshev_asymptotic : θ ~[atTop] id`, which *is* `PNT`. [`lean-pnt/`](lean-pnt/) bridges the two and obtains

$$\texttt{erdos\_200\_variants\_upper\_unconditional} \;:\; L(n)\le(1+o(1))\log n$$

with **no hypothesis** — DeepMind's statement, unconditionally, machine-checked. The port cost exactly one identifier (`Set.mem_ofPred_eq` → `Set.mem_setOf_eq`, renamed in Mathlib on 2026-07-09); `lean/` remains canonical and the shared modules are generated by `lean-pnt/scripts/sync.sh`.

Two caveats, both structural:

* PNT+ pins Lean `v4.32.2` and `formal-conjectures` pins `v4.33.1`. Two Mathlib revisions cannot coexist in one Lake project, so `lean-pnt/` cannot depend on `formal-conjectures` and the `Expr`-level fidelity check of `lean/` is unavailable there. The definitions are copied verbatim into `FCShim.lean` and `scripts/check_shim.py` re-downloads `formal-conjectures` `main` and fails on any difference. That is a weaker guarantee than `Expr` identity, and is stated rather than glossed.
* PNT+ itself contains two `sorry`s, at `Wiener.lean:323` and `:342`. Both are dead code: `prelim_decay_2` is referenced nowhere, and `prelim_decay_3` only by `decay_alt`, which is also referenced nowhere across its 27 source files. That is a textual argument; `#assert_axioms chebyshev_asymptotic` is what settles it, reporting `[propext, Classical.choice, Quot.sound]`. For contrast `#print axioms prelim_decay_3` reports `sorryAx`, so the check is not vacuous.

Consequently **no result in this repository now depends on an unformalised hypothesis.**

### 9.2 Checking the checkers

Two mechanical checks run as part of `lean/verify.sh`:

1. **Axioms.** `#assert_axioms` walks the transitive dependency graph of each of the 20 theorems with `Lean.collectAxioms` and errors if `sorryAx` appears. All 20 rest on exactly `{propext, Classical.choice, Quot.sound}`.
2. **Statement fidelity.** A `run_cmd` compares the elaborated `Expr` of the conclusion of `erdos_200_variants_upper`, with the `PNT` binder stripped, against the elaborated `Expr` of `Erdos200.erdos_200.variants.upper`. They are structurally identical. This check exists because check 1 cannot detect the failure that matters — a complete proof of the wrong statement.

A checker that has never failed is not evidence, so both were tested against planted defects. A `sorry` in `primorial_lt_two_mul` propagated `sorryAx` to six downstream theorems and failed the build; a `sorry` in `primorial_dvd` made `verify.sh` exit 1; and reassociating the conclusion to `log n * (1 + o n)` — mathematically identical, and still fully proved — was rejected by check 2 with `statement mismatch`. `lean-pnt/verify.sh` was tested the same way: altering `n < l` to `n ≤ l` in the copied `IsAPOfLengthWith`, and merging two non-contiguous declarations into one verbatim block, each made `check_shim.py` exit 1 with the offending line flagged.

### 9.3 What this does not establish

Nothing about the open question. Erdős 200 asks whether $L(N)=o(\log N)$; §5 argues the answer is conjecturally yes with decay $2/\log\log N$, and that remains a heuristic. The formalisation covers the known upper bound only. It is also not a new mathematical argument: it is Proposition 2 plus Chebyshev, transcribed, and the unconditional form now leans on PNT+, which is someone else's formalisation of someone else's theorem.

It also does **not** close the `sorry` upstream. `formal-conjectures` depends on Mathlib alone, and Mathlib has no PNT, so `lean-pnt/`'s proof cannot be submitted there as-is. It becomes submittable when either PNT+ updates to Lean `v4.33.1` or Mathlib absorbs the prime number theorem. Until then the unconditional result lives here, not upstream.

## 10. Reproducing

```sh
cc -O3 -march=native -funroll-loops -pthread -o src/apsearch src/apsearch.c
sh refs/fetch.sh                      # re-download the primary sources
python src/verify_records.py          # 26 records, 4 independent checks
python src/ref_bruteforce.py          # assumption-free reference
./src/apsearch 20 572945039351 10    # ~67 min on 10 threads
python src/heuristic.py               # singular series vs exact counts
python src/compare_counts.py          # measured vs predicted, 30 points
python src/calibration.py             # KS test of the Poisson model
python src/barrier.py                 # B(k,N) = 2^k k! C(k,N), kept for §6.2
python src/sieve_limits.py            # what sieves can and cannot say (§6, §8.1); ~70 s, ~3.1 GB peak memory
python src/analysis.py                # Tables 1 and 2
python src/final_check.py             # re-derives every claim above; exit 0 = pass; ~9 s
bash lean/verify.sh                   # Lean: axiom audit + statement fidelity; exit 0 = pass
bash lean-pnt/verify.sh               # Lean: same bound with no PNT hypothesis; exit 0 = pass
```

`sieve_limits.py` opens the tracker database `../tracker-and-solver/data/erdos.db` read-only if it is present; without it the §8.1 survey is recorded as null, and `final_check.py` says it skipped that check rather than failing.

`lean-pnt/verify.sh` builds `PrimeNumberTheoremAnd` from source. Measured cold: 3708 jobs in 3 min 17 s, and 3717 jobs for the full project, on the same machine. It also re-fetches `formal-conjectures` `main` over the network to run `scripts/check_shim.py`, so it needs connectivity; `lean/verify.sh` does not.

Timings on 10 threads (5 performance cores, Apple silicon): $k=17$ to $2\times10^{10}$, 139 s; $k=18$ to $3\times10^{10}$, 291 s; $k=19$ to $1.5\times10^{11}$, 390 s; $k=20$ to $5.73\times10^{11}$, 4027 s. Cost scales as $\varphi(W_k)M^2/(2(k-1)\cdot64)$ word-operations at roughly $4.9\times10^7$ per core-second, which is what makes $k\ge21$ impractical here. These are wall times on an otherwise idle machine; the $k\le16$ runs were re-measured later under heavy competing load and took up to $14\times$ longer, so treat any single timing as indicative only.
