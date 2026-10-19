# Notes — Erdős #1186

Technical record: normalisation, what was reproduced from the literature, the construction that produced the improvements, and what was tried and failed. Every bound in `results/claims.json` is re-derived by `src/final_check.py` in exact rational arithmetic, and exit 0 means every check passed. Other numbers here, such as UNSAT verdicts, run times and the $k \ge 9$ rows of §4, are not re-derived by the gate.

## 1. The problem and the normalisation

Erdős [Er80] asks for the optimal guaranteed coefficient $\delta_k$ such that every 2-colouring of $\{1,\dots,n\}$ contains at least $(\delta_k+o(1))n^2$ monochromatic $k$-APs, and for the analogous $\tilde\delta_k$ over $\mathbb{F}_p$. This is the largest such guaranteed coefficient, not the least; that wording is corrected here on October 7, 2026 against [LuPe12, eq. (10)]. Only upper bounds are attacked here: an explicit colouring with few monochromatic $k$-APs bounds $\delta_k$ from above. **No lower-bound work was attempted.**

Getting the normalisation right is the whole game, because a bound that looks like an improvement is usually a different quantity. Lu and Peng [LuPe12, eq. (10)] define

> $c_k = \frac{1}{2(k-1)} \lim_{n\to\infty} m_k([n])$

and state plainly, of their $m_k$,

> "Here we allow $k$-APs to be degenerated"

so the cyclic normalization counts **ordered pairs** $(a,b)$ with $a+jb$ monochromatic for $j<k$, including $b=0$. This repo's $\Phi_k$ uses that normalization for a particular construction, rather than taking the minimum over all colorings. The construction therefore gives upper bounds:

$$\delta_k \le \Phi^{\mathbb{Z}_m}_k/(2(k-1)), \qquad \tilde\delta_k \le \Phi_k/2 .$$

These inequalities correct the former equal signs; an arbitrary coloring's value is not the optimum. No stored numerical bound changes.

The $1/(2(k-1))$ is the area of $\{x,t\ge0,\ x+(k-1)t\le1\}$, the density of $k$-APs in $\{1,\dots,n\}$; the $/2$ over $\mathbb{F}_p$ is Wolf's unoriented convention. A random colouring gives $\delta_k\le1/((k-1)2^k)$ and $\tilde\delta_k\le2^{-k}$.

Degenerate progressions **must** be counted: a $\mathbb{Z}_m$-degenerate AP (one with $b\ne0$ that wraps) corresponds to a genuine non-degenerate AP in $\{1,\dots,n\}$. Dropping them answers a different question — this is why Grier's AP-free sets $G(3,2)$ are not comparable here.

The periodic-to-integer step is [LuPe12, Lemma 1]:

> "For any $k \ge 3$ and any positive integer $b$, we have $\lim_{n\to\infty} m_k([n]) \le m_k(\mathbb{Z}_b)$."

**Checked, not assumed.** `src/verify_two.c` brute-forces $\{1,\dots,n\}$ directly. The boundary term is $1/(2n)$, so (measured $-$ claimed) $\times n$ must sit at $\approx0.50$ and stay flat. The stored $k=6$ base word is the zipped quadratic-residue colouring of $\mathbb{Z}_{226}$ and is cyclic 6-AP-free, so on its own it gives $\delta_6 \le 1/2260$ (§3.2). For it the check measures **0.506, 0.501, 0.501** at $n = 20000, 60000, 120000$ (`results/verify_two_k6_226.log`). The earlier $k=6$ record ($m=3698$) gave 0.504, 0.502, 0.502.

## 2. Literature reproduced before anything was claimed

An evaluator that cannot recover a published value is not measuring the published quantity. Six reproductions, all exact:

| source | published | recomputed here |
|---|---|---|
| [PRS08, Thm 5] | $c_3 \le 117/2192$ | $117/2192$ ✓ |
| [LuPe12, Thm 1] $B_{20}$, $\mathbb{F}_p$ blocks | $m_4(\mathbb{Z}_p) \le 17/150$ | $17/150$ ✓ |
| [LuPe12] periodic, $20\mid n$ | $m_4(\mathbb{Z}_n) \le 0.09$ | $m_4(\mathbb{Z}_{20}) = 9/100$ — exactly $0.09$ ✓ |
| [LuPe12] periodic, $22\mid n$ | $m_4(\mathbb{Z}_n) \le 0.086777$ | $m_4(\mathbb{Z}_{22}) = 21/242 = 0.0867769$ ✓ |
| [LuPe12] $B_{74}$ | $m_5(\mathbb{Z}_{74}) = 73/2738$ | $73/2738$ ✓ |
| [RL12, Table 1], two colours | $W(2,k) > 3703, 11495, 41265, 103474, 193941, 638727$ for $k = 7, \dots, 12$ | all six, as $(k-1)m+1$ for the largest good $m$ (`src/rabung.c`, §4) ✓ |

The two threshold values are the tightest of these: $9/100$ is $0.09$ on the nose, and $21/242$ reproduces all six quoted digits of $0.086777$. (The best $m_4$ found here over moduli divisible by 22 is in fact $61/726 = 0.084022$ at $m=66$, better than their $\mathbb{Z}_{22}$ value — but this is an intermediate quantity, and after dividing by $2(k-1)=6$ it still does not reach $1/72$; see §6.)

The $B_{74}$ check is the sharp one, because it reproduces the *structure* and not just the value. Lu and Peng write:

> "All monochromatic 5-APs in $B_{74}$ are degenerated ones. Among them there are 74 5-APs with $d = 0$ and 72 5-APs with $d = 37$. This coloring gives $m_5(\mathbb{Z}_{74}) = 146/74^2 = 73/2738$."

`src/verify_zm.c` prints the $d$-breakdown and returns exactly 74 at $d=0$ and 72 at $d=37$.

Targets: $c_4 \le 1/72$ and $c_5 \le 1/304$ [LuPe12, eqs. (12)–(13)], and $\lim m_4(\mathbb{Z}_n)\le1/12$, $\lim m_5(\mathbb{Z}_n)\le1/38$ [LuPe12, Thm 5]. Over $\mathbb{F}_p$: $m_4(\mathbb{Z}_p)\le17/300+o(1)$ [Thm 1] and $m_5(\mathbb{Z}_n)\le3629/65712+o(1)$ for odd $n$ [Thm 4].

All quotations above were read from `refs/lupe12.pdf` directly (text extracted with `pypdf`), not from a secondary source.

## 3. The wildcard recursion

This is the ingredient behind every improvement here, and it is now a theorem (§3.1).

Lu and Peng found their $k=4$ bound by noticing that the two halves of a good $\mathbb{Z}_{22}$ colouring differ in one bit, giving an 11-bit pattern with a wildcard. Running the same comparison one level up — on the optima found at $m=88$ and $m=132$ for $k=5$ — both reduce to the **same** 44-bit base with exactly four free positions $\{3,14,25,36\}$, the coset $x\equiv3\pmod{11}$.

**Bases.** A *base* is a 2-colouring $c$ of $\mathbb{Z}_b$ together with a proper coset $F = r + g\mathbb{Z}_b$ of the subgroup of order $f<b$, where $g = b/f$. It is *valid* for $k$ when

- **(H1)** every $k$-AP $a, a+d, \dots, a+(k-1)d$ of $\mathbb{Z}_b$ with $d \not\equiv 0 \pmod b$ and no term in $F$ is non-monochromatic, and
- **(H2)** every $k$-AP with $d \not\equiv 0 \pmod b$ and terms both inside and outside $F$ is non-constant on its terms outside $F$.

The colours on $F$ enter neither condition. (H2) is exactly the test `cosetprobe.py:safe_coset`: an AP is *dangerous* when its outside terms agree, because colouring the inside to match would make it monochromatic. Both conditions are finite checks over $\mathbb{Z}_b$.

**Theorem.** Let $(c, F)$ be a valid base with $f<b$, $t \ge 1$, and $h$ any 2-colouring of $\mathbb{Z}_{ft}$. Colour $\mathbb{Z}_{bt}$ by $C(x) = c(x \bmod b)$ when $x \not\equiv r \pmod g$, and $C(r + gu) = h(u)$ otherwise. Then

$$\#\mathrm{mono}_k(\mathbb{Z}_{bt}, C) \;=\; (b-f)\,t^2 + \#\mathrm{mono}_k(\mathbb{Z}_{ft}, h),$$

where $\#\mathrm{mono}$ counts pairs $(a,d)$ with $d=0$ included, so that $\Phi = \#\mathrm{mono}/m^2$. Hence

$$m_k(\mathbb{Z}_{bt}) \;\le\; \frac{(b-f) + f^2\, m_k(\mathbb{Z}_{ft})}{b^2}, \qquad\text{fixed point } \frac{1}{b+f}, \qquad \delta_k \le \frac{1}{2(k-1)(b+f)} .$$

Earlier versions of these notes called the recursion "verified but not proved in general" and its limits conjectural. §3.1 proves it for every $k$, $b$, $f$ and $t$.

**The formula reproduces both published constants exactly, and the bases behind them are quadratic-residue colourings (§4).** "vs random" is the bound divided by the random-colouring value $1/((k-1)2^k)$, which equals $2^{k-1}/(b+f)$.

| $k$ | base | $b$ | $f$ | limit $1/(b+f)$ | $\delta_k \le$ | vs random | status |
|---|---|---|---|---|---|---|---|
| 4 | Lu–Peng $B_{11}$ = QR colouring of $\mathbb{Z}_{11}$ | 11 | 1 | $1/12$ | $1/72$ | 0.667× | [LuPe12] eq. (12) ✓ |
| 5 | Lu–Peng $B_{37}$ = QR colouring of $\mathbb{Z}_{37}$ | 37 | 1 | $1/38$ | $1/304$ | 0.421× | [LuPe12] eq. (13) ✓ |
| 5 | $B_{44}$, found by search here | 44 | 4 | $1/48$ | $\mathbf{1/384}$ | 0.333× | **new**; optimal in the family (§5) |
| 6 | $B_{86}$, found by search here | 86 | 2 | $1/88$ | $1/880$ | 0.364× | superseded |
| 6 | zipped QR colouring of $\mathbb{Z}_{226}$ | 226 | 2 | $1/228$ | $\mathbf{1/2280}$ | 0.140× | **new** |
| 7 | QR colouring of $\mathbb{Z}_{617}$ | 617 | 1 | $1/618$ | $\mathbf{1/7416}$ | 0.104× | **new** |
| 8 | zipped QR colouring of $\mathbb{Z}_{1642}$ | 1642 | 2 | $1/1644$ | $\mathbf{1/23016}$ | 0.078× | **new** |

Lu and Peng used $b=37$ because their search found $B_{37}$; the largest 5-AP-free modulus is 44, and it carries a free coset of size 4 rather than 1. At $k=6$ the step is larger. The largest cyclic 6-AP-free modulus is 226 (§5), and the zipped quadratic-residue colouring of $\mathbb{Z}_{226}$ is a valid base with $F=\{0,113\}$, so $1/880$ drops to $1/2280$.

**Stored bases.** `results/safebase_k{4..8}.json` holds one record per base, with fields `k, b, f, r, word, source` and an informational `bound`, assembled and checked by `src/bases.py`. The gate never reads `bound`; it recomputes the bound from $b$ and $f$. Each `word` also colours $F$, so that the whole word is cyclic AP-free; the extension lemma of §5 guarantees this is possible. The $k=5$ base is `10010111010001011011110110100010111010010000` with $F=\{3,14,25,36\}$. The superseded $B_{86}$ ($F=\{33,76\}$) is kept as a record.

### 3.1 Proof

Let $S = \{x \in \mathbb{Z}_{bt} : x \equiv r \pmod g\}$ be the lift of $F$. It has $ft$ elements, and $u \mapsto r + gu$ is a bijection $\mathbb{Z}_{ft} \to S$. Membership in $S$ depends only on $x \bmod b$, and so does the colour $C(x)$ off $S$. Sort the pairs $(a, d) \in \mathbb{Z}_{bt}^2$ into three cases:

1. **$d \not\equiv 0 \pmod g$.** Reducing mod $b$ gives a $k$-AP of $\mathbb{Z}_b$ with nonzero difference. It cannot lie inside $F$, because consecutive terms differ by $d \not\equiv 0 \pmod g$. If it avoids $F$, (H1) makes it non-monochromatic. If it meets $F$, (H2) makes its outside terms disagree, whatever $h$ does. None of these pairs is monochromatic.
2. **$d \equiv 0 \pmod g$ and $a \notin S$.** All terms are congruent mod $g$, so all lie outside $S$ and carry colours of $c$. If $d \not\equiv 0 \pmod b$, the reduction avoids $F$ and has nonzero difference, so (H1) forbids a monochromatic AP. If $d \equiv 0 \pmod b$, every term has colour $c(a \bmod b)$, so the AP is monochromatic. There are $(b-f)t$ choices of $a$ and $t$ of $d$, giving $(b-f)t^2$.
3. **$d \equiv 0 \pmod g$ and $a \in S$.** The map $(u, e) \mapsto (r + gu,\ ge)$ is a bijection from $\mathbb{Z}_{ft}^2$ onto these pairs, and $a + jd = r + g(u + je)$. So the AP is monochromatic under $C$ exactly when $(u, e)$ is monochromatic under $h$, which contributes $\#\mathrm{mono}_k(\mathbb{Z}_{ft}, h)$. ∎

**Chain.** Start from the stored word on $\mathbb{Z}_b$, which has $\Phi_0 = 1/b$. Apply the theorem with $t = m_s/f$ and $h$ equal to the previous level. Then $m_{s+1} = b\,m_s/f$ and $\Phi_{s+1} = (b-f)/b^2 + (f/b)^2\,\Phi_s$, which converges to $1/(b+f)$ from above with ratio $(f/b)^2$. By [LuPe12, Lemma 1], $\delta_k \le \Phi_s/(2(k-1))$ at every level, so $\delta_k \le 1/(2(k-1)(b+f))$. The bound is an infimum; every finite level lies strictly above it.

**What the gate checks.** For every stored base, check C in `src/final_check.py` verifies (H1) and (H2) by a plain loop over $\mathbb{Z}_b$ and recomputes the bound from $b$ and $f$; the stored `bound` is never read. Check H then re-tests (H2) with the independently written `cosetprobe.safe_coset`, confirms with `verify_zm` that the stored word is cyclic AP-free, and re-measures the identity itself on explicit lifts, at $t=2$ with $h \equiv 0$ and at $t=3$ with a random non-constant $h$. Planted defects are in §10.

### 3.2 Finite levels, certified by direct count

The theorem makes finite levels unnecessary for the bounds. They remain as a cross-check that does not depend on it: each is one explicit colouring whose monochromatic APs are counted directly.

| $k$ | certified by direct count | limit (theorem) |
|---|---|---|
| 5 | $m_5(\mathbb{Z}_{5324}) = 13421/644204 \Rightarrow \delta_5 \le 13421/5153632$ | $\delta_5 \le 1/384$ |
| 6 | $m_6(\mathbb{Z}_{226}) = 1/226 \Rightarrow \delta_6 \le 1/2260$ | |
| 6 | $m_6(\mathbb{Z}_{25538}) = 12657/2885794 \Rightarrow \delta_6 \le 12657/28857940$ | $\delta_6 \le 1/2280$ |

The zipped word on $\mathbb{Z}_{226}$ is cyclic 6-AP-free, so only the 226 pairs with $d=0$ are monochromatic (`verify_zm` reports "b=0 contributes 226, b!=0 contributes 0"). That gives $\delta_6 \le 1/2260$ with no recursion at all, 2.57 times below this repo's previous $k=6$ value $1807/1590140$ ($m=3698$). One level of the recursion ($t=113$, $h$ = the stored word) gives $\mathbb{Z}_{25538}$. `verify_zm` counts exactly $2860482 = 224\cdot113^2 + 226$ monochromatic pairs there, which is the theorem's count, and $2911332 = 224\cdot113^2 + 226^2$ for the control $h \equiv 0$ (`results/verify_lift_k6_25538.log`). The resulting bound sits $3.0\times10^{-10}$ above $1/2280$. The $k=5$ level and $\mathbb{Z}_{226}$ were each checked by `exact.py`, `verify_zm.c` and `verify_two.c`, which share no code. $\mathbb{Z}_{25538}$ is beyond `exact.py`'s practical range and was checked by `verify_zm` alone.

Before the proof, the recursion was confirmed against explicit constructions at $t = 2, 3, 4, 11, 121$ for $k=5$ ($m = 88, 132, 176, 484, 5324$) and at $t=43$ for $k=6$ ($m=3698$). The $t=4$ case refutes the obvious guess. The naive law $(10t+1)/(484t)$ assumes the inner colouring is 5-AP-free, and it predicts $41/1936 = 0.02117769$ at $m=176$. But the true minimum over $\mathbb{Z}_{16}$ is $m_5(\mathbb{Z}_{16}) = 5/64$, not $1/16$, and the construction attains $15/704 = 0.02130682$. The naive law claims a value **below what is achievable**, so it is not a valid bound. The theorem gives $15/704$ exactly. `src/wildcard.py:predicted` implements the theorem, and the naive law survives only in its docstring, as a warning.

## 4. Where the bases come from: quadratic residues

For a prime $p$, colour $x \in \mathbb{Z}_p^*$ by its quadratic character: $c(x) = 1$ iff $x$ is a non-residue. This is the colouring behind the cyclic van der Waerden certificates of [RL12]. Because $c(dz) = c(z) \oplus c(d)$ for $z \ne 0$, every AP $a + jd$ with $d \ne 0$ is the window $y, y+1, \dots, y+k-1$ (where $y = a/d$) rescaled by $d$, with its colours kept or all flipped. So $\mathbb{Z}_p$ has no monochromatic $k$-AP with $d \ne 0$ iff

- (i) no $k$ consecutive elements of $1, \dots, p-1$ share a colour, and
- (ii) no window of $k$ consecutive residues through 0 is constant on its nonzero elements.

[RL12, p. 5] uses the same multiplicativity to cut the check from quadratic to linear time.

**Criterion (i) is (H1), and criterion (ii) is (H2), for $F = \{0\}$.** Every prime that passes the test is therefore a valid base with $f = 1$, giving $\delta_k \le 1/(2(k-1)(p+1))$. The zipped colouring of $\mathbb{Z}_{2p}$ [HHLM07], read through $\mathbb{Z}_{2p} \cong \mathbb{Z}_p \times \mathbb{Z}_2$, is $z(x) = c(x \bmod p) \oplus [x \text{ even}]$ for $x \ne 0, p$, with $z(0) \ne z(p)$. Its criteria are (H1) and (H2) for $F = \{0, p\}$. The condition $z(0) \ne z(p)$ concerns the difference $p$; it matters for full AP-freeness but not for the base. Every good zip is therefore a base with $f = 2$, giving $\delta_k \le 1/(2(k-1)(2p+2))$. The argument is the rescaling above, applied to odd and even differences separately, and the header of `src/rabung.c` spells it out. Nothing downstream relies on it: the gate verifies (H1) and (H2) directly for $k \le 8$ (checks C and H), and `src/verify_base.c` does so for the larger words below.

**Lu–Peng's bases are quadratic-residue colourings.** $B_{11} = (1,1,1,0,1,*,0,1,0,0,0)$ and the 37-bit $B_{37}$ of [LuPe12], read cyclically from the wildcard, are exactly the colourings above for $p = 11$ and $p = 37$ (1 = non-residue), with the wildcard at 0. Their Property 2, "No matter which bit-value the '∗' takes, $B_{37}$ contains no non-degenerated monochromatic 5-APs of $\mathbb{Z}_{37}$", is (H1)+(H2) at $F = \{0\}$. Their Lemma 7, the $k=5$ step of the recursion, comes with "The proof is omitted"; Lemma 6, the $k=4$ step, is proved there. §3.1 covers every $k$ and every $f$.

**Every two-colour lower bound in [RL12, Table 1] is reproduced.** `src/rabung.c search 3 12 60000` takes about a second and finds the largest good $p$ and $2p$ below $6\times10^4$ for each $k$. [RL12] write $W(\text{colours}, \text{length})$, and their bound $W(2,k) > (k-1)m + 1$ matches their table at every $k = 7, \dots, 12$, including the $k=11$ entry they mark as correcting an earlier misreported value.

| $k$ | $m$ | kind | $f$ | [RL12] $W(2,k)$ | $\delta_k \le 1/(2(k-1)(m+f))$ | vs random | (H1)+(H2) checked by |
|---|---|---|---|---|---|---|---|
| 4 | 11 | QR | 1 | $= 35$ | $1/72$ | 0.667× | gate |
| 5 | 37 | QR | 1 | $= 178$ | $1/304$ | 0.421× | gate |
| 6 | 226 | zip | 2 | $= 1132$ | $1/2280$ | 0.140× | gate |
| 7 | 617 | QR | 1 | $> 3703$ | $1/7416$ | 0.104× | gate |
| 8 | 1642 | zip | 2 | $> 11495$ | $1/23016$ | 0.078× | gate |
| 9 | 5158 | zip | 2 | $> 41265$ | $1/82560$ | 0.050× | `verify_base` |
| 10 | 11497 | QR | 1 | $> 103474$ | $1/206964$ | 0.045× | `verify_base` |
| 11 | 19394 | zip | 2 | $> 193941$ | $1/387920$ | 0.053× | `verify_base` |
| 12 | 58066 | zip | 2 | $> 638727$ | $1/1277496$ | 0.035× | `verify_base` |

For $k = 4, 5, 6$, [RL12] list the exact value of $W$. For $k=5$ the QR family peaks at 37, and search does better (44). Rows with $k \ge 9$ are **not** in `results/claims.json` and are not gate-checked. Their words were checked by `src/verify_base` (both conditions, 0 violations) and by `src/verify_zm` (cyclic AP-free), with output in `results/verify_base_k9_12.log`. These are fixed-size constructions. The ratio to random falls from 0.667× at $k=4$ to 0.035× at $k=12$, though not monotonically: it rises from 0.045× at $k=10$ to 0.053× at $k=11$. That says nothing about asymptotics, and no asymptotic formula for $\delta_k$ is claimed.

## 5. How far the family can go

Two lemmas bound the family, and an exact SAT scan settles what is left.

**Cap lemma.** If $\mathbb{Z}_m$ has a 2-colouring with no monochromatic $k$-AP of nonzero difference, then $(k-1)m < W(k;2)$. *Proof.* Repeat the word $k-1$ times along $\{0, \dots, (k-1)m-1\}$. A $k$-AP there with difference $d \ge 1$ reduces mod $m$ to an AP of $\mathbb{Z}_m$. If $m \nmid d$, that AP is non-monochromatic; if $m \mid d$, the last term is at least $(k-1)m$, outside the range. ∎ [RL12, Table 1] lists $W(3;2)=9$, $W(4;2)=35$, $W(5;2)=178$ [SS78] and $W(6;2)=1132$ [KP08], so the caps are $4, 11, 44, 226$ for $k = 3, 4, 5, 6$. All four are attained: by `0011`, the QR colouring of $\mathbb{Z}_{11}$, $B_{44}$, and the zipped $\mathbb{Z}_{226}$.

**Extension lemma.** Any valid base $(c, F)$ on $\mathbb{Z}_b$ can be recoloured on $F$ to be cyclic $k$-AP-free, so $b$ is at most the cap. *Proof.* Take a coset $s + g\mathbb{Z}_b \ne F$ and read it as the colouring $u \mapsto c(s + gu)$ of $\mathbb{Z}_f$. Its APs have differences $ge$ with $e \not\equiv 0 \pmod f$ and avoid $F$, so by (H1) it has no monochromatic $k$-AP of nonzero difference. Copy it onto $F$. APs avoiding $F$ are then fine by (H1). APs meeting both $F$ and its complement are fine by (H2), whatever $F$'s colours. APs inside $F$ are APs of the copied coset. ∎ So every stored base word is fully AP-free, and $b \le 4, 11, 44, 226$ for $k = 3, 4, 5, 6$.

**Exhaustive scan.** (H1) and (H2) are clauses over the colours outside $F$: each outside-set must be non-constant. `src/safesat.py` therefore decides each $(k, b, f)$ exactly. A model is an explicit base, re-checked by `violations()` and `safe_coset`. UNSAT is a proof that none exists. FORCED means some AP has a single outside term, which no colouring can make non-constant. By the extension lemma, scanning $b$ up to the cap suffices. The scans for $k \le 5$ go further, as a margin; the $k=6$ scan stops at the cap. Exceeding $b + f = T$ needs $b > 2T/3$, since $f \le b/2$.

| $k$ | cap | scanned | candidates $(b,f)$ | result | largest $b+f$ | best bound from the family |
|---|---|---|---|---|---|---|
| 3 | 4 | $b \le 60$, $b+f > 4$ | 199 | 110 FORCED, 89 UNSAT | 4, at $(3,1)$ | $\delta_3 \le 1/16$ = random |
| 4 | 11 | $b \le 150$, $b+f > 12$ | 613 | 190 FORCED, 423 UNSAT | 12, at $(11,1)$ | $\delta_4 \le 1/72$ = Lu–Peng |
| 5 | 44 | $b \le 200$, $b+f > 48$ | 772 | 234 FORCED, 538 UNSAT | 48, at $(44,4)$ | $\delta_5 \le 1/384$ |
| 6 | 226 | $b \le 226$, $b+f > 228$ | 103 | 31 FORCED, 72 UNSAT | 228, at $(226,2)$ | $\delta_6 \le 1/2280$ |

The $k=3$ row is a control with a known answer. A base with $b+f \ge 5$ would give $\delta_3 \le 1/20 < 1675/32768$, contradicting the lower bound of [PRS08, Thm 4]. So the scan must find nothing, and it finds nothing. For $k=4$ the family cannot beat Lu–Peng. For $k=5$, $1/384$ is optimal within it, and for $k=6$, $1/2280$ is.

**The hard case was $b = 221 = 13 \cdot 17$.** Of the 103 candidates for $k=6$, 101 were decided without symmetry breaking. The plain scan was stopped by hand after 70 minutes, 66 of them spent on $(221,13)$ without an answer, and a separate run with a budget of $2\times10^7$ conflicts gave up on $(221,13)$ with UNKNOWN. With `--symbreak 24` (next paragraph), `cadical195` proved $(221,13)$ UNSAT in 1763 s and `kissat404` proved $(221,17)$ UNSAT in 2356 s. Each was then re-run with the other solver: `cadical195` finds $(221,17)$ UNSAT in 2593 s, and a `kissat404` run on $(221,13)$ was stopped after 72 minutes without an answer, so that instance rests on one solver. The commands and times are in the header of `results/safesat_k6.log`.

**Symmetry breaking.** The maps $x \mapsto ux + s$, with $u$ a unit mod $b$ and $s \in F$, fix $F$ and permute the clauses. They form a group of order $f\varphi(b)$, which is 2496 or 3264 at $b = 221$. `safesat.py --symbreak L` adds lex-leader constraints for every such map on the first $L$ free positions. This is sound because the lex-least member of any orbit satisfies all of them, together with the colour-swap clause. `safesat.py --self-test` checks that it changes no answer on all 588 instances $(k,b,f)$ with $k = 3, 4, 5$ and $b < 30, 60, 90$, 24 of them SAT. The test has teeth: lex-leader constraints for random non-symmetries cut 13 of those 24 SAT instances, and substituting random permutations for the group drops agreement to 578 of 588, so the self-test fails. The $(k,b,f)$ of every stored base with $k \le 6$ stays SAT under symmetry breaking, with both `cadical153` and `kissat404`. All three are in `results/safesat_selftest.log`.

## 6. Why $\delta_4$ and $\delta_3$ were **not** improved

Stated plainly because it is the main negative result.

**$\delta_4$: reproduced exactly, and the family cannot go further.** The recursion on the quadratic-residue base $\mathbb{Z}_{11}$ ($f=1$) gives exactly Lu–Peng's $\delta_4 \le 1/72$, now as a theorem. The best finite level certified by direct count is $2917/209088 = 0.01395106$ (two-scale, $m=66$, $B=4$), 0.45 % above it. No other base does better: the extension lemma forces $b \le 11$, and `safesat.py` finds no base with $b+f > 12$ for any $b \le 150$ (§5). At $m = 11t$ the optimum is $\Phi = (10t+1)/(121t)$, which is the theorem's count with a 4-AP-free inner colouring and matches [LuPe12, Lemma 6]. So the finite levels approach $1/12$ strictly from above, as in [LuPe12, Thm 5]: $m_4(\mathbb{Z}_{11^s}) \le 1/12 + 1/(12\cdot11^{2s-1})$. This is consistent with [LuPe12, Conjecture 1]:

> "$\inf\{m_4(\mathbb{Z}_n) : n$ is not divisible by $4\} = 1/12$."

Beating $1/72$ needs a construction outside this family.

**$\delta_3$: reproduced, not improved.** $117/2192$, equal to [PRS08, Thm 5]. Both periodic families are weak at $k=3$. The best periodic $\mathbb{Z}_m$ value is $\Phi_3 = 1/16$ at $m=12$, which is the random value. No base has $b+f > 4$ (§5), so the wildcard family also stops at random.

## 7. Other negative results

Recorded so they are not re-run.

- **Larger free cosets at $b=44$ do not exist.** $f$ must divide 44. Direct search at $m=88$ gives $f=4 \to 21/968$ (best), $f=11 \to 87/3872$ and $f=22 \to 43/1936$. The exact structural test agrees: only $f=2$ and $f=4$ are safe, and $f=4$ wins. §5 extends this to every $b$, so $1/48$ is the ceiling of the family for $k=5$.
- **At $b=86$ for $k=6$, $f=2$ is the best coset.** The divisors are $1, 2, 43, 86$. The structural test is exhaustive over all of them and finds $f=2$ (at $r=33$) the only safe coset of size $>1$. So the ceiling at $b=86$ is $1/88$; the gain to $1/228$ came from a larger base, not a larger coset.
- $f=1$ cosets *do* lower $\Phi$ at finite $t$ (e.g. $87/3872 < 1/44$) without being safe. Improving a finite instance is not the same as supporting the recursion; do not confuse the two.
- Random restarts are useless past $m\approx100$. Tiling a $\mathbb{Z}_m$ word to $\mathbb{Z}_{tm}$ preserves $\Phi$ exactly, so seeding from the tiled record is what makes large $m$ reachable (`src/bigm.py`). SAT made this moot for $k=6$: `src/apfree.py 6 226` finds a cyclic 6-AP-free $\mathbb{Z}_{226}$ in 0.9 s.
- **Cyclic 6-AP-free moduli, $87 \le m \le 240$** (`results/apfree_k6.json`; `src/apfree.py` with a budget of $3\times10^6$ conflicts per modulus). The range starts above the earlier record 86. The scan finds SAT at $m = 127, 139, 226$ and UNSAT at 79 moduli, and leaves 72 undecided. It does far better on even moduli. Of the 70 even $m \le 226$, 61 are UNSAT, 226 is SAT, and the 8 undecided are exactly $2p$ for the eight primes $79 \le p \le 109$. Of the 70 odd ones, 58 are undecided and 127 and 139 are SAT. The 10 odd UNSAT moduli up to 226 are all multiples of 3, but 14 other odd multiples of 3 are undecided. Above 226 it finds 8 UNSAT, 6 undecided and no model. That no modulus above 226 works follows from the cap lemma (§5), not from this scan.
- **Row lemma; it did not settle $b=221$.** If $\gcd(f, g) = 1$, then $\mathbb{Z}_b \cong \mathbb{Z}_g \times \mathbb{Z}_f$ with $F = \{0\} \times \mathbb{Z}_f$. APs with difference $\equiv 0 \pmod f$ stay in a row $\mathbb{Z}_g \times \{y\}$, which meets $F$ once, so every row of a base is a $(g,1)$ base. The stored bases conform. $B_{44}$'s rows are $(11,1)$ bases for $k=5$, and the rows of the $k=6$ bases $\mathbb{Z}_{86}$ and $\mathbb{Z}_{226}$ and of the $k=8$ base $\mathbb{Z}_{1642}$ are $(43,1)$, $(113,1)$ and $(821,1)$ bases, with 0 violations in all 10 rows. At $b=221$ it does not help: for $k=6$ there are 888 $(13,1)$ bases and 3810 $(17,1)$ bases. A cube-and-conquer attempt built on the lemma fixed one row to each orbit representative of the row words under $x \mapsto vx$ and colour swap, giving 41 cubes for $(221,17)$ and 123 for $(221,13)$. It was stopped when the first of the 41 $(221,17)$ cubes was still undecided after about 8.5 minutes, against 39 to 43 minutes for the whole instance with symmetry breaking. Its script was not added to the repo.

## 8. Provenance caveats

- **[Wo10] and [BCG10] are paywalled and were not read.** Their numbers ($c_4 < 0.0172202$, $c_5 < 0.005719619$) appear here only as quoted by [LuPe12]. They are weaker than the bounds used as targets, so nothing depends on them, but they are second-hand and are marked as such.
- [CCS07]'s typeset Theorem 4.4 says "4-coloring" while its proof uses two colour classes. Almost certainly a typo; the two-colour reading is the one used here.
- **$k\ge6$: no published bound was located, for either problem.** The $k \ge 6$ entries are therefore reported as "no published bound found", not as records. That the sources read contain nothing is not proof that nothing exists. The recursion theorem, the quadratic-residue reading of Lu–Peng's bases, and the values for $\delta_6$, $\delta_7$, $\delta_8$ may already be known to specialists; [BCG10] in particular was not read.
- **The optimality statements of §5 rest on computer results; the bounds do not.** The scans cover $b \le 150$, $200$ and $226$ directly, and larger $b$ are excluded by the cap. The cap uses $W(5;2)=178$ [SS78] and $W(6;2)=1132$ [KP08], neither of which was read. Both values are taken from [RL12, Table 1], and [KP08]'s from its title as well.
- **Most UNSAT verdicts behind §5 are now externally proof-checked, but the k=6 optimality statement still has four unchecked cases.** `src/prove_unsat.py` rebuilds each CNF with `safesat.encode_cnf`, runs Homebrew CaDiCaL 3.0.1 with `--lrat`, and checks the proof with `lrat-check` from `drat-trim` commit `2e3b2dc0ecf938addbd779d42877b6ed69d9a985`. `results/lrat_check.json` records 1118 accepted LRAT proofs and no rejected completed proofs: 89/89 UNSAT instances for $k=3$, 423/423 for $k=4$, 538/538 for $k=5$, and 68/72 for $k=6$. The unchecked $k=6$ instances are $(203,29)$, $(217,31)$, $(221,13)$ and $(221,17)$, each timed out under a 300 s per-instance cap in the LRAT run. The two $b=221$ encodings include lex-leader symmetry-breaking clauses; LRAT checks those encoded CNFs, but soundness of adding those clauses rests on the symmetry argument above rather than on LRAT. FORCED answers need no solver. No bound depends on an UNSAT answer: every bound in `results/claims.json` comes from an explicit word that the gate re-checks.

## 9. A subagent's numbers were wrong — do not trust unverified ones

A literature agent returned two colouring certificates and several claimed values. On checking: its "88-bit" string had **92** characters and its "176-bit" string had **140**, so neither certificate could be evaluated at all. Its claim that $m_5(\mathbb{Z}_{44}) = 7/242$ is false for the word used here, which is antipodal ($c(x+22) = 1-c(x)$), so no difference divisible by 11 can be monochromatic. Its best claimed values ($85/3872$ at $m=88$, $333/15488$ at $m=176$) are both beaten by direct search here ($21/968$, $15/704$).

Its one useful contribution was pointing out that Lu–Peng's bound comes from a recursion rather than a fixed modulus — a structural hint, which was then verified against the PDF.

On a second request it re-sent the colourings as explicit integer lists. Those *did* have the stated lengths and were evaluable, and its two best candidates were checked here in exact arithmetic before publication:

| agent's claim | verified $\Phi$ | this repo | outcome |
|---|---|---|---|
| $k=5$, $m=220$ | $51/2420 = 0.02107438$ | $13421/644204 = 0.02083346$ | does not beat |
| $k=6$, $m=344$ | $169/14792 = 0.01142509$ | $1807/159014 = 0.01136378$ | does not beat |

It independently arrived at $\mathbb{Z}_{86}$ as the best 6-AP-free modulus, which corroborates the base used in §3. The standing lesson is unchanged: every number it supplied was re-derived here before being believed, and the ones that mattered were checked before anything was published.

*Later correction.* $\mathbb{Z}_{86}$ is not the best 6-AP-free modulus. SAT finds cyclic 6-AP-free colourings of $\mathbb{Z}_{127}$, $\mathbb{Z}_{139}$ and $\mathbb{Z}_{226}$ (§7), and none exists above 226 (§5). Agreement between two searches was therefore no evidence that 86 was optimal. The $k=6$ comparison above is against the superseded $m=3698$ value and is kept as it stood.

## 10. Validation

`src/final_check.py` re-derives every claim from the stored artifacts and asserts the contents of `results/claims.json` against them; it does not read any cached score. The full run, `python src/final_check.py`, passes all 117 checks with exit 0 in 179 s (`results/final_check_full.log`).

**The gate was tested by planting defects, not merely observed to pass.** A checker that has only ever seen good input is untested. Each defect below was planted in a stored file, the gate was run as `final_check.py --fast --no-brute`, and the file was put back. The clean run and the run after restoring both pass all 100 checks with exit 0, and the three files touched were byte-identical to the originals afterwards.

| planted defect | exit | FAIL lines |
|---|---|---|
| bit 1 of the $\mathbb{Z}_{226}$ base flipped (outside $F$) | 1 | 5: re-certification (`H1/H2 or shape`), `safe_coset` (6 dangerous progressions), AP-freeness (284 monochromatic pairs), and the identity at $t=2$ ($1128 \ne 224\cdot2^2+16$) and at $t=3$ ($2532 \ne 224\cdot3^2+12$) |
| bit 0 of the same base flipped (inside $F=\{0,113\}$) | 1 | 1: AP-freeness (228 monochromatic pairs) |
| the same base's coset moved to $r=1$ | 1 | 4: re-certification, `safe_coset` (112 dangerous progressions), and the identity at $t=2$ ($1136 \ne 224\cdot2^2+16$) and at $t=3$ ($2532 \ne 224\cdot3^2+12$) |
| $\delta_6$ claim changed to $1/2300$ | 1 | `[FAIL] delta_6: certified 1/2280 <= claimed 1/2300` |
| $\delta_5$ claim changed to $1/400$ | 1 | `[FAIL] delta_5: certified 1/384 <= claimed 1/400` |
| one bit of the stored $m=484$ word in `results/wildcard_k5_deep.json` flipped, its recorded `phi` left intact | 1 | `[FAIL] zm   k=5: 12 words re-certify   1 bad: [('484', '633/234256 != 111/42592')]` |

Two rows need comment. Colours on $F$ enter neither (H1) nor (H2), so the bit-0 defect is invisible to every base check except AP-freeness. That check was added when a reading of the gate showed it was missing, and this defect was planted to test it. The $m=484$ row shows that $\Phi$ is recomputed from the colouring rather than read from the file, so a corrupted certificate cannot pass by carrying a correct-looking number.

**The SAT tooling was tested the same way** (`results/safesat_selftest.log`). `safesat.py --self-test` decides all 588 instances $(k,b,f)$ with $k = 3, 4, 5$ and $b < 30, 60, 90$ with and without symmetry breaking, and every answer agrees; 24 are SAT. Lex-leader constraints for random non-symmetries turn 13 of the 24 UNSAT, so the comparison can see a changed answer. Replacing the symmetry group by random permutations drops agreement to 578 of 588, and the self-test exits 1.

**External LRAT proof checking now covers most UNSAT scan answers** (`results/lrat_check.json`). The checker path was validated on $(k,b,f)=(4,12,2)$: `lrat-check` accepts the correct CaDiCaL LRAT proof, rejects the proof after corrupting a hint/clause, and rejects the same proof against a CNF with one clause dropped. The checked corpus is 89 $k=3$ UNSAT instances, 423 $k=4$ instances, 538 $k=5$ instances, and 68 $k=6$ instances; the four $k=6$ timeouts are $(203,29)$, $(217,31)$, $(221,13)$ and $(221,17)$. The proof generator records Homebrew CaDiCaL 3.0.1 and `drat-trim` commit `2e3b2dc0ecf938addbd779d42877b6ed69d9a985`, and `src/final_check.py` checks the ledger counts and the three validation verdicts.

**`verify_base` agrees with Python** (`results/verify_base_crosscheck.log`). On each of the 7 stored bases, as stored, with one bit outside $F$ flipped, and with the coset moved by one, the C program and `safesat.violations` return identical $(h_1, h_2)$ in all 21 cases, 14 of which have violations. On the $k \ge 9$ words of §4 it reports none, and its two controls fail as they must: the $\mathbb{Z}_{226}$ base against the coset $r=1$ gives $h_2 = 112$, and the $k=9$ word with bit 1000 flipped gives $h_1 = 188$, $h_2 = 2$, both with exit 1 (`results/verify_base_k9_12.log`).

Known gap in coverage: `results/zm_k4.json` starts at $m=12$, so no $k=4$ periodic record is marked AP-free even though $\mathbb{Z}_{11}$ is 4-AP-free. This affects reporting only. $\mathbb{Z}_{11}$ is the $k=4$ base in `results/safebase_k4.json`, check H confirms that its word is cyclic 4-AP-free, and the $\delta_4$ bound comes from it.

The `.gitignore` rules for the compiled C binaries were checked with `git check-ignore -v`. All seven binaries in `src/` match a rule, and none of the new sources, the results, `refs/README.md` or `refs/fetch.sh` does.

## References

- [Er80] P. Erdős, *A survey of problems in combinatorial number theory*, Ann. Discrete Math. 6 (1980), 89–115. MR 593525.
- [LuPe12] L. Lu, X. Peng, *Monochromatic 4-term arithmetic progressions in 2-colorings of $\mathbb{Z}_n$*, [arXiv:1107.2888](https://arxiv.org/abs/1107.2888).
- [PRS08] P. Parrilo, A. Robertson, D. Saracino, *On the asymptotic minimum number of monochromatic 3-term arithmetic progressions*, JCTA 115 (2008).
- [CCS07] Cameron, Cilleruelo, Serra, *On monochromatic solutions of equations in groups*, Rev. Mat. Iberoam. 23 (2007).
- [Wo10] J. Wolf, *The minimum number of monochromatic 4-term progressions in $\mathbb{Z}_p$*, J. Comb. 1 (2010). **Not read — paywalled.**
- [BCG10] Butler, Costello, Graham. **Not read — paywalled.**
- [RL12] J. Rabung, M. Lotts, *Improving the use of cyclic zippers in finding lower bounds for van der Waerden numbers*, Electron. J. Combin. 19(2) (2012), #P35. [DOI 10.37236/2363](https://doi.org/10.37236/2363).
- [HHLM07] P. R. Herwig, M. J. H. Heule, P. M. van Lambalgen, H. van Maaren, *A new method to construct lower bounds for van der Waerden numbers*, Electron. J. Combin. 14 (2007), #R6. [DOI 10.37236/925](https://doi.org/10.37236/925).
- [SS78] R. S. Stevens, R. Shantaram, *Computer-generated van der Waerden partitions*, Math. Comp. 32 (1978), 635–636. [DOI 10.1090/S0025-5718-1978-0491468-X](https://doi.org/10.1090/S0025-5718-1978-0491468-X). **Not read**; $W(5;2)=178$ is taken from [RL12, Table 1].
- [KP08] M. Kouril, J. L. Paul, *The van der Waerden number $W(2,6)$ is 1132*, Exp. Math. 17 (2008), 53–61. [DOI 10.1080/10586458.2008.10129025](https://doi.org/10.1080/10586458.2008.10129025). **Not read**; the value is in its title and in [RL12, Table 1].

Fetch with `refs/fetch.sh`; nothing is redistributed here.
