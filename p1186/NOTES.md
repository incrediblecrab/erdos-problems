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

An evaluator that cannot recover a published value is not measuring the published quantity. Seven reproductions, all exact:

| source | published | recomputed here |
|---|---|---|
| [PRS08, Thm 5] | $c_3 \le 117/2192$ | $117/2192$ ✓ |
| [LuPe12, Thm 1] $B_{20}$, $\mathbb{F}_p$ blocks | $m_4(\mathbb{Z}_p) \le 17/150$ | $17/150$ ✓ |
| [LuPe12] periodic, $20\mid n$ | $m_4(\mathbb{Z}_n) \le 0.09$ | $m_4(\mathbb{Z}_{20}) = 9/100$ — exactly $0.09$ ✓ |
| [LuPe12] periodic, $22\mid n$ | $m_4(\mathbb{Z}_n) \le 0.086777$ | $m_4(\mathbb{Z}_{22}) = 21/242 = 0.0867769$ ✓ |
| [LuPe12] $B_{74}$ | $m_5(\mathbb{Z}_{74}) = 73/2738$ | $73/2738$ ✓ |
| [RL12, Table 1], two colours | $W(2,k) > 3703, 11495, 41265, 103474, 193941, 638727$ for $k = 7, \dots, 12$ | all six, as $(k-1)m+1$ for the largest good $m$ (`src/rabung.c`, §4) ✓ |
| [BGL14, Table 1 and Theorem 6], $r=2$ | $\delta_k \le 1/(2(m+1)(k-1))$ for $m = 11, 37, 139, 617, 1069$, $k = 4, \dots, 8$ | each QR colouring of $\mathbb{Z}_m$ avoids nontrivial $k$-APs with 0 coloured either way, giving $1/72$, $1/304$, $1/1400$, $1/7416$, $1/14980$ (gate check B; added November 2, 2026) ✓ |

The two threshold values are the tightest of these: $9/100$ is $0.09$ on the nose, and $21/242$ reproduces all six quoted digits of $0.086777$. (The best $m_4$ found here over moduli divisible by 22 is in fact $61/726 = 0.084022$ at $m=66$, better than their $\mathbb{Z}_{22}$ value — but this is an intermediate quantity, and after dividing by $2(k-1)=6$ it still does not reach $1/72$; see §6.)

The $B_{74}$ check is the sharp one, because it reproduces the *structure* and not just the value. Lu and Peng write:

> "All monochromatic 5-APs in $B_{74}$ are degenerated ones. Among them there are 74 5-APs with $d = 0$ and 72 5-APs with $d = 37$. This coloring gives $m_5(\mathbb{Z}_{74}) = 146/74^2 = 73/2738$."

`src/verify_zm.c` prints the $d$-breakdown and returns exactly 74 at $d=0$ and 72 at $d=37$.

Targets: $c_4 \le 1/72$ and $c_5 \le 1/304$ [LuPe12, eqs. (12)–(13)], and $\lim m_4(\mathbb{Z}_n)\le1/12$, $\lim m_5(\mathbb{Z}_n)\le1/38$ [LuPe12, Thm 5]. Over $\mathbb{F}_p$: $m_4(\mathbb{Z}_p)\le17/300+o(1)$ [Thm 1] and $m_5(\mathbb{Z}_n)\le3629/65712+o(1)$ for odd $n$ [Thm 4].

All quotations above were read from `refs/lupe12.pdf` directly (text extracted with `pypdf`), not from a secondary source.

## 3. The wildcard recursion

This is the ingredient behind every improvement here, and it is now a theorem (§3.1).

**Prior work.** The case $f = 1$, $F = \{0\}$ is [BGL14, Theorem 6], second case: Butler, Graham and Lu colour $\ell$ by the first nonzero base-$m$ digit of $\ell$, which is the chain below with $t = b^s$, and derive $1/(2(m+1)(k-1))$, including Lu–Peng's $1/72$ from $\mathbb{Z}_{11}$. The first case of the same theorem, $1/(2m(k-1))$ for any cyclic $k$-AP-free $\mathbb{Z}_m$, is the periodic transfer of [`PERIODIC-TRANSFER.md`](PERIODIC-TRANSFER.md). Their Table 1 gives $\delta_6 \le 1/1400$, $\delta_7 \le 1/7416$ and $\delta_8 \le 1/14980$ from quadratic-residue colourings of $\mathbb{Z}_{139}$, $\mathbb{Z}_{617}$ and $\mathbb{Z}_{1069}$. These notes did not cite [BGL14] before November 2, 2026, and called the $k \ge 6$ values unpublished; that was wrong. What remains this repository's own, as far as the literature searched shows, is the recursion with cosets of size $f > 1$ and its proof for every $f$, the base $B_{44}$, and zipped colourings used as bases. "Not found in the literature searched" is the claim, not "new".

Lu and Peng found their $k=4$ bound by noticing that the two halves of a good $\mathbb{Z}_{22}$ colouring differ in one bit, giving an 11-bit pattern with a wildcard. Running the same comparison one level up — on the optima found at $m=88$ and $m=132$ for $k=5$ — both reduce to the **same** 44-bit base with exactly four free positions $\{3,14,25,36\}$, the coset $x\equiv3\pmod{11}$.

**Bases.** A *base* is a 2-colouring $c$ of $\mathbb{Z}_b$ together with a proper coset $F = r + g\mathbb{Z}_b$ of the subgroup of order $f<b$, where $g = b/f$. It is *valid* for $k$ when

- **(H1)** every $k$-AP $a, a+d, \dots, a+(k-1)d$ of $\mathbb{Z}_b$ with $d \not\equiv 0 \pmod b$ and no term in $F$ is non-monochromatic, and
- **(H2)** every $k$-AP with $d \not\equiv 0 \pmod b$ and terms both inside and outside $F$ is non-constant on its terms outside $F$.

The colours on $F$ enter neither condition. (H2) is exactly the test `cosetprobe.py:safe_coset`: an AP is *dangerous* when its outside terms agree, because colouring the inside to match would make it monochromatic. Both conditions are finite checks over $\mathbb{Z}_b$.

**Theorem.** Let $(c, F)$ be a valid base with $f<b$, $t \ge 1$, and $h$ any 2-colouring of $\mathbb{Z}_{ft}$. Colour $\mathbb{Z}_{bt}$ by $C(x) = c(x \bmod b)$ when $x \not\equiv r \pmod g$, and $C(r + gu) = h(u)$ otherwise. Then

$$\#\mathrm{mono}_k(\mathbb{Z}_{bt}, C) \;=\; (b-f)\,t^2 + \#\mathrm{mono}_k(\mathbb{Z}_{ft}, h),$$

where $\#\mathrm{mono}$ counts pairs $(a,d)$ with $d=0$ included, so that $\Phi = \#\mathrm{mono}/m^2$. Hence

$$m_k(\mathbb{Z}_{bt}) \;\le\; \frac{(b-f) + f^2\, m_k(\mathbb{Z}_{ft})}{b^2}, \qquad\text{fixed point } \frac{1}{b+f}, \qquad \delta_k \le \frac{1}{2(k-1)(b+f)} .$$

Earlier versions of these notes called the recursion "verified but not proved in general" and its limits conjectural. §3.1 proves it for every $k$, $b$, $f$ and $t$; for $f = 1$ the statement and a proof are in [BGL14, Theorem 6].

**The formula reproduces both published constants exactly, and the bases behind them are quadratic-residue colourings (§4).** "vs random" is the bound divided by the random-colouring value $1/((k-1)2^k)$, which equals $2^{k-1}/(b+f)$.

| $k$ | base | $b$ | $f$ | limit $1/(b+f)$ | $\delta_k \le$ | vs random | status |
|---|---|---|---|---|---|---|---|
| 4 | Lu–Peng $B_{11}$ = QR colouring of $\mathbb{Z}_{11}$ | 11 | 1 | $1/12$ | $1/72$ | 0.667× | [LuPe12] eq. (12) ✓ |
| 5 | Lu–Peng $B_{37}$ = QR colouring of $\mathbb{Z}_{37}$ | 37 | 1 | $1/38$ | $1/304$ | 0.421× | [LuPe12] eq. (13) ✓ |
| 5 | $B_{44}$, found by search here | 44 | 4 | $1/48$ | $\mathbf{1/384}$ | 0.333× | beats $1/304$ by 20.83 %; not found in the literature searched; optimal in the family (§5) |
| 6 | QR colouring of $\mathbb{Z}_{139}$ [BGL14] | 139 | 1 | $1/140$ | $1/1400$ | 0.229× | [BGL14] Table 1 ✓ (gate check B) |
| 6 | $B_{86}$, found by search here | 86 | 2 | $1/88$ | $1/880$ | 0.364× | superseded |
| 6 | zipped QR colouring of $\mathbb{Z}_{226}$ | 226 | 2 | $1/228$ | $\mathbf{1/2280}$ | 0.140× | beats $1/1400$ [BGL14] by 38.60 % |
| 7 | QR colouring of $\mathbb{Z}_{617}$ | 617 | 1 | $1/618$ | $1/7416$ | 0.104× | $=$ [BGL14] Table 1, **not improved** |
| 8 | QR colouring of $\mathbb{Z}_{1069}$ [BGL14] | 1069 | 1 | $1/1070$ | $1/14980$ | 0.120× | [BGL14] Table 1 ✓ (gate check B) |
| 8 | zipped QR colouring of $\mathbb{Z}_{1642}$ | 1642 | 2 | $1/1644$ | $\mathbf{1/23016}$ | 0.078× | beats $1/14980$ [BGL14] by 34.91 % |

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

For a prime $p$, colour $x \in \mathbb{Z}_p^*$ by its quadratic character: $c(x) = 1$ iff $x$ is a non-residue. This is the colouring behind the cyclic van der Waerden certificates of [RL12], and its use for this problem is [BGL14]: their Theorem 7 states that such a colouring with $0$ coloured arbitrarily avoids nontrivial $k$-APs once the residues do, and their search over primes produced Table 1. Because $c(dz) = c(z) \oplus c(d)$ for $z \ne 0$, every AP $a + jd$ with $d \ne 0$ is the window $y, y+1, \dots, y+k-1$ (where $y = a/d$) rescaled by $d$, with its colours kept or all flipped. So $\mathbb{Z}_p$ has no monochromatic $k$-AP with $d \ne 0$ iff

- (i) no $k$ consecutive elements of $1, \dots, p-1$ share a colour, and
- (ii) no window of $k$ consecutive residues through 0 is constant on its nonzero elements.

[RL12, p. 5] uses the same multiplicativity to cut the check from quadratic to linear time.

**Criterion (i) is (H1), and criterion (ii) is (H2), for $F = \{0\}$.** Every prime that passes the test is therefore a valid base with $f = 1$, giving $\delta_k \le 1/(2(k-1)(p+1))$. The zipped colouring of $\mathbb{Z}_{2p}$ [HHLM07], read through $\mathbb{Z}_{2p} \cong \mathbb{Z}_p \times \mathbb{Z}_2$, is $z(x) = c(x \bmod p) \oplus [x \text{ even}]$ for $x \ne 0, p$, with $z(0) \ne z(p)$. Its criteria are (H1) and (H2) for $F = \{0, p\}$. The condition $z(0) \ne z(p)$ concerns the difference $p$; it matters for full AP-freeness but not for the base. Every good zip is therefore a base with $f = 2$, giving $\delta_k \le 1/(2(k-1)(2p+2))$. The argument is the rescaling above, applied to odd and even differences separately, and the header of `src/rabung.c` spells it out. Nothing downstream relies on it: the gate verifies (H1) and (H2) directly for $k \le 8$ (checks C and H), and `src/verify_base.c` does so for the larger words below.

**Lu–Peng's bases are quadratic-residue colourings.** [BGL14] made this observation for $\mathbb{Z}_{11}$: "the nonzero elements colored white, $\{1,3,4,5,9\}$, are precisely the nonzero quadratic residues of $\mathbb{Z}_{11}$". $B_{11} = (1,1,1,0,1,*,0,1,0,0,0)$ and the 37-bit $B_{37}$ of [LuPe12], read cyclically from the wildcard, are exactly the colourings above for $p = 11$ and $p = 37$ (1 = non-residue), with the wildcard at 0. Their Property 2, "No matter which bit-value the '∗' takes, $B_{37}$ contains no non-degenerated monochromatic 5-APs of $\mathbb{Z}_{37}$", is (H1)+(H2) at $F = \{0\}$. Their Lemma 7, the $k=5$ step of the recursion, comes with "The proof is omitted"; Lemma 6, the $k=4$ step, is proved there. §3.1 covers every $k$ and every $f$.

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

For $k = 4, 5, 6$, [RL12] list the exact value of $W$. For $k=5$ the QR family peaks at 37, and search does better (44). [BGL14, Table 1] lists, for $k = 4, \dots, 12$, $m = 11, 37, 139, 617, 1069, 3389, 11497, 17863, 58013$, all with $f = 1$. Against those, the rows above have a larger $m + f$ at $k = 6, 8, 9, 11, 12$ (the zipped words, $f = 2$), the same at $k = 4, 5, 7, 10$, and only $k \le 8$ is claimed. Rows with $k \ge 9$ are **not** in `results/claims.json` and are not gate-checked. Their words were checked by `src/verify_base` (both conditions, 0 violations) and by `src/verify_zm` (cyclic AP-free), with output in `results/verify_base_k9_12.log`. These are fixed-size constructions. The ratio to random falls from 0.667× at $k=4$ to 0.035× at $k=12$, though not monotonically: it rises from 0.045× at $k=10$ to 0.053× at $k=11$. That says nothing about asymptotics, and no asymptotic formula for $\delta_k$ is claimed.

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

**Free breakpoints on $\{1,\dots,n\}$ do not help either (November 2, 2026).** Freeing the colour changes from the grid is what improved the $\mathbb{F}_p$ bounds (§11). `src/twofree.py` tries the same idea on $\{1,\dots,n\}$ in two families. One is block colourings with free boundaries, the family of [PRS08]'s word, with $B = 2, \dots, 30$ blocks. The other is two-scale colourings $i \mapsto C[i \bmod m][\text{coarse block of } i/n]$ with free coarse boundaries, which `src/twoscale.py` fixes on an equal grid. Each start alternates continuous optimisation of the boundaries with single-bit flips of $C$, then tries random kicks. The evaluator is the floating-point C port `src/twof.c`. It agrees with `exact.psi_two_scale` on 24 stored two-scale words and with `blocks.psi_blocks` on 30 random block colourings to $7.1\times10^{-17}$. Nothing below is certified, and nothing was stored, because nothing reached its target.

| $k$ | family | starts | best $\Psi_k$ (float) | target | gap | log |
|---|---|---|---|---|---|---|
| 3 | blocks, $B = 2..30$, free boundaries | 290 | 0.053375912409 | $117/2192 = 0.053375912409$ | 0 (67 starts end on [PRS08]'s word) | `results/twofree_k3_m1.log` |
| 3 | two-scale, $m = 2, 3, 4$, $B = 2..12$ | 190 | 0.053571428571 | $117/2192$ | $+0.37$ % | `results/twofree_k3.log` |
| 4 | two-scale, $m = 11, 22, 33, 44, 66$, $B = 2..8$ | 153 | 0.013927774461 ($m=66$, $B=8$) | $1/72 = 0.013888888889$ | $+0.28$ % | `results/twofree_k4.log` |

For $k=3$, 67 of the 290 block starts converge to [PRS08]'s twelve lengths over 548, and none goes lower. This is consistent with their word being a local optimum among free-boundary block colourings, and with the result of Greenwood, Kariv and Williams [GKW24] that it is optimal among antisymmetric colourings with at most 12 blocks; it is evidence, not a proof. The two-scale runs did not even recover it from random colour tables. For $k=4$ the free coarse boundaries do lower the two-scale family, from the equal-grid $2917/209088 = 0.0139511$ ($m=66$, $B=4$) to 0.0139278 in floating point at $B = 8$. Its best value was still falling as $B$ grew, but it stays 0.28 % above the recursion limit $1/72$. A plausible reason, not a proof: the recursion reaches $1/72$ by recolouring the free coset at every scale, so the progressions with difference $\equiv 0 \pmod{11}$ are handled recursively. The two-scale family handles them with one coarse block colouring per residue, and block colourings are poor at $k=4$: the best stored one is $411/23104 = 0.01779$, 28 % above $1/72$.

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
- **$k \ge 6$ over $\{1,\dots,n\}$: corrected November 2, 2026.** Earlier versions said no published bound was located. [BGL14, Table 1 with Theorem 6] gives $\delta_6 \le 1/1400$, $\delta_7 \le 1/7416$ and $\delta_8 \le 1/14980$, so $1/7416$ is their value, and $1/2280$ and $1/23016$ are compared with theirs (§3). The source was found by the coordinating agent and its Table 1 and Theorems 6 and 7 were re-read here from the arXiv source and `refs/bgl14.pdf`; the gate reproduces the five $k \le 8$ entries (check B). Over $\mathbb{F}_p$, no published bound for $k \ge 6$ was located; [BGL14] treats $\{1,\dots,n\}$ only.
- **Records, no claim.** A proof claim by Carlos Toledo posted on the erdosproblems.com forum on October 5, 2026 ([thread](https://www.erdosproblems.com/forum/thread/1186/proof-claims), unrefereed, written with Claude) asserts $\delta_3 = 117/2192$, i.e. that [PRS08]'s colouring is optimal. It was not checked here. [formal-conjectures](https://github.com/google-deepmind/formal-conjectures/blob/main/FormalConjectures/ErdosProblems/1186.lean) added Erdős 1186 on October 20, 2026 (PR #7467), with `erdos_1186.variants.prs_conjecture`, $\delta_3 = 117/2192$, tagged open. Both were read on November 2, 2026.
- **Superseded, kept for the record (it was right to hedge):** "**$k\ge6$: no published bound was located, for either problem.** The $k \ge 6$ entries are therefore reported as "no published bound found", not as records. That the sources read contain nothing is not proof that nothing exists. The recursion theorem, the quadratic-residue reading of Lu–Peng's bases, and the values for $\delta_6$, $\delta_7$, $\delta_8$ may already be known to specialists; [BCG10] in particular was not read."
- **The optimality statements of §5 rest on computer results; the bounds do not.** The scans cover $b \le 150$, $200$ and $226$ directly, and larger $b$ are excluded by the cap. The cap uses $W(5;2)=178$ [SS78] and $W(6;2)=1132$ [KP08], neither of which was read. Both values are taken from [RL12, Table 1], and [KP08]'s from its title as well.
- **Every UNSAT verdict behind §5 is now externally proof-checked.** `src/prove_unsat.py` rebuilds each CNF with `safesat.encode_cnf`, runs Homebrew CaDiCaL 3.0.1 with `--lrat`, and checks the proof with `lrat-check` from `drat-trim` commit `2e3b2dc0ecf938addbd779d42877b6ed69d9a985`. `results/lrat_check.json` records 1122 accepted LRAT proofs and no rejections or timeouts: 89/89 UNSAT instances for $k=3$, 423/423 for $k=4$, 538/538 for $k=5$, and 72/72 for $k=6$. The last four $k=6$ instances, $(203,29)$, $(217,31)$, $(221,13)$ and $(221,17)$, timed out under the first run's 300 s per-instance cap. On November 1, 2026 they were re-run with `--only`, a 10800 s cap, two at a time, and `--no-plain`, which keeps CaDiCaL's default inprocessing; the other 1118 proofs used `--plain`, and each record stores its `solver_args`. They also used `--stream`, which pipes the proof through a FIFO into `lrat-check` as it is written, so no proof file was stored. The four proofs were 2.3 to 4.1 GB of text, and the solver and checker ran together for 1669 to 7104 s of wall time on a shared, heavily loaded machine. Those times include the checker's back-pressure, so they are not comparable with the solve times above. Free disk space was below the 20 GiB floor set for this run when it started, because of other workloads. Since streaming writes only the CNF, the run went ahead with a 2 GiB kill threshold instead. Each record's variable and clause counts equal those of the timed-out records it replaced, so the CNFs were rebuilt identically. The two $b=221$ encodings include lex-leader symmetry-breaking clauses; LRAT checks those encoded CNFs, but soundness of adding those clauses rests on the symmetry argument above rather than on LRAT. FORCED answers need no solver. No bound depends on an UNSAT answer: every bound in `results/claims.json` comes from an explicit word that the gate re-checks.

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

`src/final_check.py` re-derives every claim from the stored artifacts and asserts the contents of `results/claims.json` against them; it does not read any cached score. The full run, `python src/final_check.py --threads 4`, passes all 155 checks with exit 0 in 260 s (`results/final_check_full.log`, November 2, 2026, on a shared machine with load average about 20 on 11 cores). The November 1 run passed 139 in 175 s; the new checks are those of §11.1, check C′ and the [BGL14] reproduction, and `verify_zpx` on all 60 residues of the two $m \approx 970$ words took 56 and 59 s when timed on its own (§11.1).

**The gate was tested by planting defects, not merely observed to pass.** A checker that has only ever seen good input is untested. Each defect below was planted in a stored file, the gate was run as `final_check.py --fast --no-brute`, and the file was put back. The clean run and the run after restoring both pass all 100 checks with exit 0, and the three files touched were byte-identical to the originals afterwards.

| planted defect | exit | FAIL lines |
|---|---|---|
| bit 1 of the $\mathbb{Z}_{226}$ base flipped (outside $F$) | 1 | 5: re-certification (`H1/H2 or shape`), `safe_coset` (6 dangerous progressions), AP-freeness (284 monochromatic pairs), and the identity at $t=2$ ($1128 \ne 224\cdot2^2+16$) and at $t=3$ ($2532 \ne 224\cdot3^2+12$) |
| bit 0 of the same base flipped (inside $F=\{0,113\}$) | 1 | 1: AP-freeness (228 monochromatic pairs) |
| the same base's coset moved to $r=1$ | 1 | 4: re-certification, `safe_coset` (112 dangerous progressions), and the identity at $t=2$ ($1136 \ne 224\cdot2^2+16$) and at $t=3$ ($2532 \ne 224\cdot3^2+12$) |
| $\delta_6$ claim changed to $1/2300$ | 1 | `[FAIL] delta_6: certified 1/2280 <= claimed 1/2300` |
| $\delta_5$ claim changed to $1/400$ | 1 | `[FAIL] delta_5: certified 1/384 <= claimed 1/400` |
| one bit of the stored $m=484$ word in `results/wildcard_k5_deep.json` flipped, its recorded `phi` left intact | 1 | `[FAIL] zm   k=5: 12 words re-certify   1 bad: [('484', '633/234256 != 111/42592')]` |
| bit 100 of the $m=496$ $\mathbb{F}_p$ word in `results/zp_k4_c.json` flipped, its recorded `value` left intact | 1 | 3: re-certification (`40499/369024 != 20237/184512`), `circle.py` agreement, and `verify_zpx` (check D′); check D still passed, at $|\Delta\Phi| = 6.79\times10^{-5}$ inside its $2.00\times10^{-4}$ tolerance |
| $\tilde\delta_4$ claim changed to $20237/369025$ | 1 | `[FAIL] delta~_4: certified 20237/369024 <= claimed 20237/369025` |
| bit 100 of the $m=974$ $k=6$ word in `results/zp_k6_c.json` flipped, its recorded `value` left intact | 1 | 4: re-certification by `exact.phi_bits` (`43609/1897352 != 81632/3557535`), check C′, `circle.py` agreement, and `verify_zpx` on 60 of 60 residues; check D still passed |
| bit 1000 of the $m=2998$ $k=4$ word in `results/zp_k4_c.json` flipped, its recorded `value` left intact | 1 | 4: re-certification by `exact.phi_bits` (`246457/2247001 != 1478233/13482006`), check C′, `circle.py` agreement, and `verify_zpx`; check D still passed |
| $\tilde\delta_6$ claim changed to $40816/3557536$ | 1 | `[FAIL] delta~_6: certified 40816/3557535 <= claimed 2551/222346` |
| [BGL14] Table 1 entry for $k=6$ changed from 139 to 138 in `src/final_check.py` | 1 | 2: the QR re-test (`Z_138`, 2602 and 2688 monochromatic pairs) and the Table 1 quotation (`missing or inconsistent: [6]`) |
| $\delta_6$'s quoted improvement changed from 38.60 to 38.61 in `results/claims.json` | 1 | `[FAIL] delta_6: quoted improvement 38.61 % equals the exact 38.60 %` |

Two rows need comment. Colours on $F$ enter neither (H1) nor (H2), so the bit-0 defect is invisible to every base check except AP-freeness. That check was added when a reading of the gate showed it was missing, and this defect was planted to test it. The $m=484$ row shows that $\Phi$ is recomputed from the colouring rather than read from the file, so a corrupted certificate cannot pass by carrying a correct-looking number.

The last two rows were planted later, for the $\mathbb{F}_p$ words of §11, and the gate was run as `final_check.py --fast`, with brute force, because check D′ is part of it. Both files were byte-identical to the originals afterwards. These two runs were made while the LRAT ledger was being extended, so each also printed the then-stale LRAT count failure, which is unrelated to the plant and is not listed. The bit-100 row is the reason check D′ exists: a one-bit change of a 496-block word moves $\Phi$ by less than check D's tolerance, so only an exact comparison catches it.

The last five rows were planted on November 2, 2026, for §11.1 and the [BGL14] comparison (`results/planted_defects_2026-11-02.log`). The two bit flips ran as `final_check.py --fast --threads 4`, with brute force; the other three ran with `--no-brute` as well. Every planted file was restored from a byte copy, and its SHA-256 matched the value taken before the plants.

**The SAT tooling was tested the same way** (`results/safesat_selftest.log`). `safesat.py --self-test` decides all 588 instances $(k,b,f)$ with $k = 3, 4, 5$ and $b < 30, 60, 90$ with and without symmetry breaking, and every answer agrees; 24 are SAT. Lex-leader constraints for random non-symmetries turn 13 of the 24 UNSAT, so the comparison can see a changed answer. Replacing the symmetry group by random permutations drops agreement to 578 of 588, and the self-test exits 1.

**External LRAT proof checking now covers every UNSAT scan answer** (`results/lrat_check.json`). The checker path was validated on $(k,b,f)=(4,12,2)$; the stored validation record is from the November 1, 2026 run and used the `--no-plain` flags, and the gate reads its three verdicts, not its flags: `lrat-check` accepts the correct CaDiCaL LRAT proof, rejects the proof after corrupting a hint/clause, and rejects the same proof against a CNF with one clause dropped. The checked corpus is 89 $k=3$ UNSAT instances, 423 $k=4$ instances, 538 $k=5$ instances, and 72 $k=6$ instances. The proof generator records Homebrew CaDiCaL 3.0.1 and `drat-trim` commit `2e3b2dc0ecf938addbd779d42877b6ed69d9a985`. `src/final_check.py` checks the ledger counts, the three validation verdicts, and that every UNSAT record in `safesat_k{3..6}.json` has an accepted proof under its own key. That last check was tested on a scratch copy of `results/` with one extra UNSAT record planted in `safesat_k6.json`: the counts still matched, and the check failed, naming the planted key. The streaming path was tested before use: on two $k=5$ instances its byte and line counts equal those of the stored file-mode proofs, a wrapper that corrupts the proof before streaming it is rejected, and a 20 s cap is recorded as a timeout with no proof left on disk.

**`verify_base` agrees with Python** (`results/verify_base_crosscheck.log`). On each of the 7 stored bases, as stored, with one bit outside $F$ flipped, and with the coset moved by one, the C program and `safesat.violations` return identical $(h_1, h_2)$ in all 21 cases, 14 of which have violations. On the $k \ge 9$ words of §4 it reports none, and its two controls fail as they must: the $\mathbb{Z}_{226}$ base against the coset $r=1$ gives $h_2 = 112$, and the $k=9$ word with bit 1000 flipped gives $h_1 = 188$, $h_2 = 2$, both with exit 1 (`results/verify_base_k9_12.log`).

Known gap in coverage: `results/zm_k4.json` starts at $m=12$, so no $k=4$ periodic record is marked AP-free even though $\mathbb{Z}_{11}$ is 4-AP-free. This affects reporting only. $\mathbb{Z}_{11}$ is the $k=4$ base in `results/safebase_k4.json`, check H confirms that its word is cyclic 4-AP-free, and the $\delta_4$ bound comes from it.

The `.gitignore` rules for the compiled C binaries were checked with `git check-ignore -v`. All eight binaries in `src/`, `verify_zpx` included, match a rule, and none of the sources, the results, `refs/README.md` or `refs/fetch.sh` does (rechecked November 1, 2026).

## 11. $\mathbb{F}_p$: colour changes off the equal-block grid

The $\tilde\delta_k$ bounds use the block colouring $x \mapsto c[\lfloor mx/p \rfloor]$ of $\mathbb{F}_p$. Its limit is $\Phi_k(S)$ for the set $S \subset \mathbb{R}/\mathbb{Z}$ that $c$ marks, with the normalisation of §1: ordered pairs $(x,d)$, difference zero included (a null set in the limit), and $\tilde\delta_k \le \Phi_k(S)/2$. A colouring of the circle that changes colour at arbitrary rational points of common denominator $M$ is an equal-block word of length $M$, so freeing the breakpoints adds nothing to the family in principle. What it changes is the search. The stored equal-block searches cover 118 values of $m$ up to 153 for $k=4$ and 80 values up to 114 for $k=5$, and at each $m$ every colour change sits on a multiple of $1/m$.

**Method.** `src/circle.py` computes $\Phi_k$ directly from the breakpoints $0 = a_0 < \dots < a_B = 1$ by clipping the unit square of $(x,d)$ with the strips $a_i + w \le x + jd < a_{i+1} + w$, $0 \le w \le j$, kept only for one colour at a time. It uses `blocks.py`'s clipper and none of `exact.py`'s carry regions. `src/circleopt.py` starts from the best stored words, minimises $\Phi_k$ over the arc lengths in floating point (L-BFGS-B, then Powell, then 8 random perturbations of the incumbent per start, from 4 starts), rounds the breakpoints onto every grid $100 \le M \le 600$, polishes the 15 best grids by moving single cuts by one block, and re-derives the best grid word in exact arithmetic with `exact.phi`. Transcripts, including the floating-point optima, are in `results/circleopt_k4.log` and `results/circleopt_k5.log`. Those floating-point values are search output, not certificates.

**Result.** Started from the stored $m=152$ word for $k=4$ and $m=114$ word for $k=5$, the continuous search moves to breakpoints with a lower floating-point $\Phi_k$ (first line of each transcript). That the stored words are not local minima in breakpoint space rests on this search leaving them, not on a derivative test. The rounded words give

| | $m$ | $\Phi_k$ (exact) | $\tilde\delta_k \le$ | file |
|---|---|---|---|---|
| $k=4$ | 496 | $20237/184512$ | $20237/369024$ | `results/zp_k4_c.json` |
| $k=5$ | 600 | $20741/432000$ | $20741/864000$ | `results/zp_k5_c.json` |

against $159/2888$ and $1283/51984$ from the equal-block grids. The $k=4$ word, like the $m=152$ word it came from, is antipodal: $c(x + 1/2) = 1 - c(x)$. The $k=5$ word is not. It came from the $m=112$ start and has 30 arcs where the $m=114$ word has 26. The four starts ended at three different local optima for each of $k=4$ and $k=5$, so the search is local and better colourings of this kind may exist. No optimality is claimed. §11.1 widens the search and adds $k = 6, 7$; the words in this table remain stored and certified, but are no longer the best.

**Certification.** Each word is certified the way the earlier $\mathbb{F}_p$ words are, plus two new checks. `exact.phi` is the certificate (check C). `src/circle.py` must give the same rational from the word's breakpoints, a different decomposition (check F). `src/verify_zp.c` must measure it within $12/p$ at $p = 60013$ (check D). The new `src/verify_zpx.c` recovers $\Phi_k$ *exactly* from brute-force counts, sharing no code with either evaluator (check D′). It colours $\mathbb{Z}_{Mt}$ with blocks of size $t$ and counts the monochromatic pairs $(x,d)$, $d=0$ included. That count is the number of lattice points of $tP$ for a half-open polygon $P \subset [0,M)^2$ cut out by the sides of the square and the lines $x + jd \in \mathbb{Z}$, whose vertices have denominators dividing some $j - j' < k$. By Ehrhart's theorem it is therefore a quasi-polynomial of degree 2 in $t$ with period dividing $L = \mathrm{lcm}(1,\dots,k-1)$ and leading coefficient $\mathrm{area}(P) = M^2\Phi_k$. The second difference at step $L$ thus gives $\Phi_k$ exactly. The program tests the structure rather than assuming it: for every residue $t_0 \le L$, the second difference must be the same and the third must vanish. Until November 2, 2026 the gate ran only two of the 60 residues for $k = 6, 7$; it now runs all of them (§11.1). It reproduces Lu–Peng's $17/150$ and every stored best $\mathbb{F}_p$ word for $k = 4, \dots, 7$, and it must equal `exact.phi` digit for digit. That matters here because check D's tolerance at $p=60013$ is of the same order as the improvement over the $m=152$ word.

The improvement is over this repository's previous equal-block values. Against the published bounds the gate (check E) reported $\tilde\delta_4 \le 20237/369024$ as 3.22 % below $17/300$ and $\tilde\delta_5 \le 20741/864000$ as 13.06 % below $3629/131424$; §11.1 has the current values, 3.25 % and 13.19 %. These are explicit colourings, so they are upper bounds only and do not determine $\tilde\delta_4$ or $\tilde\delta_5$.

### 11.1 November 2, 2026: a wider search, and $k = 6, 7$

**Search.** `src/circleopt.py --wide` uses `src/circlef.c`, a C port of `circle.phi_circle` in doubles, for search only. At the start of every run it is compared with the Python evaluator on 200 random colourings; the largest difference in the four runs was $2.8\times10^{-15}$. At about 30 arcs it takes 0.09 to 0.25 ms per evaluation for $k = 4, \dots, 7$, against 2.7 to 7.2 ms for the Python evaluator, about 30 times faster. A run has three stages. Stage A optimises locally from the 6 best stored words and from random Dirichlet arc lengths for every even number of arcs $B$ from 2 to 30 ($k = 4, 5$; 10 starts per $B$ and symmetry class, 466 starts in all) or to 36 ($k = 6, 7$; 6 per class, 330 starts). The random starts are unconstrained or confined to one of three symmetric subspaces: antipodal, $c(x + 1/2) = 1 - c(x)$; mirror, $c(-x) = c(x)$; or both. Stage B basin-hops the 12 ($k=4,5$) or 10 ($k=6,7$) best distinct stage-A optima without symmetry constraints, for 80 or 60 hops. Each hop is a log-normal perturbation, a split of one arc into three, or a merge of three arcs into one, followed by a full local optimisation, and only improvements are kept. Finally the three best basins are rounded onto every grid $100 \le M \le 3000$ ($k = 4, 5$) or $\le 1000$ ($k = 6, 7$) and polished. The polish moves single cuts by one or two blocks and shifts single arcs by one block. Transcripts are `results/circleopt_k{4,5,6,7}_wide.log`, with every start's seed and result. Their floating-point values are search output.

| | previous best | best continuous $\Phi_k$ (float) | arcs | $m$ | $\Phi_k$ (exact) | $\tilde\delta_k \le$ | decimal | change |
|---|---|---|---|---|---|---|---|---|
| $k=4$ | $20237/184512$, $m=496$ (§11) | 0.109644025537 | 26 | 2998 | $1478233/13482006$ | $1478233/26964012$ | 0.05482244 | $-0.03$ % |
| $k=5$ | $20741/432000$, $m=600$ (§11) | 0.047942290967 | 30 | 2969 | $1267880/26444883$ | $633940/26444883$ | 0.02397212 | $-0.14$ % |
| $k=6$ | $2821/108000$, $m=60$ (grid) | 0.022901316461 | 46 | 974 | $81632/3557535$ | $40816/3557535$ | 0.01147311 | $-12.15$ % |
| $k=7$ | $1517/115320$, $m=62$ (grid) | 0.010805738918 | 36 | 973 | $615023/56803740$ | $615023/113607480$ | 0.00541358 | $-17.69$ % |

The words are in `results/zp_k{4,5,6,7}_c.json`. Against the published bounds, $\tilde\delta_4$ is 3.25 % below $17/300$ and $\tilde\delta_5$ is 13.19 % below $3629/131424$; no published bound for $k \ge 6$ over $\mathbb{F}_p$ was found.

**What did not improve.** For $k = 4$ and $5$ the wider search found no better basin. The best optimum reached from any random start was 0.1097950 for $k=4$ (22 arcs, antipodal) and 0.0485199 for $k=5$, both above the basins reached from the stored words, and no hop left those basins except, for $k = 4$, to a point 1.0e-6 lower in the same 26-arc antipodal basin. The $k=4, 5$ gains are therefore almost entirely the finer grids: the rounded words sit $8.6\times10^{-7}$ and $1.96\times10^{-6}$ above the continuous values, against $3.3\times10^{-5}$ and $6.9\times10^{-5}$ at $m = 496$ and $600$.

**$k = 6, 7$ are not converged.** Both best basins started from the stored $m = 62$ word, and most improving hops were splits; at $k = 6$ the best basin went from 42 to 46 arcs in three improving hops, two of them splits. More arcs, more hops, or grids beyond $m = 1000$ (where rounding still costs $4.5\times10^{-5}$ and $2.1\times10^{-5}$) would very likely lower these values. They are what one finite run found.

**Certification.** `exact.phi` tabulates about $m^2$ times the number of carry cells (8, 14, 24, 36 for $k = 4, \dots, 7$), and its memory grows accordingly. It was measured at 2.9 GB and 19 s for $m = 1200$, $k = 4$; 3.8 GB and 29 s for $m = 1000$, $k = 5$; and 4.7 GB and 28 s for $m = 1500$, $k = 4$ (`results/exact_phi_oneoff.log`). At 180 to 270 bytes per cell, the $m = 2998$ and $2969$ words would need about 18 and 31 GB, more than the machine has. So `exact.py` gained `phi_bits`. It is the same carry-region formula with the same exact areas, summed as integer popcounts of ANDed rotations of the word instead of through the table, and it takes 0.01 s where `exact.phi` takes 20 to 30 s. The gate certifies with `exact.phi` when $m^2$ times the cell count is at most `exact.PHI_TABLE_CELLS` $= 10^7$ ($m \le 1118, 845, 645, 527$ for $k = 4, \dots, 7$), and with `exact.phi_bits` above it (check C). `phi_bits` must equal `exact.phi` on every 0/1 word of length below 10 for $k = 3, \dots, 7$ (5110 pairs, check A), and must reproduce every stored $\mathbb{F}_p$ value (check C′). All four new words are above the cap, so in the gate their certificate is `phi_bits`. In addition:

- `circle.phi_circle` gives the same rational from the breakpoints (check F);
- `verify_zpx` gives the same rational by brute-force counting on every residue mod $L$, including all 60 for $k = 6, 7$, in 0.5, 3.4, 56 and 59 s with 4 threads (check D′);
- `verify_zp` agrees within $12/p$ at $p = 60013$ (check D);
- for $k = 6, 7$, `exact.phi` with its table was run once by hand and gives the same value, in 81 s at 5.5 GB and 145 s at 6.2 GB (`results/exact_phi_oneoff.log`). For $k = 4, 5$ this was not possible.

**`verify_zpx` was rewritten for speed.** It now counts with 64-bit bitsets and uses the reflection $(x, d) \mapsto (x + (k-1)d, -d)$ to scan only $d \le N/2$. At $t = 1, 2$ it also runs the old plain double loop over every $(x, d)$ and exits 1 on any difference. On the stored words compared, it agrees with the committed plain-loop program (`results/verify_zpx_selftest.log`), and two planted defects make it exit 1. One counts $d = N/2$ twice: its interpolated value was still right, and only the plain-loop comparison caught it. On an antipodal word the same plant is invisible, because the count at $d = N/2$ is zero there.

**A second exact evaluator for arbitrary breakpoints.** `circle.phi_sweep` integrates over $d$. For fixed $d$ the set of $x$ whose $k$ terms share a colour is a union of intervals with endpoints $a_i - jd$. Their cyclic order changes only at $d = (a_{i'} - a_i + w)/q$ for $q < k$, and between such points the measure is linear in $d$, so the trapezoid rule over them is exact. It shares no geometry with `phi_circle` and needs no common grid. It agrees with `phi_circle` on 60 random colourings, half of them with unrelated denominators near $10^6$, and on the $k = 5$, $m = 600$ word.

**The continuum optimum is rational but is not certified here.** $\Phi_k$ is piecewise quadratic in the breakpoints. The polygons it is made of have vertices that are rational-linear in the breakpoints, and their combinatorial type is constant on the cells of a polyhedral subdivision. On one cell $\Phi_k$ is therefore an exact quadratic with rational coefficients, and its stationary point is rational. `src/circlestat.py` relaxes the stored word to the floating-point minimum and snaps it to a base point of denominator $2^{30}$. There it takes the exact gradient and Hessian by central differences, with `phi_circle` in Fractions at steps $2^{-34}$ and $2^{-35}$. The two must be identical, which tests that the stencil lies on one quadratic piece. It then solves for the stationary point $a^*$ exactly and evaluates $\Phi_k(a^*)$ with both `phi_circle` and `phi_sweep`. Logs are `results/circlestat_k4.log` and `results/circlestat_k5.log`.

| | arcs | derivatives at the two steps | Hessian | $\Phi_k(a^*)$ = quadratic's prediction | `phi_circle` = `phi_sweep` | denominator of $a^*$ | $\Phi_k(a^*)$ | below the certified word by |
|---|---|---|---|---|---|---|---|---|
| $k=4$ | 26 | identical | positive definite | yes | yes | 17 digits | $322879924836886/2944801809830641 = 0.109644025536461$ | $8.60\times10^{-7}$ |
| $k=5$ | 30 | identical | positive definite | yes | yes | 64 digits | 0.047942290966596 (64-digit denominator) | $1.96\times10^{-6}$ |

So, within the cell found, the continuum local optimum is a single colouring with rational breakpoints, and $\Phi_k(a^*)/2$ is an upper bound on $\tilde\delta_k$ whose exact value two independent evaluators agree on. It is **not claimed**. As an equal-block word it has $4.7\times10^{16}$ or about $10^{64}$ blocks, so neither `exact.phi`, `phi_bits` nor the brute-force `verify_zpx` can evaluate it, and two decomposition-based evaluators fall short of the gate's standard. That $a^*$ lies in the cell is inferred from the exact match with the quadratic's prediction, not proved. The certified words are within $2\times10^{-6}$ of these values anyway, so little is at stake: rounding, not the continuum, is the last step.

## References

- [Er80] P. Erdős, *A survey of problems in combinatorial number theory*, Ann. Discrete Math. 6 (1980), 89–115. MR 593525.
- [LuPe12] L. Lu, X. Peng, *Monochromatic 4-term arithmetic progressions in 2-colorings of $\mathbb{Z}_n$*, [arXiv:1107.2888](https://arxiv.org/abs/1107.2888).
- [PRS08] P. Parrilo, A. Robertson, D. Saracino, *On the asymptotic minimum number of monochromatic 3-term arithmetic progressions*, JCTA 115 (2008).
- [CCS07] Cameron, Cilleruelo, Serra, *On monochromatic solutions of equations in groups*, Rev. Mat. Iberoam. 23 (2007).
- [Wo10] J. Wolf, *The minimum number of monochromatic 4-term progressions in $\mathbb{Z}_p$*, J. Comb. 1 (2010). **Not read — paywalled.**
- [BCG10] Butler, Costello, Graham. **Not read — paywalled.**
- [GKW24] Greenwood, Kariv, Williams, *From discrete to continuous: monochromatic 3-term arithmetic progressions*, Math. Comp. (2024). [arXiv:2301.00336](https://arxiv.org/abs/2301.00336), [DOI 10.1090/mcom/3970](https://doi.org/10.1090/mcom/3970).
- [BGL14] S. Butler, R. Graham, L. Lu, *Unrolling residues to avoid progressions*, Math. Mag. 87(2) (2014), 83–94. [arXiv:1209.2687](https://arxiv.org/abs/1209.2687), [DOI 10.4169/math.mag.87.2.83](https://doi.org/10.4169/math.mag.87.2.83). Theorem and table numbers are those of arXiv:1209.2687v1, which was read; the journal version was not read, and only its metadata (title, authors, volume, pages) was checked, against Crossref.
- [RL12] J. Rabung, M. Lotts, *Improving the use of cyclic zippers in finding lower bounds for van der Waerden numbers*, Electron. J. Combin. 19(2) (2012), #P35. [DOI 10.37236/2363](https://doi.org/10.37236/2363).
- [HHLM07] P. R. Herwig, M. J. H. Heule, P. M. van Lambalgen, H. van Maaren, *A new method to construct lower bounds for van der Waerden numbers*, Electron. J. Combin. 14 (2007), #R6. [DOI 10.37236/925](https://doi.org/10.37236/925).
- [SS78] R. S. Stevens, R. Shantaram, *Computer-generated van der Waerden partitions*, Math. Comp. 32 (1978), 635–636. [DOI 10.1090/S0025-5718-1978-0491468-X](https://doi.org/10.1090/S0025-5718-1978-0491468-X). **Not read**; $W(5;2)=178$ is taken from [RL12, Table 1].
- [KP08] M. Kouril, J. L. Paul, *The van der Waerden number $W(2,6)$ is 1132*, Exp. Math. 17 (2008), 53–61. [DOI 10.1080/10586458.2008.10129025](https://doi.org/10.1080/10586458.2008.10129025). **Not read**; the value is in its title and in [RL12, Table 1].

Fetch with `refs/fetch.sh`; nothing is redistributed here.
