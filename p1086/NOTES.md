# Erdős #1086: working notes

The initial computations here were done on September 23, 2026, on an Apple M3 Pro (11 cores) under macOS 27.0, with Apple clang 21.0.0, Python 3.14.7, sympy 1.14.0 and numpy 2.5.3. The $g(7)\le20$ check was added on September 28, 2026. The exact values $g(7)=18$ and $g(8)=28$ were added on October 19, 2026, computed in a 4-core x86-64 Linux container with gcc 13.3.0, Python 3.13.16 and sympy 1.14.0, where the grid histograms were also regenerated and matched the stored SHA-256 hashes. The exact value $g(9)=36$ and the searches of §5 that go beyond the small boxes were added on November 1, 2026, on the Apple M3 Pro under macOS 27.2, with Apple clang 21.0.0, Python 3.14.7, sympy 1.14.0 and z3 5.1.0. The exact value $g(10)=48$ and the constructions of §5 were added on November 2, 2026, on the same machine with Python 3.15.2, sympy 1.14.0, python-sat 1.9.dev17, networkx 3.7, kissat 4.0.4, cadical 3.0.1 and drat-trim's `lrat-check` at commit 2e3b2dc0. Sources are cited by the tags in [`refs/README.md`](refs/README.md).

**Notation.** The *doubled area* of a triangle $pqr$ is $|\det(q-p,\,r-p)|$, twice its area; for lattice points it is a positive integer. $G_m=\{0,\dots,m-1\}^2$ is the $m\times m$ grid, $n=m^2$, and $T_N(m)$ is the number of triangles of $G_m$ with doubled area $N$. The code and `results/claims.json` write the doubled area as $D$. $\sigma(N)$ is the sum and $d(N)$ the number of divisors of $N$; $\gamma$ is Euler's constant; $\log$ is natural.

## 1. The record

The problem page, fetched September 23, 2026 (it reports "last edited 16 October 2025"), states: "Let $g(n)$ be minimal such that any set of $n$ points in $\mathbb{R}^2$ contains the vertices of at most $g(n)$ many triangles with the same area. Estimate $g(n)$." It is labelled **open**, "This is open, and cannot be resolved with a finite computation", with no comments, no proof claims and "Formalised statement? No". The google-deepmind/formal-conjectures repository had no file for 1086 when checked the same day; the commit was not recorded. It has since added `FormalConjectures/ErdosProblems/1086.lean` (commit a7b87c8d, PR #7483, October 22, 2026; read on November 2, 2026 at HEAD 7bdc4a42). There `triangleAreaCount P a` is the number of 3-element subsets $\{p,q,r\}$ of a finite set $P\subset\mathbb R^2$ with `|EuclideanGeometry.triangle_area p q r| = a`, and `g n` is `sSup {m | ∃ P a, P.card = n ∧ 0 < a ∧ triangleAreaCount P a = m}`. That is the reading used here: only positive areas count (item 1 below), and $g(n)$ is the largest such count over $n$-point sets and areas. The file states the problem and the [ErPu71] lower and [RaSh17] upper bounds, with `sorry` proofs.

The recorded bounds are $n^2\log\log n \ll g(n) \ll n^{20/9}$. The lower bound is [ErPu71, Theorem 2], read from the page image of p. 249 because the scan's OCR layer is unusable. It takes $a=[\sqrt{\log n}]$ rows of about $n/a$ integral points each and counts triangles of area $\tfrac12 a!$, i.e. doubled area $a!$. The result is $g^{(2)}_2(n)\geqslant cn^2\log\log n$ with an unspecified $c$. The upper bound is [RaSh17, Theorem 1]. The chain before it is $O(n^{5/2})$ [ErPu71], $O(n^{7/3})$ [PaSh92], $O(n^{44/19})$ [DST09], $O(n^{9/4+\varepsilon})$ [ApSh10] and $O(n^{9/4})$ [Ap13]. [PaSh92] and [Ap13] were not read; their bounds are as [DST09] and [RaSh17] report them.

Three discrepancies between the page and its sources:

1. **Area zero.** Read literally, "triangles with the same area" includes degenerate ones: $n$ collinear points give $\binom n3$ triples of area $0$, so $g(n)=\binom n3$. The page's next sentence ("how many triangles of area $1$") and [Er75f, p. 104], which says "non zero", fix the intended reading. Everything here counts nonzero areas only.
2. **"Apfaulbaum [Ap13]".** The author is Roel Apfelbaum, as on [ApSh10] and in [RaSh17]'s reference [1].
3. **[ErPu71, p. 248].** The page image reads "The same method shows that $G^{(k)}_{2k}(kn)\geqslant n^k$ and $g^{(k)}_{2k+2}(kn+1)\geqslant n^{k+1}$", followed by the conjecture $g^{(k)}_{2k+2}(n)=\frac{n^{k+1}}{(k+1)^{k+1}}(1+o(1))$. As printed, the first bound gives $g^{(k)}_{2k+2}(N)\gtrsim N^{k+1}/k^{k+1}$, which exceeds the conjectured value, so the two cannot both hold. The problem page states the bound as $\bigl(\frac1{(k+1)^{k+1}}+o(1)\bigr)n^{k+1}$, consistent with the conjecture. Which line is misprinted is inferred, not stated in the source.

**Literature checked for the constant in §3.** [ErPu71], [Er75f], [DST09], [ApSh10] and [RaSh17] each state the lower bound only as $\Omega(n^2\log\log n)$ or $c\,n^2\log\log n$ with an unspecified constant (a loose grep of each text for "log log", with every hit read). Two web searches on September 23, 2026 found no explicit constant either. The first answer only restated the query and cited nothing that supports it; the second returned only the bound history above. OEIS and Semantic Scholar sweeps earlier the same day found nothing that improves either bound. Not read: [PaSh92], [Pu74], [Ap13], and Brass, Moser and Pach's *Research Problems in Discrete Geometry* (2005). A further search on November 2, 2026, of the grid literature ([DST09] treats only minimum-area triangles) and of OEIS (queries such as 1, 4, 7, 12, 18, 28, 36), found neither the constant $6e^\gamma/\pi^2$ nor a table of small values of $g$. Absence from what was read or searched is not evidence that the constant is new, and no novelty is claimed.

## 2. Reproducing the published grid count first

[DST09, Theorem 3] states that the $\lfloor\sqrt n\rfloor\times\lfloor\sqrt n\rfloor$ integer grid spans $(6/\pi^2-o(1))n^2$ minimum-area triangles, about $0.6079n^2$. This is a paraphrase: the arXiv text layer scrambles the formula, and the page was not rendered. In $\mathbb Z^2$ the minimum nonzero doubled area is $1$, so this count is $T_1(m)$. `src/grid_count` gives:

| $m$ | 25 | 50 | 100 | 200 | 400 | 800 |
|---|---|---|---|---|---|---|
| $T_1/n^2$ | 0.60334 | 0.60705 | 0.60758 | 0.60787 | 0.60791 | 0.60793 |

$6/\pi^2=0.607927$, and at $m=800$, $T_1=249{,}006{,}655{,}884$ and $T_1/n^2=0.607926$. This is the evidence that the engine measures the published quantity. The six histogram totals also equal OEIS A045996, the number of triangles in the $m\times m$ grid, and $\binom{m^2}3$ minus the collinear triples.

## 3. The square grid: Theorems A and B

**Status: informal proofs, written for this directory. Not machine-checked and not reviewed by anyone else.** The upper-bound step reuses the charging of [DST09, Lemma 1], which assigns each triangle to a longest side.

For $u\in\mathbb Z^2\setminus\{0\}$ let $\gcd(u)=\gcd(|u_x|,|u_y|)$. Let $P_a(m)$ be the number of unordered pairs $\{p,q\}\subset G_m$ with $\gcd(q-p)=a$.

**Lemma 1 (apexes).** Let $p\neq q$ in $\mathbb Z^2$, $u=q-p$, $a=\gcd(u)$ and $N\ge1$. The $w\in\mathbb Z^2$ with $|\det(u,w-p)|=N$ exist iff $a\mid N$, and then they fill two lattice lines parallel to $u$. Put $t(w)=(w-p)\cdot u/|u|^2$, so $p+t(w)u$ is the foot of the perpendicular from $w$. On each line exactly $a$ points have $t\in[0,1)$, and at most $a+1$ have $t\in[0,1]$.

*Proof.* Write $u=au'$ with $u'$ primitive. Then $\det(u,w-p)=a\det(u',w-p)$, and $w\mapsto\det(u',w-p)$ maps $\mathbb Z^2$ onto $\mathbb Z$ with fibres that are translates of $\mathbb Zu'$. So the solutions are the fibres over $\pm N/a$, and these exist iff $a\mid N$. Moving along a fibre by $u'$ changes $t$ by $u'\cdot u/|u|^2=1/a$, so on each line the $t$-values form a coset of $\frac1a\mathbb Z$. ∎

**Lemma 2 (pairs).** $P_a(m)=\frac3{\pi^2}\frac{m^4}{a^2}+O(m^3/a)$ uniformly for $1\le a<m$, and $P_a(m)=0$ for $a\ge m$.

*Proof.* Let $Q_b$ count pairs with $b\mid\gcd(q-p)$, that is, pairs in the same residue class mod $b$. Write $m=kb+\rho$ with $0\le\rho<b$. Each coordinate class then has $k$ or $k+1$ members, and $s_b=\sum_j c_j^2=m^2/b+\theta$ with $0\le\theta\le b/4$. So $Q_b=\frac12(s_b^2-m^2)=\frac{m^4}{2b^2}+E_b$, with $|E_b|\le m^2/2$ for $b\le m$, and $Q_b=0$ for $b\ge m$. Möbius inversion over multiples gives $P_a=\sum_{k<m/a}\mu(k)Q_{ak}=\frac{m^4}{2a^2}\sum_{k<m/a}\frac{\mu(k)}{k^2}+O(m^3/a)$, and the sum is $\frac6{\pi^2}+O(a/m)$. ∎

**Proposition 3 (upper).** For all $m,N\ge1$: $T_N(m)\le\sum_{a\mid N}2a\,P_a(m)=\frac6{\pi^2}\,m^4\,\frac{\sigma(N)}N+O(m^3d(N))$.

*Proof.* Charge each triangle to a longest side $\{p,q\}$. The angles at $p$ and $q$ are acute, so the third vertex has $t\in(0,1)$. By Lemma 1 there are at most $2a$ such vertices, with $a=\gcd(q-p)\mid N$. Summing Lemma 2 over $a\mid N$ gives $\frac6{\pi^2}m^4\sum_{a\mid N}\frac1a+O(m^3d(N))$, and $\sum_{a\mid N}\frac1a=\sigma(N)/N$. ∎

**Proposition 4 (lower).** If $3N<m$, then $T_N(m)\ge\frac6{\pi^2}\,m^4\,\frac{\sigma(N)}N-O(m^3\sigma(N))$.

*Proof.* Let $G'=\{N,\dots,m-1-N\}^2$, a copy of $G_{m'}$ with $m'=m-2N>N$. Orient each pair of $G'$ whose gcd $a$ divides $N$. By Lemma 1, exactly $2a$ lattice points $w$ have doubled area $N$ and $t(w)\in[0,1)$. Each lies within $N/|u|\le N$ of the segment $pq\subset[N,m-1-N]^2$, so it is in $G_m$. The number of such (pair, $w$) is $I=\sum_{a\mid N}2aP_a(m')$. Every divisor of $N$ is below $m'$, so Lemma 2 gives $I=\frac6{\pi^2}m'^4\frac{\sigma(N)}N+O(m^3d(N))$.

A triangle arises at most three times, once per side. An obtuse triangle arises at most once: for either side at the obtuse vertex, the foot of the third vertex falls outside the closed side. So $T_N\ge I-2\,\mathrm{NO}_N$, where $\mathrm{NO}_N$ counts non-obtuse triangles of doubled area $N$.

In such a triangle, let $\ell$ be the shortest side, opposite angle $A$, and let $B\ge C$ be the other two angles. The height onto $\ell$ is $h=\ell\sin B\sin C/\sin A\ge\ell\sin B$, because $A\le C\le90^\circ$. Since $60^\circ\le B\le90^\circ$, this gives $h\ge\frac{\sqrt3}2\ell$ and so $\ell^2\le2N/\sqrt3$. Orient the shortest side. Both base angles are at most $90^\circ$, so $t\in[0,1]$, and Lemma 1 allows at most $2(\gcd(s)+1)$ third vertices. Hence $\mathrm{NO}_N\le m^2\sum_{0<|s|^2\le2N/\sqrt3}2(\gcd(s)+1)=O(m^2N(1+\log N))$.

Finally $m'^4\ge m^4-8Nm^3$, $d(N)\le\sigma(N)$, and $N(1+\log N)\le m\,\sigma(N)$. ∎

**Theorem A.** For fixed $N$, $T_N(m)/m^4\to\frac6{\pi^2}\frac{\sigma(N)}N$ as $m\to\infty$. This follows from Propositions 3 and 4. For $N=1$ it is the limit behind [DST09, Theorem 3].

**Theorem B.** As $m\to\infty$, $\max_{N\ge1}T_N(m)=\bigl(\frac{6e^\gamma}{\pi^2}+o(1)\bigr)\,n^2\log\log n$ with $n=m^2$. Hence

$$g(n)\ \ge\ \Bigl(\frac{6e^\gamma}{\pi^2}-o(1)\Bigr)\,n^2\log\log n,\qquad \frac{6e^\gamma}{\pi^2}=1.08276\ldots$$

*Proof.* **Lower.** Take $N=\operatorname{lcm}(1,\dots,y)$ with $y=\lfloor\frac12\log m\rfloor$. By the prime number theorem, $\log N=\psi(y)=(1+o(1))y$, so $N=m^{1/2+o(1)}$ and $3N<m$ for large $m$. Since $\sigma(N)/N\le1+\log N$, the error in Proposition 4 is $o(m^4)$. Each prime $p\le y$ divides $N$ to the largest power $p^e\le y$. So $\frac{\sigma(N)}N=\prod_{p\le y}\frac{1-p^{-e-1}}{1-1/p}\ge(1-\frac1y)^{\pi(y)}\prod_{p\le y}(1-\frac1p)^{-1}$. Mertens' product [Gr13, p. 117] gives $\prod_{p\le y}(1-1/p)^{-1}\sim e^\gamma\log y$, and $(1-\frac1y)^{\pi(y)}\to1$. With $\log y=\log\log n-\log4+o(1)$ this gives the lower bound. **Upper.** $T_N(m)=0$ for $N>(m-1)^2$. Gronwall's theorem [Gr13, eq. (25)], $\limsup\sigma(N)/(N\log\log N)=e^\gamma$, gives $\max_{N<n}\sigma(N)/N\le(e^\gamma+o(1))\log\log n$. Also $d(N)=O_\varepsilon(N^\varepsilon)$, so the error in Proposition 3 is $o(m^4)$. **Non-square $n$.** $g$ is nondecreasing, since adding a point destroys no triangle. With $m=\lfloor\sqrt n\rfloor$, $m^4=n^2(1-O(n^{-1/2}))$ and $\log\log m^2=\log\log n+o(1)$. ∎

**What this does and does not say.** The order $n^2\log\log n$ is Erdős and Purdy's. Theorem B supplies an explicit constant for the square grid. Its upper half says that no single area does better in the square grid, so beating $1.0828\,n^2\log\log n$ needs a different point set. Nothing here touches the upper bound $n^{20/9}$.

**Convergence is very slow (measured, §6).** $T/(n^2\log\log n)$ at the best $N$ rises from 0.70 to 0.86 over $m=25,\dots,800$, still well below 1.0828. There are two visible reasons. First, at the best $N$ the count is 5 to 12 % below $\frac6{\pi^2}\frac{\sigma(N)}Nn^2$; in all 41 cells of the fixed-$N$ table with $N\ge2$ and $N/n\le0.1$, $1-\text{ratio}$ lies between $3.3\,N/n$ and $5.9\,N/n$. That is an empirical fit, not proved, and it fails for larger $N/n$. Second, $\sigma(N)/N$ at the best $N$ lags $e^\gamma\log\log n$: at $m=800$ it is 3.838 at $N=5040$, against 4.618.

## 4. Exact values: $g(5)=7$ to $g(10)=48$

**$g(5)=7$, by hand.** *Lower bound:* $R=\{(0,0),(1,0),(0,1),(1,1),(1,-1)\}$ spans 7 triangles of doubled area 1. Its full profile is 7 of doubled area 1, 2 of doubled area 2 and 1 collinear triple.

*Upper bound:* suppose 5 points span 8 triangles of one nonzero area. Each triangle lies in 2 of the 5 four-point subsets, and $2\cdot8>5\cdot3$. So some four points $A,B,C,D$ span 4 triangles of that area. No three of them are collinear, and none lies inside the triangle of the other three, whose area would then be a sum of three others. So they are in convex position, say in the order $ABCD$. Then $[ABC]=[ABD]$ with $C,D$ on the same side of the edge $AB$ forces $CD\parallel AB$, and likewise $AD\parallel BC$, so $ABCD$ is a parallelogram.

An affine map takes it to $(0,0),(1,0),(1,1),(0,1)$ with doubled area 1. A fifth point $(x,y)$ needs 4 of the six values $|y|,|y-1|,|x|,|x-1|,|x-y|,|x+y-1|$ to equal 1. No two in the pairs $\{|y|,|y-1|\}$ and $\{|x|,|x-1|\}$ can both be 1, so $x-y=\pm1$ and $x+y-1=\pm1$ must both hold. The four solutions all have $x=\frac12$ or $y=\frac12$, which rules out a whole axis pair. So at most 3 of the six equal 1, giving at most $4+3=7$, a contradiction. ∎

*Independent certificate:* `src/g5_certificate.py 8`. It maps one of the triangles to $(0,0),(1,0),(0,1)$. That leaves $\binom97=36$ choices of seven more triangles and $2^7$ sign vectors: 4608 polynomial systems in the other two points. Every Gröbner basis over $\mathbb Q$ is $[1]$, so none has a solution even over $\mathbb C$. The control `g5_certificate.py 7` finds consistent systems and a real configuration spanning 7.

**$g(6)=12$, computer-assisted** (`src/g6_exact.py`; exact arithmetic, 3 s).

1. *Classification.* The $K=7$ systems number $\binom96\cdot2^6=5376$, and 63 are consistent. For each, sympy's `solve` and `solve_poly_system` must return the same finite set. Its size must also equal $\dim_{\mathbb Q}\mathbb Q[x_3,y_3,x_4,y_4]/I$, read off a grevlex Gröbner basis, so no complex solution is missed. The 42 real solutions with five distinct points are each mapped onto $R$ by an explicit affine bijection. So every 5-point set with 7 triangles of one nonzero area is an affine copy of $R$.
2. *Reduction.* Take six points with $t\ge13$ such triangles. Each triangle lies in 3 of the 6 five-subsets, and $3t\ge39>36$, so some five-subset has 7 and is an affine copy of $R$ with doubled area 1. At least 6 of the 10 triangles through the sixth point $X$ then have $\det=\pm1$. No direction occurs among more than 4 pairs of $R$, so $X$ lies on two non-parallel such lines.
3. *Enumeration.* All 50 intersection points are tried in exact rationals. The largest count of doubled-area-1 triangles in $R\cup\{X\}$ is 12, at $X=(0,-1)$ and $(0,2)$; both make an affine copy of the $2\times3$ grid. That grid has 12, so $g(6)=12$.

`src/g6_brute.py` shares no code with `g6_exact.py`. It lists the possible denominators $\{1,2,3\}$ and scans every $X\in\frac16\mathbb Z^2$ with $|x|,|y|\le16$ (37,244 points) in integer arithmetic. It finds the same maximum at the same two points. During development three planted defects each made the check fail: a wrong $R$, one dropped solution (caught by the quotient-dimension test), and a shrunken brute-force box (maximum 10). The result rests on sympy's Gröbner bases and solvers. The cross-checks guard against missed solutions, but this is not a formal proof.

**$g(7)\le20$, the first step** (`src/g7_no21.py`, September 28, 2026). If seven points had 21 triangles of one nonzero area, every six-subset would have exactly 12, since each triangle is contained in four six-subsets and $g(6)=12$. Every five-subset has at most 7; when it has 7, the classification above says its equal-area triple system is $K_5^{(3)}$ with exactly the three triples through one pair missing. Enumerating the resulting 3-uniform hypergraphs on seven vertices leaves two isomorphism classes. For each class and each possible marked equal-area triple, an affine map sends that triple to $(0,0),(1,0),(0,1)$, and the polynomial equations $\det(t)^2=1$ for the other marked triples have Gröbner basis $[1]$ over $\mathbb Q$. Thus no such 21-triangle configuration exists, even over $\mathbb C$. This method does not reach 19 or 20, because the averaging no longer forces every six-subset to be extremal, and it was feared that five-point sets with 6 equal-area triangles form positive-dimensional families. They do not, as the next paragraph shows.

**$g(7)=18$ and $g(8)=28$, computer-assisted** (`src/g7_exact.py`; exact rational arithmetic, under a second; October 19, 2026). The argument is one written lemma and three finite enumerations. Throughout, a *triangle* means one of doubled area exactly 1. If the common doubled area is $D$, one affine map of determinant $1/D$, such as $(x,y)\mapsto(x/D,y)$, first makes it 1; after that, equivalence and the canonical form use only affine maps of determinant $\pm1$, which preserve areas.

*Lemma (four points).* Four distinct points span at most 4 triangles of one nonzero area $N$. If they span 4 they form a parallelogram. If they span exactly 3, three of them form a triangle of area $3N$ and the fourth is its centroid. *Proof.* If three of the points are collinear, the three triangles through the fourth point have areas proportional to the three segments on that line, the longest of which is the sum of the other two, so at most two are equal. Otherwise the four points are either in convex position, where with cyclic order $ABCD$ the two ways of cutting the quadrilateral give $[ABC]+[ACD]=[ABD]+[BCD]$, so three equal areas force the fourth, and four equal areas give $CD\parallel AB$ and $AD\parallel BC$ as in the $g(5)$ proof; or one point $D$ lies inside the triangle $ABC$, where $[ABD]+[BCD]+[CAD]=[ABC]$, so $[ABC]$ cannot equal any of the three, and three equal small triangles mean that $D$ has equal barycentric coordinates, that is, $D$ is the centroid. ∎

*Step 1: five points with at least 6 triangles.* Each triangle lies in two of the five four-subsets, so the four-subsets carry $2\cdot6=12>5\cdot2$ triangles between them and one of them has at least 3. By the lemma it is a parallelogram, taken to the unit square, or a triangle with its centroid, taken to $(0,0),(1,0),(0,3)$ with centroid $(\tfrac13,1)$. The fifth point $E$ forms a triangle with a pair $\{p,q\}$ of the four exactly when it lies on one of the two lines $\{X:\det(q-p,X-p)=\pm1\}$, so $E$ lies on at least 2, respectively 3, of these 12 lines. The script checks that no single line serves that many pairs, so two of the lines through $E$ are distinct and $E$ is their intersection: finitely many candidates, each counted exactly. Every candidate set is reduced to a canonical form, the smallest image under the affine maps that send one of its triangles to $(0,0),(1,0),(0,1)$; these maps have determinant $\pm1$, and the canonical form is a complete invariant of area-preserving affine equivalence. **Result: exactly four configurations.** $R$, with 7, and three with 6. In the canonical coordinates the script prints, the three are the parallelogram $(0,0),(1,0),(0,1),(-1,1)$ together with $(-4,2)$, with $(-3,2)$, or with $(-1,\tfrac12)$; the first two have degree sequence $(2,3,4,4,5)$, the last has $(2,4,4,4,4)$ and is a regular hexagon minus a vertex up to affine maps.

*Step 2: six points with at least 11.* The five-subsets carry $3\cdot11=33>6\cdot5$ triangles, so one has at least 6 and is one of the four configurations of Step 1. The sixth point lies on at least $11-7=4$ of the 20 lines of that configuration, no line serves 4 pairs, and the intersections are tried. **Result: exactly four configurations**, two with 12 and two with 11. The two with 12 are the $2\times3$ grid, with degree sequence $(5,5,5,5,8,8)$, and the affine regular hexagon $\{\pm u,\pm v,\pm(u-v)\}$, in which every point has degree 6. In the hexagon every five-subset has exactly 6 triangles, which is why `g6_exact.py`, which starts from a five-subset with 7, could not see it; the value $g(6)=12$ it proves is unaffected, and `small_search` found both configurations in the $5\times5$ box without distinguishing them.

*Step 3: seven points with at least 19.* If seven points span $e\ge19$ triangles, the degrees sum to $3e$, so some point has degree at most $\lfloor3e/7\rfloor$, and the six others span at least $e-\lfloor3e/7\rfloor\ge19-8=11$, since that quantity does not decrease in $e$: one of the four configurations of Step 2. The seventh point lies on at least $e-12\ge19-12=7$ of its 30 lines, no line serves 7 pairs, and the 411 intersection points give at most 18. **So $g(7)\le18$**, and the $2\times3$ grid plus $(1,2)$ gives 18: $g(7)=18$. Averaging then gives $g(8)\le\lfloor 8\cdot18/5\rfloor=28$, and the $3\times3$ grid minus its centre has 28: $g(8)=28$.

*Cross-checks.* `src/g7_box_check.py` shares no code with `g7_exact.py`. It verifies the lemma on every four-subset of the $7\times7$ box for every area; it scans every five-subset of the $6\times6$ box with at least 6 triangles of one area and every six-subset with at least 11, and finds that each has the area profile, the multiset of all 20 doubled areas in units of the common one, of a configuration in `results/g7_classes.json`, with all four of each kind present and nothing else; and it identifies the two 12-triangle configurations as a $2\times3$ grid and a centrally symmetric hexagon. A box scan can refute the classification but not confirm it, so the proof is the lemma, the three averaging steps and the exact enumeration. The lemma and the averaging steps are written, not machine-checked. The values agree with `g5_certificate.py` and `g6_exact.py`, which use Gröbner bases instead. No table of small values of $g$ was found in the sources read (§5), so whether $g(7)=18$ and $g(8)=28$ are new is not known, and no novelty is claimed.

**$g(9)=36$, computer-assisted** (`src/g9_exact.py`; exact arithmetic, about 10 s; November 1, 2026). The conventions are those of the $g(7)$ argument: a *triangle* has doubled area exactly 1 after the scaling map of determinant $1/D$, and $C(n,t)$ is the set of $n$-point sets with at least $t$ triangles, up to affine maps of determinant $\pm1$. The steps extend Steps 1 to 3 to lower thresholds. The averaging steps and the Type II reduction below are written proofs, not machine-checked; the script checks the arithmetic in them and carries out the enumerations.

*Averaging (written).* If $n$ points span $e$ triangles, the degrees sum to $3e$, so some point has degree at most $\lfloor3e/n\rfloor$, and the other $n-1$ points span at least $e-\lfloor3e/n\rfloor$ triangles. That quantity does not decrease as $e$ grows, and it is $16-6=10$, $25-9=16$ and $37-12=25$ for $(n,e)=(7,16),(8,25),(9,37)$. So every member of $C(7,16)$, $C(8,25)$ and $C(9,37)$ is a member of $C(6,10)$, $C(7,16)$ and $C(8,25)$, respectively, plus one point. If the base has $c$ triangles, the new point lies on at least $t-c$ of its lines. The script checks that no line of any base serves that many pairs, so the new point is the intersection of two distinct lines, as in Steps 1 to 3.

*Step 4: six points with at least 10.* If $e\ge11$, some five-subset has at least 6, since $3e>30$. If $e=10$ and some five-subset has at least 6 (Type I), the sixth point lies on at least $10-c$ lines of a Step 1 configuration, no line serves that many pairs, and the intersections are tried. Otherwise (Type II) every five-subset has at most 5. The five-subset counts sum to $3\cdot10=30$, so each is exactly 5 and every point has degree 5. The script lists the 3-uniform hypergraphs on six vertices with 10 edges and all degrees 5, 7 up to isomorphism. For each, one edge is sent to $(0,0),(1,0),(0,1)$, and the other nine edges give $\det=\pm1$. Of the $7\cdot2^9=3584$ sign systems, 4 have a Gröbner basis over $\mathbb Q$ other than $[1]$. All four are zero-dimensional, and each has exactly one distinct solution, which equals its quotient dimension of 1, so no complex solution is missed. Each solution is rational, satisfies its system exactly on substitution, and has all ten edges confirmed again in rational arithmetic. A solution would be discarded only if a coordinate were provably non-real, and the run stops if realness cannot be decided. Each solution is the affine regular hexagon, which has 12 triangles, more than the 10 that Type II allows, so no six-point set is of Type II. **Result: $C(6,10)$ has 15 configurations:** the four of Step 2 and eleven with exactly 10.

*Steps 5 to 7.* $C(7,16)$ has 19 configurations: 3 with 18, 4 with 17 and 12 with 16. The three with 18 include the $2\times3$ grid plus $(1,2)$, and none has 19, as Step 3 requires. $C(8,25)$ has 6: one with 28, one with 26 and four with 25. The only one with 28 is the $3\times3$ grid minus its centre, so **every eight-point set with 28 triangles of one area is an affine image of the $3\times3$ grid minus its centre.** No intersection point extends any of the six beyond 36, so $C(9,37)$ is empty. **So $g(9)\le36$**, and the stored configuration of §5, the $3\times3$ grid minus its centre plus $(3,0)$, has 36: $g(9)=36$. Averaging then gave $g(10)\le\lfloor36\cdot10/7\rfloor=51$, which the next step closes.

*Cross-checks.* `src/g9_box_check.py` shares no code with `g7_exact.py` or `g9_exact.py`, and its box scans use a new C program, `src/box_scan.c`. It makes four checks:

1. It recounts every listed configuration in integer arithmetic.
2. It scans every 6-subset of the $6\times6$ box with at least 10 triangles of one area, every 7-subset with at least 16, every 8-subset with at least 25 and every 9-subset with at least 37. Every area profile found is listed, and every listed profile is found. Two of the 19 seven-point configurations share a profile, so the scan sees 18 profiles there and cannot tell those two apart. No 9-subset reaches 37.
3. It repeats Steps 4 (Type I) to 7 from the listed configurations of each previous level with its own integer code. The profiles it reaches equal the listed ones at every level, and it finds no line serving as many pairs as the new point needs.
4. It asks z3's nonlinear real arithmetic, for each of the 7 Type II hypergraphs, for six distinct real points whose triangles of doubled area 1 are exactly those edges. All 7 are unsatisfiable. As a control, the hexagon's 12-edge hypergraph is satisfiable.

Checks 2 and 4 can refute the classification but not confirm it, and check 3 guards against implementation errors, not against an error in the argument. The result rests on sympy's Gröbner bases and solver, cross-checked by z3 for Type II. It is not a formal proof. No table of small values of $g$ was found in the sources read, so whether $g(9)=36$ is new is not known, and no novelty is claimed.

**$g(10)=48$, computer-assisted** (`src/g10_exact.py`; exact arithmetic and the SAT solver kissat, about 45 s; November 2, 2026). The conventions are those of the $g(9)$ argument. The obstacle recorded after $g(9)$ was that the averaging chain for ten points runs down to $C(6,9)$, which is infinite. The argument below never classifies $C(6,9)$ or $C(7,15)$. It needs $C(8,24)$, and it handles the one case of $C(8,24)$ that averaging does not reduce to a finite class with a lemma about hypergraphs.

*Growth rules.* A rule starts from a class $C(m,t_m)$ that is already classified completely and adds points one at a time; at each step the new point must form at least $a$ triangles with pairs of the earlier points, and the enlarged set must have at least $t$. Before every step the script checks that no line of any base serves as many pairs as the new point needs, $\max(a,t-c)$ for a base with $c$ triangles, so the new point is the intersection of two distinct lines, as in Steps 1 to 3. The four rules:

| rule | start | steps, as (points: new point $\ge a$, total $\ge t$) |
|---|---|---|
| a2 | $C(7,16)$ | 8: total $\ge24$ |
| a3 | $C(6,10)$ | 7: total $\ge15$; 8: total $\ge24$ |
| b5 | $C(5,6)$ | 6: $a\ge3$, $t\ge9$; 7: total $\ge15$; 8: total $\ge24$ |
| b4 | $C(4,3)$ | 5: $a\ge2$, $t\ge5$; 6: $a\ge3$, $t\ge8$; 7: $a\ge4$, $t\ge15$; 8: total $\ge24$ |

$C(4,3)$ is the unit square and the triangle $(0,0),(1,0),(0,3)$ with its centroid, by the four-point lemma.

*Eight points (written).* An eight-point set with $e\ge24$ triangles is in $C(8,25)$ if $e\ge25$. If $e=24$, the degrees sum to 72. A point of degree at most 8 leaves seven points with at least 16 (rule a2); otherwise every degree is exactly 9.

*Lemma (machine-checked with certified proofs, not written).* Every 3-uniform hypergraph on 8 vertices with 24 edges and every degree 9 contains a chain of rule a3, b5 or b4: nested vertex sets whose edge counts and new-vertex degrees meet the rule's bounds. The lemma uses no geometry. `g10_exact.py` encodes its negation as CNF with python-sat, after the symmetry reduction below, and kissat reports it unsatisfiable. A DRAT proof of the same CNF from kissat is accepted by drat-trim, and the LRAT proof drat-trim derives from it is accepted by `lrat-check` (*Certificates*, below).

*Symmetry reduction (written).* The encoded formula fixes $\{0,1,2\}$, $\{0,1,3\}$ and $\{0,1,4\}$ as edges. This is sound for three reasons. First, whether a hypergraph has a chain of a rule depends only on edge counts and degrees inside vertex subsets, so relabelling the vertices maps chains to chains, and a relabelled hypergraph without a chain is still without one. Second, every hypergraph in question has a pair in at least 3 edges, because the 28 codegrees sum to $3\cdot24=72>2\cdot28$, and the at least 3 edges through that pair have distinct third vertices. Third, relabelling the pair as 0 and 1 and three of those third vertices as 2, 3 and 4 gives an isomorphic hypergraph that satisfies the fixed clauses. These are unit clauses that force three edges and forbid nothing, so every other edge stays free. So a hypergraph without a chain would give a satisfying assignment of the reduced formula; the reduced formula is unsatisfiable, so there is none. The reduction keeps at least one labelling of every isomorphism class. Unlike lex-leader clauses, it does not have to keep exactly one, so no canonical form is assumed.

*Why the encoding cannot lose a counterexample (written).* The CNF has one variable per triple, and three kinds of auxiliary variables:

- totalizer outputs, one for each vertex set $X$ (or each pair of a set and a vertex $y$) and threshold $t$;
- exact-count variables for the total and the degrees;
- a chain variable $c[X]$ for each rule and vertex set.

python-sat's totalizer clauses only force an output true when at least $t$ of its inputs are true; they never force one false. The chain clauses only say that $c[X]$ is true when its base condition holds, or when $c[X\setminus y]$ holds and the step's two thresholds hold; a final clause says $c$ of all 8 vertices is false. Take any genuine 9-regular 24-edge hypergraph with no chain that satisfies the relabelling. Give each triple variable its edge value, each totalizer output the value of its actual threshold predicate, each exact-count variable its actual value, and each $c[X]$ the value "some prefix chain of the rule ends at $X$". Every clause is then true: the totalizer and counter clauses state facts about the actual counts, the chain clauses state that a prefix chain extends by one valid step, and the final clause holds because the hypergraph has no full chain. So an UNSAT answer, certified or not, rules out every such hypergraph. The one-sidedness can only make the formula harder to satisfy, never easier. The cross-check's own counters are two-sided for the exact counts and one-sided for the thresholds, and the same argument applies to them.

Three controls guard the encoding:

- the degree constraints alone are satisfiable, with a model checked to be 9-regular with 24 edges;
- with rule b4 left out the formula is satisfiable, and a direct search over vertex subsets, written separately from the CNF, finds that the model has no a3 or b5 chain and does have a b4 chain;
- every eight-point configuration with exactly 24 triangles found below has a chain of one of the four rules under that search.

Without the relabelling kissat also reports the formula unsatisfiable, in 95 s (`g10_exact.py --full-lemma`, not part of the gate; that run has no certificate). z3, on a separate encoding, and cadical, on this one, did not finish the unsplit lemma in 5 and 3 minutes respectively, which is why the cross-check below splits it into cases.

So every member of $C(8,24)$ is reached by one of the four rules. **Result: $C(8,24)$ has 24 configurations**, the six of $C(8,25)$ and 18 with exactly 24, and each rule on its own reaches the same 24. None of the 18 is 9-regular, so the case the lemma handles does not occur in the plane.

*Nine and ten points (written).* A nine-point set with $e=35$ or $36$ triangles has a point of degree at most $\lfloor3e/9\rfloor$, and deleting it leaves at least $e-\lfloor3e/9\rfloor=24$. The new point then needs at least $35-28=7$ triangles, and no line of a member of $C(8,24)$ or $C(9,35)$ serves more than 4 pairs. **$C(9,35)$ has 3 configurations**: two with 36 and one with 35. So, up to affine maps, exactly two nine-point sets span 36 triangles of one area, and the stored configuration of §5 is one of them. A ten-point set with $e\ge49$ has a point of degree at most $\lfloor3e/10\rfloor$, and $e-\lfloor3e/10\rfloor\ge35$ for $e\ge49$. So it is a member of $C(9,35)$ plus a point forming at least $49-36=13$ triangles with its pairs, and none of the 1,293 intersection points does. **So $g(10)\le48$**, and the stored ten-point configuration of §5 has 48: $g(10)=48$. The largest count among all intersection-point extensions of $C(9,35)$ is 48, attained by one configuration up to affine maps, the stored one. Ten-point sets with 48 none of whose nine-subsets has 35 triangles are outside this enumeration, so it does not show that the stored configuration is the only one with 48. Averaging then gives $g(11)\le\lfloor48\cdot11/8\rfloor=66$.

*Cross-checks.* `src/g10_check.py` shares no code with `g7_exact.py`, `g9_exact.py` or `g10_exact.py`; it reuses the integer helpers of `g9_box_check.py`. It makes three checks:

1. It recounts every listed configuration in integer arithmetic.
2. It repeats the four rules in integer arithmetic from the stored classes and its own four-point bases. Every rule reaches exactly the listed eight-point area profiles; two of the 24 configurations share a profile, so it sees 23. Extending those gives exactly the three nine-point profiles, no extension of them reaches 49, the largest reaches 48, and no line serves as many pairs as a new point needs.
3. It checks the lemma with its own CNF encoding, with its own sequential counters and chain clauses, and a different solver, cadical. The hypergraphs are split into 21 cases, a reduction of the same kind as the one above (written). If some pair has codegree $k\ge4$, then $k\le6$, since only 6 other vertices remain; the pair is relabelled $\{0,1\}$ and its third vertices $2,\dots,k+1$, giving the three cases $k=4,5,6$. Otherwise every codegree is at most 3. Then the link of vertex 0 (the pairs $\{u,v\}$ with $\{0,u,v\}$ an edge) is a graph on the other 7 vertices with 9 edges, since vertex 0 has degree 9, and maximum degree at most 3. Relabelling vertices 1 to 7 makes it equal to one of 18 representatives, which come from an enumeration of all labelled graphs with those properties. All 21 cases are unsatisfiable. With rule b4 left out, 17 are satisfiable.

`g10_check.py --box` also compares every eight-point subset of the $6\times6$ box that has at least 24 triangles of one area with the list. It is not part of the gate; its one run is reported in §7. Check 2 guards against implementation errors, not against an error in the written steps.

*Certificates* (`src/g10_certify.py`, `results/g10_lrat.json`). Each of the 22 CNFs, the gate's and the 21 of the cross-check, is written in DIMACS with the clauses in the order its builder emits them, and its SHA-256 is recorded. The 21 cross-check CNFs are certified by CaDiCaL 3.0.1 with `-q --no-binary --lrat`, whose LRAT proof is streamed through a pipe into `lrat-check` from drat-trim at commit 2e3b2dc0, built under `refs/` by `sh refs/fetch.sh drat-trim`; no proof is stored. In every case CaDiCaL exited 20 (unsatisfiable), and `lrat-check` exited 0 and printed `c VERIFIED`. These proofs have 923,291,610 bytes together: 409,141,626, 213,316,952 and 151,020,723 for codegree 4, 5 and 6, and 5.9 to 11.9 million for each link case. CaDiCaL with proof logging had not finished the gate's CNF after 16 min 51 s of CPU time, so that CNF goes by a second route. kissat 4.0.4 (`-q --no-binary`) writes a DRAT proof of 96,619,305 bytes in 23 s. drat-trim, at the same commit, checks it by backward checking in 17 s, printing `s VERIFIED`, and writes the LRAT form of the checked proof, 76,827,961 bytes. `lrat-check` accepts that LRAT proof in 2 s. The proof files are deleted afterwards. As a control, the proof of one link case, written to a file, is accepted, and the same proof with one addition that a later step uses deleted is rejected. The gate does not repeat the proofs. `g10_certify.py --verify` rebuilds all 22 CNFs from `g10_exact.py` and `g10_check.py` in 8 s and requires each SHA-256 to match an accepted record.

*Trust base.* $g(10)=48$ rests on three things:

- exact rational arithmetic for all the geometry, that is, the growth enumerations and every count, redone in integer arithmetic by the cross-check;
- SAT solving with LRAT proofs checked by `lrat-check` for the hypergraph lemma, with the CNF encodings trusted to state the lemma (the controls above and the two independent encodings are the evidence that they do);
- written steps: the four-point lemma, the averaging steps and the Type II reduction behind $g(7)=18$ and $g(9)=36$ (Type II also rests on sympy's Gröbner bases, cross-checked by z3); the fact that a new point needing more triangles than any line serves lies on two distinct lines; the eight-point case split; the symmetry reduction and the 21-case split; and the nine- and ten-point averaging steps.

It is not a formal proof. No table of small values of $g$ was found (§1), so whether $g(10)=48$ is new is not known, and no novelty is claimed.

## 5. Small $n$

`src/small_table.py` ran in 2 min 18 s on 10 threads. On November 1, 2026 it was rerun on 3 threads in 8 min 2 s, after $g(9)=36$ was added. Every configuration came out identical, and only the upper bounds and the exact flag at $n=9$ changed. On November 2, 2026, after $g(10)=48$ was added, it was rerun again on 3 threads in 5 min 13 s; again every configuration was identical, and only the upper bounds and the exact flag at $n=10$ changed. It uses two methods:

- `src/small_search` searched exhaustively over every $n$-subset of the $5\times5$ box ($n\le10$) and the $6\times6$ box ($n\le12$).
- Simulated annealing ran in the $10\times10$ box for $n=4,\dots,16$: 1M iterations, 4 restarts, seeds 1 to 3, and target doubled areas $\{1,2,3,4,6,8,12\}$. A second annealing pass used the $12\times12$ box with 4M iterations and seeds 4 to 6 for $n=11,\dots,16$. That makes 399 annealing runs.

Every run's configuration was recounted by `src/recount.py` before it could be kept. As a positive control, both exhaustive boxes rediscover $g(4)=4$, $g(5)=7$, $g(6)=12$, $g(7)=18$, $g(8)=28$, $g(9)=36$ and $g(10)=48$. Annealing reaches every exhaustive box optimum, and the second pass reproduced every first-pass value without improving any.

The upper bounds for $n\ge11$ come from averaging: each triangle of an $n$-set lies in $n-3$ of its $(n-1)$-subsets, so $g(n)\le\lfloor g(n-1)\,n/(n-3)\rfloor$, started from $g(10)=48$ (§4). Started from $g(7)=18$ it gives $g(8)\le28$, which meets the lower bound, but only $g(9)\le42$, which is why $g(9)$ and $g(10)$ needed the classifications of §4. Beyond $n=10$ the bounds are weak, and $g(11)\in[60,66]$ is the first open value.

| $n$ | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|
| lower bound | 18 | 28 | 36 | 48 | 60 | 74 | 91 | 108 | 128 | 148 |
| its doubled area | 1 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| found by | both | both | both | both | both | both | anneal | anneal | anneal | anneal |
| upper bound | 18 | 28 | 36 | 48 | 66 | 88 | 114 | 145 | 181 | 222 |

For $n=7$ the configuration is the $2\times3$ grid plus $(1,2)$; for $n=8$ it is the $3\times3$ grid minus its centre. The rest are in `results/small_n.json`. The exhaustive values are optimal only within their box. An optimal set may need a larger box or non-lattice coordinates, so every lower bound for $n\ge11$ is only that; for $n=9$ and $n=10$ the box optima 36 and 48 are optimal by §4. No published table of small values was found; an OEIS search for 1, 4, 7, 12, 18 returned four unrelated sequences.

**Searches beyond the small boxes** (`src/lower_search.py`, November 1, 2026; 25 min on 3 threads; `results/lower_search.json`). None of them improves any lower bound in the table. A triangular-lattice or affine-hexagon configuration is an affine image of a subset of $\mathbb Z^2$, and affine maps preserve ratios of areas, so square boxes already contain those configurations. What the new searches add is shape, size and unbounded rational coordinates:

- *Exhaustive strips and a wider box*, with `src/box_scan`, at every doubled area. No 10-subset of the $3\times12$ box, the $4\times9$ box or the $7\times7$ box ($5{,}529{,}123{,}446$ subsets up to translation) has more than 48 triangles of one area. No 11-subset of the $3\times12$ or $4\times9$ box has more than 60, and no 12-subset of either has more than 74.
- *Annealing in larger boxes.* 196 runs of `src/small_search` in the $16\times16$ and $20\times20$ boxes for $n=10,\dots,16$: 2M iterations, 4 restarts, seeds 7 and 8, and targets $\{1,2,3,4,6,8,12\}$. They reach the table's values and none higher.
- *A box-free exact beam search.* This starts from the 15 configurations of $C(6,10)$ (§4). At each step it adds every point on two lines of the current set, in exact rationals and with no bound on the coordinates, and keeps the 150 best classes up to area-preserving affine maps. It reaches 18, 28, 36, 48, 60 and 74 at $n=7,\dots,12$ and nothing higher.

Every kept configuration was recounted by `src/recount.py`. These are searches, and they show nothing about optimality for $n\ge11$.

**Structured constructions** (`src/constructions.py`, November 2, 2026; 28 s on 3 processes; `results/constructions.json`). Two families were scanned exhaustively over small parameters and counted exactly at every doubled area. *Rows*: 2, 3 or 4 parallel rows of consecutive integers, at heights $0,1,\dots$ or $0,1,3,(4)$, with lengths summing to $n$ and offsets 0, 1 or 2; these are the unions of parallel lines with arithmetic structure. The [ErPu71] construction takes $\lfloor\sqrt{\log n}\rfloor$ rows, which is a single collinear row for every $n\le16$, so it says nothing at this size. *Lattice sections*: the $n$ points of $\mathbb Z^2$ smallest under a positive definite quadratic form with coefficients up to 5, centred at a point of $(\frac14\mathbb Z)^2$, which covers round triangular and hexagonal patches up to affine maps. Neither family reaches the table at any $n$:

| $n$ | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|
| rows | 44 | 56 | 70 | 86 | 105 | 124 | 146 |
| lattice sections | 44 | 56 | 70 | 85 | 104 | 124 | 146 |
| table | 48 | 60 | 74 | 91 | 108 | 128 | 148 |

Non-lattice constructions, such as regular polygons or sets with coordinates in $\mathbb Q(\sqrt5)$, were not tried.

## 6. Finite lattice data

`src/lattice_table.py` writes `results/lattice.json`, and `src/verify_lattice.py` rechecks it. One full run took 12 min 41 s wall-clock on 10 threads on a shared machine. Every count is exact, and each best count is a certified lower bound on $g(n)$ at that $n$. For example, $g(640{,}000)\ge911{,}639{,}642{,}036$.

**Square grids at the best doubled area $N$:**

| grid | $n$ | best $N$ | $\sigma(N)/N$ | $T_N/n^2$ | ratio to $\frac6{\pi^2}\frac{\sigma(N)}N$ | $T_N/(n^2\log\log n)$ |
|---|---|---|---|---|---|---|
| 25×25 | 625 | 12 | 2.333 | 1.3072 | 0.922 | 0.7020 |
| 50×50 | 2,500 | 60 | 2.800 | 1.5014 | 0.882 | 0.7298 |
| 100×100 | 10,000 | 120 | 3.000 | 1.7110 | 0.938 | 0.7706 |
| 200×200 | 40,000 | 360 | 3.250 | 1.8797 | 0.951 | 0.7963 |
| 400×400 | 160,000 | 2,520 | 3.714 | 2.0572 | 0.911 | 0.8284 |
| 800×800 | 640,000 | 5,040 | 3.838 | 2.2257 | 0.954 | 0.8584 |

**Theorem A at fixed $N$.** The table gives $T_N/(\frac6{\pi^2}\frac{\sigma(N)}Nn^2)$, which should tend to 1:

| $N$ | $m=25$ | 50 | 100 | 200 | 400 | 800 |
|---|---|---|---|---|---|---|
| 1 | 0.9925 | 0.9986 | 0.9994 | 0.9999 | 1.0000 | 1.0000 |
| 2 | 0.9864 | 0.9965 | 0.9991 | 0.9997 | 0.9999 | 1.0000 |
| 6 | 0.9617 | 0.9898 | 0.9974 | 0.9994 | 0.9998 | 1.0000 |
| 12 | 0.9215 | 0.9799 | 0.9949 | 0.9987 | 0.9997 | 0.9999 |
| 60 | 0.6027 | 0.8820 | 0.9721 | 0.9934 | 0.9984 | 0.9996 |
| 120 | 0.3473 | 0.7720 | 0.9382 | 0.9857 | 0.9967 | 0.9992 |
| 360 | 0.0241 | 0.4522 | 0.8172 | 0.9514 | 0.9882 | 0.9974 |
| 840 | 0.0000 | 0.1468 | 0.6228 | 0.8857 | 0.9706 | 0.9928 |
| 2520 | 0.0000 | 0.0000 | 0.2414 | 0.6970 | 0.9111 | 0.9773 |

**Shape at $n\approx40{,}000$, and the [ErPu71] strip.** Best $T_N/n^2$ for rectangles:

| grid | 200×200 | 100×400 | 40×1000 | 20×2000 | 10×4000 | 5×8000 | 3×13333 |
|---|---|---|---|---|---|---|---|
| best $N$ | 360 | 360 | 360 | 360 | 120 | 12 | 2 |
| $T_N/n^2$ | 1.8797 | 1.8714 | 1.8046 | 1.6791 | 1.4401 | 1.1059 | 0.7777 |

At every $n$ computed here, $\lfloor\sqrt{\log n}\rfloor=3$, so the [ErPu71] construction is a 3-row strip with doubled area $3!=6$. Its $T_6/n^2$ is 0.7768 at $n=9{,}999$, 0.7775 at $39{,}999$ and 0.7777 at $159{,}999$. The best area in the strip, $N=2$, is marginally higher. The square beats the strip by 2.42× at $n\approx40{,}000$ and 2.65× at $n\approx160{,}000$. This is a finite-$n$ comparison only. The strip's asymptotic constant was not derived, and its row count grows with $n$. Reading $\log$ as $\log_2$ would give 4 rows at $n=159{,}999$. No 4-row strip was computed. The 3-, 5- and 10-row rectangles at $n\approx40{,}000$ all do far worse than the square: 0.78, 1.11 and 1.44 against 1.88.

## 7. How the numbers are checked

`src/final_check.py` re-derives every entry of `results/claims.json`, runs 104 checks, and exits 0 only if all pass. On November 1, 2026, with 77 checks, it took 3 min 34 s on the shared M3 Pro at a load average between 10 and 17, and 7 min 0 s at a load average near 20; the five runs of `g9_exact.py` and `g9_box_check.py` account for about 75 s of it. On November 2, 2026, with the $g(10)$ checks and the constructions added, it took 7 min 19 s at a load average between 8 and 10, and with the certificate check and the `box_scan` guard added, 6 min 55 s at a load average near 4. The checks:

- `grid_count` must match the independent triple enumeration `src/brute_check.py` bin by bin at 12 sizes. These go up to 30×30 (120,464,616 triangles) and 17×53 (120,727,840).
- `src/brute_count.c`, a second triple enumeration, written in C and sharing no code or method with `grid_count`, must reproduce the stored 25×25, 50×50 and 100×100 histograms byte for byte.
- `grid_count` must rebuild five of the histograms to the SHA-256 that `lattice.json` records. `verify_lattice.py` must pass on all 16 runs, which includes checking each raw file against that hash.
- The $g(5)$ certificate and its control must pass, as must `g6_exact.py`, `g6_brute.py`, `g7_no21.py`, `g7_exact.py` (whose output must equal the stored `results/g7_classes.json`), `g7_box_check.py`, `g9_exact.py` (whose output must equal the stored `results/g9_classes.json`), `g9_box_check.py`, `g10_exact.py` (whose output must equal the stored `results/g10_classes.json`), `g10_check.py` and the recount of every small-$n$ configuration. Every configuration kept by `lower_search.py` or `constructions.py` must also recount, and none may exceed the claimed lower bound at its $n$.
- Each claimed value is then compared with what those steps produced. Theorems A and B are reported as NOT MACHINE-CHECKED.

The gate plants eighteen defects of its own, and each must make the component exit 1:

- 6 extra triangles in one bin, against `brute_check`;
- 6 extra triangles at $D=12$ in the 25×25 histogram with its hash updated, against `verify_lattice`;
- a wrong reference configuration, against `g6_exact`;
- a wrong centroid in a four-point base, against `g7_exact`;
- a classification with one six-point configuration dropped, against `g7_box_check`;
- an averaging step claimed for 9-point sets with 36 triangles instead of 37, against `g9_exact`;
- the Type II solutions dropped, against `g9_exact` (caught by the quotient-dimension test);
- one Type II solution coordinate perturbed by $1/1000$, against `g9_exact` (caught by substituting each solution back into its system);
- a classification with one eight-point configuration dropped, against `g9_box_check`;
- rule b4 left out of the hypergraph lemma, against `g10_exact` (kissat then finds the formula satisfiable);
- the ten-point enumeration run at threshold 48 instead of 49, against `g10_exact`;
- a classification with one 24-triangle eight-point configuration dropped, against `g10_check`;
- one LRAT proof record removed from `g10_lrat.json`, against `g10_certify --verify`;
- a recorded CNF hash that does not match the generated formula, against `g10_certify --verify`;
- a box size and threshold passed to `box_scan` as one argument, and a threshold below the subset size, against `box_scan`'s argument guard;
- an overstated small-$n$ count, against `recount`;
- an overstated count in `lower_search.json`, against `recount`.

Separately, before `g7_no21.py` was added, a copy of the directory was given three false claims in `claims.json`: $g(6)=13$, $19\le g(7)$, and one extra triangle at $n=640{,}000$. The same copy had 6 triangles added at $D=7$ in the stored 100×100 histogram. The old gate exited 1, with 44 of 49 checks passing. It named the three claims, the `brute_count` comparison, and `verify_lattice`, which found that file's hash, its total, and its total against A045996 all wrong. The rebuild check passed, correctly: it compares `grid_count`'s output with the recorded hash, not with the file.

`box_scan` refuses arguments that are not single integers and thresholds below the subset size, and aborts with exit 2 after 1,000,000 output lines. These guards were added on November 2, 2026, after a shell invocation whose size and threshold reached `box_scan` as one argument ran a scan at threshold 1 and wrote 31 GB of output before it was stopped and deleted; legitimate runs print a few thousand lines at most. On November 2, 2026, `g10_check.py --box`, which the gate does not run, passed in 2 min 39 s. The $6\times6$ box holds 23 area profiles of eight-point sets with at least 24 triangles of one area, and all 23 are listed profiles, so every listed eight-point profile occurs in the box.

**Known blind spot.** `verify_lattice.py` checks totals, the best and fixed-$N$ bins, and hashes. A planted move of 6 triangles between two unchecked bins (7 and 11), with the hash updated, **passed**. Five of the 14 stored grids have every bin confirmed independently. The gate has `brute_count` reproduce the 25×25, 50×50 and 100×100 grids. One-off runs logged in `results/brute_count.log`, and not gate-checked, reproduced the 200×200 grid (28 min) and the 3×3333 strip byte for byte. The other nine were not brute-forced: 400×400, 800×800, the five other rectangles at $n=40{,}000$, and the strips at $n=39{,}999$ and $159{,}999$. The cost grows as $n^3$, so 400×400 would take an estimated 21 to 30 hours at the rates measured here. For those nine, per-bin correctness rests on `grid_count`'s agreement at the smaller sizes, on every raw count being divisible by 6, which `grid_count` enforces, and on the totals.

## 8. Not done

- No Lean formalization. formal-conjectures now states the problem (§1), and the small values here could be stated against `Erdos1086.g`; none were.
- Theorems A and B are informal and unreviewed. The constant has not been checked against the unread literature in §1.
- The asymptotics of the [ErPu71] strip were not derived.
- The method of §4 stops at $g(10)$. To bring $g(11)$ below 66 it must rule out 66 triangles on 11 points, and the averaging chain for that runs $66-\lfloor198/11\rfloor=48$, so $C(10,48)$; then $48-14=34$, $34-11=23$, $23-8=15$ and $15-6=9$, ending at $C(6,9)$, which is infinite: $R=\{(0,0),(1,0),(0,1),(1,1),(1,-1)\}$ together with any point $(2,t)$ has at least 9 triangles of doubled area 1, and exactly 9 unless $(2,t)$ also lies on one of the finitely many other lines of $R$, because the line $x=2$ serves two pairs of $R$ (`g9_exact.py` checks this at 101 values of $t$). Whether growth rules and a hypergraph lemma can avoid it, as they did for $C(8,24)$ in the proof of $g(10)=48$, was not tried. Parametrised families (a new point confined to one line, counted as a function of its position) were not needed for $g(10)$ and were not implemented.
- $g(11)$ remains between 60 and 66, and the small-$n$ lower bounds for $n\ge11$ come from searches (§5). The four-point lemma, the averaging steps, the Type II reduction and the eight-point case split behind $g(7)=18$, $g(9)=36$ and $g(10)=48$ are written proofs, and the hypergraph lemma behind $g(10)=48$ rests on SAT answers certified by LRAT proofs, with the CNF encodings trusted to state it.
- The problem itself remains open.
