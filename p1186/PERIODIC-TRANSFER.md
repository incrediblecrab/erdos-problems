# A direct periodic transfer from cyclic colorings

**October 7, 2026.** The safe-coset obstruction in the initial OpenAI transfer study does not obstruct the direct periodic construction. That construction already appears in this repository and in [Lu and Peng, Lemma 1 and equation (10)](https://arxiv.org/abs/1107.2888). The follow-up formalizes its divisibility bridge, derives an exact finite count, and applies it to OpenAI's stated cyclic-coloring growth theorem. The resulting asymptotic statement is conditional on that upstream theorem, not a new unconditional resolution of problem 1186.

## Statement and count

Let $k\ge2$, let $m\ge1$, and let $c:\mathbb Z/m\mathbb Z\to\{0,1\}$ have no monochromatic $k$-term cyclic progression with nonzero difference. Extend it to the nonnegative integers by $C(a)=c(a\bmod m)$.

**Divisibility bridge.** For every $a,d\ge0$, the sequence $a,a+d,\ldots,a+(k-1)d$ is monochromatic under $C$ if and only if $m\mid d$. If $m\nmid d$, reduction modulo $m$ contradicts the cyclic hypothesis. If $m\mid d$, all terms have the same residue, so they have the same color. This direction does not require any wildcard property.

Count progressions contained in $\{0,\ldots,N-1\}$ with $d\ge0$, as the repository's existing finite-count convention does. Write

$$D=(k-1)m,\qquad q=\lfloor N/D\rfloor,\qquad r=N-qD,\quad 0\le r<D.$$

The bridge permits exactly the differences $d=mj$. For each such difference the number of starts is $\max(N-Dj,0)$. Including a zero contribution at the upper endpoint when $D\mid N$ is harmless. Thus

$$M(N)=\sum_{j=0}^{q}(N-Dj)=N(q+1)-\frac{Dq(q+1)}2.$$

Substituting $N=Dq+r$ gives the exact identity

$$M(N)=\frac{N^2}{2D}+\frac N2+\frac{r(D-r)}{2D}.$$

For $N>0$, therefore,

$$\frac{M(N)}{N^2}=\frac1{2D}+\frac1{2N}+\frac{r(D-r)}{2DN^2},\qquad 0\le\frac{r(D-r)}{2DN^2}\le\frac{D}{8N^2}.$$

If only increasing progressions are counted, subtract the $N$ zero-difference progressions. The finite boundary term becomes $-1/(2N)$ rather than $+1/(2N)$, and the limit is unchanged. Since $\delta_k$ is the optimal guaranteed coefficient over all two-colorings, this one coloring proves

$$\delta_k\le\frac1{2(k-1)m}.$$

This is an upper bound from a construction, not equality with the optimum. The former normalization wording in `NOTES.md` is corrected accordingly.

The interval limit is taken at fixed $k$ and $m$. Only afterward does the conditional corollary let $k$ grow. The error depends on $D=(k-1)m$, so this is not a uniform finite-interval estimate for a progression length growing with $N$.

## Conditional consequence of the new source

The [pinned OpenAI source, theorem `perturb:cyclic`](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/Quantitative-Superexponential-Bounds-for-van-der-Waerden-Numbers-September-23-2026/build/sections/06-perturbation.tex), states that for every sufficiently large $k$ there is such a cyclic two-coloring with $m\ge k^{ck}$, where $c=10^{-5}$. Its existence claim, not the later digit-product construction, is the imported premise.

Combining the premise with the periodic bound yields

$$\delta_k\le\frac1{2(k-1)k^{k/100000}}\quad\text{for all sufficiently large }k.$$

Taking $k$th roots bounds $\delta_k^{1/k}$ by $[2(k-1)]^{-1/k}k^{-1/100000}$, which tends to zero. This is a conditional superexponential upper bound as progression length grows, not an asymptotic formula for $\delta_k$ and not an improvement to any small-length record in `results/claims.json`.

The constant is small and the source threshold has not been made effective here. The derived expression is eventually smaller than the random-coloring bound, but this does not provide a useful numerical crossover in the range of the existing search. No claim is made that this corollary is the first such statement in the literature. The periodic transfer itself is prior work.

The source's [selected Comparator challenge](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/lean/ComparatorChallenges/QuantitativeVanDerWaerden.lean) states an interval van der Waerden lower bound. A long interval coloring does not automatically give a cyclic coloring. Consequently that challenge, even if accepted, would not by itself discharge the stronger cyclic premise used here. The source may have internal formal lemmas beyond that interface; none was verified in this study.

## What is checked locally

`src/PeriodicTransfer.lean` proves the divisibility bridge for every base satisfying the stated cyclic hypothesis. Its `fourColor_ap_free` theorem proves that the hypothesis holds for the word `0011` modulo four at length three, so the generic result is not supported only by an impossible assumption. The finite summation, limit argument, and application to the imported growth theorem above are written arguments, not claimed as Lean formalizations.

`src/openai_transfer.py` enumerates finite interval progressions independently of the closed formula, including intervals shorter than one spacing and intervals on both sides of divisibility boundaries. It checks the exact remainder identity and compares the resulting density with the existing `exact.psi_zm` evaluator. The modulo-four and quadratic-residue modulo-eleven cases reproduce their periodic coefficients; neither beats the corresponding best stored construction.

`src/check_periodic_lean.py` uses the existing shared Erdős probe, parser, and axiom policy. It compiles the source, inspects the compiled declarations through a probe elaborated outside the candidate environment, and replays the candidate through Lean's kernel. It borrows an already installed core toolchain and needs no Mathlib. The compiled `Std` imports remain trusted and execution is unisolated, matching the stated honest-source scope. The shared verifier is an external workspace dependency, recorded by path and hash; it is not copied into this repository.

The numerical and formal results are in `results/openai-transfer.json` and `results/periodic-lean.json`. Their checks are included in `src/final_check.py`. The latter's successful exit does not certify the upstream large-coloring theorem.

## Why the earlier negative result still matters

The word `0011` has no safe proper wildcard coset, so it still refutes automatic transfer into the wildcard recursion. Its periodic extension nevertheless has exactly the count above. A failed stronger interface does not rule out a weaker, adequate interface. The new source can be useful through the periodic route without satisfying the extra hypothesis needed by our sharper finite-base recursion.
