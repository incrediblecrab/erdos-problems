# refs/

Primary sources for every external claim in `NOTES.md`. **Nothing here is committed** — the Roney-Dougal–Tracey preprint is copyright its authors, and the OEIS and GroupNames pages carry their own terms. Run `sh fetch.sh` to populate this directory locally; `.gitignore` excludes the downloads.

| File | Supports |
|---|---|
| `roney-dougal_tracey_2503.05416.pdf` | [RoTr25]. p. 1: "our logarithms are to the base 2" and $\xi=\tfrac16\log 24$; abstract: "settling a conjecture of Pyber from 1993"; Theorem 1 ($2^{n^2/16+\alpha n\log n} \le \|\mathrm{Sub}(S_n)\| \le 2^{n^2/16+\beta n^{3/2}}$); Lemma 2.4(i) and the constant $\zeta_p<4$; Theorems 4, 5, 6; the proof of Theorem 6 on p. 35 with the step "Since $\nu(1-\nu) < 1/16$" |
| `oeis_A005432.txt` | number of subgroups of $S_n$, $n\le18$; the comment giving Pyber's brackets with $c=2^{1/16}$, $d=24^{1/6}$ — an independent witness that the base is 2 |
| `oeis_A000638.txt` | conjugacy classes of subgroups of $S_n$, $n\le20$ |
| `oeis_A000085.txt` | involutions in $S_n$; $a(n,1) = A000085(n)-1$ |
| `oeis_A006116.txt` | Galois numbers $G_r = \sum_k \binom{r}{k}_2$ |
| `arxiv_2503.05416.xml` | arXiv API record: v1 only, submitted March 7, 2025, no journal reference |
| `semanticscholar_citations.json` | citations of 2503.05416 — three at the last fetch, none bearing on the order of a random subgroup (see `NOTES.md` §8); this count grows, so re-fetch before relying on it |
| `oeis_search_2subgroups.txt` | the OEIS search for 1, 2, 4, 20, 76, 631, 3417, 65054, the number of 2-subgroups of $S_n$ for $n \le 8$; it returned "No results." on November 1, 2026 (`NOTES.md` §9.8). A negative search result, so re-run it before relying on it |
| `oeis_search_transitive2.txt`, `oeis_search_transitive2_short.txt` | the OEIS searches for 1, 1, 7, 12945, 20134311622425, the number of transitive 2-subgroups of $S_b$ for $b = 1, 2, 4, 8, 16$, and for its last two terms; both returned "No results." on November 2, 2026 (`NOTES.md` §9.8). Negative results, so re-run them before relying on them |
| `oeis_A218939.txt` | nilpotent subgroups of $S_n$, counting conjugates as distinct, $n \le 13$; `src/nilpotent8.py` checks the 2-subgroup counts for $n \le 8$ against it |
| `erdosproblems_1162_proof_claims.html` | [Sa26]: the proof claim by Amir Sarid on the forum page of problem 1162, submitted October 1, 2026, and its summary ("all but an exponentially small proportion of the subgroups are subdirect products of copies of four permutation groups") |
| `sarid_SPEC_b1f33ddf.md`, `sarid_README_b1f33ddf.md` | [Sa26]: the statement file and README of the linked repository at commit b1f33ddf: the main term $L_n$, theorem T1 and "This statement specifies neither an optimized numerical c nor an effective N0"; the Lean proof is "relative to seven explicitly named published mathematical inputs". Only these two files and the forum page were read |
| `groupnames_T15.html` | [GN]: the orders of the transitive groups 8T1–8T35 and the entry for 8T22 ("Extraspecial group; = D4○D4"), used in `NOTES.md` §9.1 and parsed by `src/final_check.py` |

## Read the base-2 remark from a render, not the text layer

`pdftotext` in its default and `-layout` modes sorts text by position and splices the superscripts of the surrounding displayed math into the sentence, producing

```
... (our logarithms are to the 2 2 base 2). Pyber conjectured ...
```

so a verbatim search against them **fails** even though the sentence is on the page. Two ways to read it correctly:

```sh
pdftotext -raw roney-dougal_tracey_2503.05416.pdf -   # content-stream order
pdftoppm -png -r 150 -f 1 -l 1 -x 180 -y 1275 -W 1300 -H 105 \
         roney-dougal_tracey_2503.05416.pdf quote     # render and read it
```

Both were done. `src/final_check.py` uses `-raw`; the rendered crop was read by eye and reads

> … where ξ = ⅙ log 24 (our logarithms are to the base 2). Pyber conjectured, however, that |Sub(S_n)| = 2^{n²/16+o(n²)}.

Note also that `-raw` loses the minus signs in the exponents of Lemma 2.4's $\zeta_p = \prod_{i\ge1}(1-p^{-i})^{-1}$, rendering it as `(1 − pi)−1`. That formula is therefore *described* in `NOTES.md`, not quoted; the numerical claim $\zeta_2 = 3.46275\ldots < 4$ is checked independently.
