# refs/

Primary sources for every external claim in `NOTES.md`. **Nothing here is committed** — the Monthly article, the Green–Tao note and the five arXiv papers of NOTES §7.1 are copyright their authors and publishers, and the OEIS and Luhn pages carry their own terms. Run `sh fetch.sh` to populate this directory locally; `.gitignore` excludes the downloads.

| File | Supports |
|---|---|
| `oeis_A005115.txt` | $a(k)$, the least last term of a $k$-AP of primes, $k\le26$ |
| `oeis_A113827.txt`, `oeis_A093364.txt`, `oeis_A133277.txt` | first terms, common differences, full progressions |
| `b005115.txt` | A005115 b-file |
| `luhn.html` | OEIS-independent record table; minimality proofs for $a(23)$–$a(26)$ (Perrenet, Petukhov); the $a(27)$ upper bound |
| `granville_PrimePatterns.pdf` | eq. (2.1) $\bigl(e^{1-\gamma}k/2\bigr)^{k/2}$; the quotient range $(2/5,2)$ for $15\le n\le21$; the tower rendered with **eight** 2's |
| `green_tao_envelope.pdf` | "an exponential tower in $k$ of height **seven**"; final display with 7 arrows; "$100k$ will certainly suffice" |
| `erdos_readme.md` | problem 200 recorded as open |
| `ErdosProblem200.lean` | the Lean statement — both declarations are `sorry` |
| `tao_teravainen_2107.02158.pdf` | arXiv v4. NOTES §7.1: Theorem 1.6 and Example 1.7 (p. 5), "effective" (p. 1), Theorem 6.2 (p. 28), Theorems 8.1 and 8.3 (pp. 41, 44), the §9 sketch and footnote 13 (pp. 49–50), Theorem A.3 and its $C_d$ (p. 52) |
| `manners_1811.00718.pdf` | arXiv v2. Theorem 1.1.2 (p. 2), with constants $O_s(1)$; read from a render, since the text layer flattens $\exp\exp(O_s(\delta)^{O_s(1)})$ |
| `leng_2312.10772.pdf` | arXiv v5. Theorem 5 (p. 6), the log-power count, ineffective by Siegel's theorem; Theorem 7 (p. 32), the effective Gowers-norm bound |
| `lss_2402.17994.pdf` | arXiv v3. Theorem 1.2 (p. 1), constants $O_s(1)$; its restatement of Manners' bound as $\exp\exp(\delta^{-O_s(1)})$ (p. 4) |
| `teravainen_wang_2603.09281.pdf` | arXiv v1. Theorem 1.1 (p. 2), "for some $c_k>0$" |

Two of these must be read from a **render, not the PDF text layer**, which mangles stacked exponents:

```sh
pdftoppm -png -r 400 -f 3 -l 3 granville_PrimePatterns.pdf gran
```

The Granville text layer renders eq. (2.1) as `(e1<0x00>γk/2)k/2` and the tower as the flat string `22222222100k`; `green_tao_envelope.pdf` has a broken font encoding in which `(` is the glyph `2`. The readings in `NOTES.md` were taken from 400 dpi page images.
