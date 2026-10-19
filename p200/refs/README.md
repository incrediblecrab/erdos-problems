# refs/

Primary sources for every external claim in `NOTES.md`. **Nothing here is committed** — the Monthly article and the Green–Tao note are copyright their authors and publishers, and the OEIS and Luhn pages carry their own terms. Run `sh fetch.sh` to populate this directory locally; `.gitignore` excludes the downloads.

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

Two of these must be read from a **render, not the PDF text layer**, which mangles stacked exponents:

```sh
pdftoppm -png -r 400 -f 3 -l 3 granville_PrimePatterns.pdf gran
```

The Granville text layer renders eq. (2.1) as `(e1<0x00>γk/2)k/2` and the tower as the flat string `22222222100k`; `green_tao_envelope.pdf` has a broken font encoding in which `(` is the glyph `2`. The readings in `NOTES.md` were taken from 400 dpi page images.
