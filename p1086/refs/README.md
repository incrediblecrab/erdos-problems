# refs/

**Objective.** Record the primary sources behind `NOTES.md`, which cites them by tag. **Inputs.** `sh refs/fetch.sh` downloads them here, checks each PDF against its SHA-256 from September 23, 2026, and carries the full citations. Nothing is redistributed: `.gitignore` excludes everything here but this file and `fetch.sh`.

| tag | link | used for |
|---|---|---|
| [ErPu71] Erdős, Purdy | [DOI](https://doi.org/10.1016/0097-3165(71)90028-8) | the $n^2\log\log n$ lower bound (Thm 2, p. 249), $O(n^{5/2})$ (Thm 1), the p. 248 misprint |
| [Er75f] Erdős | [DOI](https://doi.org/10.1007/BF02414146) | "non zero" area (p. 104) |
| [DST09] Dumitrescu, Sharir, Tóth | [arXiv:0710.4109](https://arxiv.org/abs/0710.4109) | the grid count reproduced first (Thm 3); $O(n^{44/19})$ |
| [ApSh10] Apfelbaum, Sharir | [arXiv:1001.4764](https://arxiv.org/abs/1001.4764) | $O(n^{9/4+\varepsilon})$ |
| [RaSh17] Raz, Sharir | [arXiv:1501.00379](https://arxiv.org/abs/1501.00379) | $O(n^{20/9})$ (Thm 1), the current upper bound |
| [Gr13] Gronwall | [DOI](https://doi.org/10.1090/S0002-9947-1913-1500940-6) | Mertens' product (p. 117) and $\limsup\sigma(N)/(N\log\log N)=e^\gamma$ (eq. 25) |

Not read, so their bounds are as later papers report them: [PaSh92], [Pu74], and [Ap13], for which no copy was found. `fetch.sh` also saves the problem page as `erdos1086.html`, unhashed because the page embeds the access date.
