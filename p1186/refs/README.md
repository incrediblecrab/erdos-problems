# refs/

Primary sources for Erdős problem #1186. **Nothing in this directory is redistributed** — the PDFs are copyright their authors and publishers. Run `sh refs/fetch.sh` to populate it locally; the repository `.gitignore` keeps everything here out of git except this file and `fetch.sh`.

| Tag | Source | Used in NOTES.md for | Obtainable |
|---|---|---|---|
| [PRS08] | Parrilo, Robertson, Saracino, *JCTA* **115** (2008) 185–192. [DOI](https://doi.org/10.1016/j.jcta.2007.03.006) · [arXiv:math/0609532](https://arxiv.org/abs/math/0609532) | δ₃ ≤ 117/2192 and the 12-block word re-certified in `final_check.py` | yes |
| [LuPe12] | Lu, Peng, *JCTA* **119** (2012) 1048–1065. [DOI](https://doi.org/10.1016/j.jcta.2011.12.004) · [arXiv:1107.2888](https://arxiv.org/abs/1107.2888) | δ̃₄ ≤ 17/300, the B₂₀ bit string, δ̃₄ ≥ 7/192, k=5 and integer k=4,5 bounds | yes |
| [Wo10] | Wolf, *J. Combin.* **1** (2010) 53–68. [DOI](https://doi.org/10.4310/JOC.2010.v1.n1.a4) | δ̃₃ = 1/8, δ̃₄ ≥ 1/32, δ̃₄ ≤ (1−1/259200)/16, the `m_4 = 2·unoriented/p²` convention | **no** — paywalled, no arXiv version |
| [CCS07] | Cameron, Cilleruelo, Serra, *Rev. Mat. Iberoam.* **23** (2007) 385–395. [DOI](https://doi.org/10.4171/RMI/499) | earlier δ̃₄ lower bounds 1/40, 1/33 | author preprint, may move |
| [BCG10] | Butler, Costello, Graham, *Exp. Math.* **19** (2010) 399–411. [DOI](https://doi.org/10.1080/10586458.2010.10390631) | pre-LuPe integer δ₄, δ₅ bounds | **no** — paywalled, no preprint found |
| [GKW24] | Greenwood, Kariv, Williams, *Math. Comp.* (2024). [DOI](https://doi.org/10.1090/mcom/3970) · [arXiv:2301.00336](https://arxiv.org/abs/2301.00336) | PRS08's word is optimal among antisymmetric ≤12-block colourings; the 117/548 normalisation trap | yes |
| [Ve21] | Versteegen. [arXiv:2106.06846](https://arxiv.org/abs/2106.06846) | δ̃₄ < 1/16 strictly in ℤ_p, qualitative only | yes |
| [RS26] | Rué, Spiegel, *FFA* **111** (2026) 102782. [DOI](https://doi.org/10.1016/j.ffa.2025.102782) · [arXiv:2304.00400](https://arxiv.org/abs/2304.00400) | flag algebras, but over 𝔽_pⁿ not ℤ_p | yes |
| [YM26] | Yang, Mao. [arXiv:2604.02115](https://arxiv.org/abs/2604.02115) | Butler–Costello–Graham conjecture; existence, no constant | yes |
| [RL12] | Rabung, Lotts, *Electron. J. Combin.* **19**(2) (2012) #P35. [DOI](https://doi.org/10.37236/2363) | Table 1, W(2,k) for k ≤ 12, reproduced by `src/rabung.c`; the caps from W(5;2) = 178 and W(6;2) = 1132 | yes, open access |
| [HHLM07] | Herwig, Heule, van Lambalgen, van Maaren, *Electron. J. Combin.* **14** (2007) #R6. [DOI](https://doi.org/10.37236/925) | the cyclic zipper behind the f = 2 bases | yes, open access |
| [SS78] | Stevens, Shantaram, *Math. Comp.* **32** (1978) 635–636. [DOI](https://doi.org/10.1090/S0025-5718-1978-0491468-X) | W(5;2) = 178, taken from [RL12, Table 1] | not fetched, not read |
| [KP08] | Kouril, Paul, *Exp. Math.* **17** (2008) 53–61. [DOI](https://doi.org/10.1080/10586458.2008.10129025) | W(6;2) = 1132, taken from [RL12, Table 1] and the title | not fetched, not read |
| [SW17] | Saad, Wolf, *Q. J. Math.* **68** (2017) 125–140. [DOI](https://doi.org/10.1093/qmath/haw011) | the "commonness" thread | abstract only |

`fetch.sh` also saves `erdos1186.html`, the problem statement as posted.

Two sources could not be obtained legitimately ([Wo10], [BCG10]). Every number taken from them is marked in NOTES.md with how it was corroborated — in both cases by two or more independent citing papers that agree on the value. Sci-Hub and similar were deliberately not used.

For UNSAT proof checking, `src/prove_unsat.py` expects a local clone of Marijn Heule's `drat-trim` in `refs/drat-trim/`, which `fetch.sh` clones at the pinned commit and builds, and which is ignored with the rest of `refs/*`. The run recorded in `results/lrat_check.json` used commit `2e3b2dc0ecf938addbd779d42877b6ed69d9a985`, built with `clang drat-trim.c -std=c99 -O2 -o drat-trim` and `clang lrat-check.c -std=c99 -DLONGTYPE -O2 -o lrat-check`.
