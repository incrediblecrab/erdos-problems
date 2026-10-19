#!/bin/sh
# Re-download the primary sources cited in NOTES.md.
#
# Nothing fetched here is redistributed in this repository: all of it is
# copyright its authors/publishers.  Run this to populate refs/ locally.  The
# repository .gitignore excludes everything in refs/ except this script and
# refs/README.md.
#
# Two of the load-bearing sources are paywalled and cannot be fetched:
#   [Wo10]  Wolf, J. Combin. 1 (2010) 53-68.  No arXiv version exists; the only
#           copy is behind International Press.  DOI 10.4310/JOC.2010.v1.n1.a4
#   [BCG10] Butler-Costello-Graham, Exp. Math. 19 (2010) 399-411.  No preprint
#           was found on any author page.  DOI 10.1080/10586458.2010.10390631
# Their numbers are quoted in NOTES.md with that provenance stated explicitly.
#
# Two more are cited but were not fetched or read; their values are taken from [RL12, Table 1]:
#   [SS78]  Stevens-Shantaram, Math. Comp. 32 (1978) 635-636, W(5;2) = 178.  DOI 10.1090/S0025-5718-1978-0491468-X
#   [KP08]  Kouril-Paul, Exp. Math. 17 (2008) 53-61, W(6;2) = 1132.  DOI 10.1080/10586458.2008.10129025
set -e
cd "$(dirname "$0")"

# [PRS08] Parrilo, Robertson, Saracino, "On the asymptotic minimum number of
# monochromatic 3-term arithmetic progressions", JCTA 115 (2008) 185-192.
# Theorem 5 is the delta_3 <= 117/2192 upper bound; section 4 displays the
# 12-block word that src/final_check.py re-certifies exactly.
curl -sL "https://arxiv.org/pdf/math/0609532" -o prs08.pdf

# [LuPe12] Lu, Peng, "Monochromatic 4-term arithmetic progressions in
# 2-colorings of Z_n", JCTA 119 (2012) 1048-1065.  Theorem 1 is the
# m_4(Z_p) <= 17/150 bound this project beats; B_20 is the bit string that
# src/final_check.py re-derives independently as the m = 20 block optimum.
curl -sL "https://arxiv.org/pdf/1107.2888" -o lupe12.pdf

# [GKW24] Greenwood, Kariv, Williams, "From discrete to continuous:
# monochromatic 3-term arithmetic progressions", Math. Comp. (2024).
# Proves PRS08's word is optimal among antisymmetric <= 12-block colourings.
curl -sL "https://arxiv.org/pdf/2301.00336" -o gkw24.pdf

# [Ve21] Versteegen, "Linear configurations containing 4-term arithmetic
# progressions are uncommon".  Qualitative delta~_4 < 1/16 for Z_p; no constant.
curl -sL "https://arxiv.org/pdf/2106.06846" -o versteegen21.pdf

# [RS26] Rue, Spiegel, "The Rado multiplicity problem in vector spaces over
# finite fields", FFA 111 (2026) 102782.  Flag algebras, but over F_p^n, not Z_p.
curl -sL "https://arxiv.org/pdf/2304.00400" -o rue_spiegel.pdf

# [YM26] Yang, Mao, resolution of the Butler-Costello-Graham conjecture.
# Existence of a beats-random colouring for every rational pattern; no constant.
curl -sL "https://arxiv.org/pdf/2604.02115" -o yang_mao26.pdf

# [CCS07] Cameron, Cilleruelo, Serra, "On monochromatic solutions of equations
# in groups", Rev. Mat. Iberoam. 23 (2007) 385-395.  Author-hosted preprint.
# Note: the typeset Theorem 4.4 reads "4-coloring"; the proof uses two colour
# classes throughout.  See NOTES.md.
curl -sL "https://webspace.maths.qmul.ac.uk/p.j.cameron/preprints/mono.pdf" \
  -o ccs07.pdf || echo "ccs07: author copy moved; see DOI 10.4171/RMI/499"

# [RL12] Rabung, Lotts, "Improving the use of cyclic zippers in finding lower bounds for van der Waerden numbers", Electron. J. Combin. 19(2) (2012) #P35. Table 1 gives W(2,k) for k <= 12 from quadratic-residue colourings, which are the bases of NOTES.md section 4. Open access.
curl -sL "https://www.combinatorics.org/ojs/index.php/eljc/article/download/v19i2p35/pdf/" -o rl12.pdf

# [HHLM07] Herwig, Heule, van Lambalgen, van Maaren, "A new method to construct lower bounds for van der Waerden numbers", Electron. J. Combin. 14 (2007) #R6. The cyclic zipper behind the f = 2 bases. Open access.
curl -sL "https://www.combinatorics.org/ojs/index.php/eljc/article/download/v14i1r6/pdf/" -o hhlm07.pdf

# The problem statement itself.
curl -sL "https://www.erdosproblems.com/1186" -o erdos1186.html

# drat-trim's lrat-check, the external UNSAT proof checker used by src/prove_unsat.py (results/lrat_check.json). Pinned to the commit that run used.
if [ ! -d drat-trim ]; then git clone -q https://github.com/marijnheule/drat-trim drat-trim; fi
git -C drat-trim checkout -q 2e3b2dc0ecf938addbd779d42877b6ed69d9a985
(cd drat-trim && clang drat-trim.c -std=c99 -O2 -o drat-trim && clang lrat-check.c -std=c99 -DLONGTYPE -O2 -o lrat-check)

echo "fetched into $(pwd)"
