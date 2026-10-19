#!/bin/sh
# Re-download the primary sources cited in NOTES.md into refs/.
#
# Nothing fetched here is redistributed in this repository: the PDFs are copyright their authors and publishers. The repository .gitignore excludes everything in refs/ except this script and refs/README.md.
#
# Three sources the problem page cites were not fetched or read, so NOTES.md takes nothing from them except what the papers below say about them:
#   [PaSh92] Pach, Sharir, JCTA 59 (1992) 12-22.  DOI 10.1016/0097-3165(92)90094-B
#   [Pu74]   Purdy, Discrete Math. 7 (1974) 305-315.  DOI 10.1016/0012-365X(74)90041-7
#   [Ap13]   Apfelbaum, Ph.D. thesis, Tel Aviv University (2013).  No link found.
set -e
cd "$(dirname "$0")"

# [ErPu71] Erdős, Purdy, "Some extremal problems in geometry", JCTA 10 (1971) 246-252, DOI 10.1016/0097-3165(71)90028-8. Scan from the Rényi Institute's Erdős archive. Theorem 2 (p. 249) is the n^2 log log n lattice construction; p. 248 has the Lenz construction and the misprinted conjecture.
curl -sL "https://users.renyi.hu/~p_erdos/1971-20.pdf" -o erpu71.pdf

# [Er75f] Erdős, "On some problems of elementary and combinatorial geometry", Ann. Mat. Pura Appl. (4) 103 (1975) 99-108, DOI 10.1007/BF02414146. Same archive. p. 104 states the problem with "non zero" volume.
curl -sL "https://users.renyi.hu/~p_erdos/1975-25.pdf" -o er75f.pdf

# [DST09] Dumitrescu, Sharir, Tóth, JCTA 116 (2009) 1177-1198, DOI 10.1016/j.jcta.2009.03.008. Theorem 3: the floor(sqrt n) x floor(sqrt n) grid spans (6/pi^2 - o(1)) n^2 minimum-area triangles, the value src/ reproduces first.
curl -sL "https://arxiv.org/pdf/0710.4109v1" -o dst09.pdf

# [ApSh10] Apfelbaum, Sharir, Discrete Comput. Geom. 44 (2010) 753-761, DOI 10.1007/s00454-010-9265-0. O(n^{9/4+eps}); history of the upper bound.
curl -sL "https://arxiv.org/pdf/1001.4764v1" -o apsh10.pdf

# [RaSh17] Raz, Sharir, Combinatorica 37 (2017) 1221-1240, DOI 10.1007/s00493-016-3440-8. Theorem 1: O(n^{20/9}), the current upper bound.
curl -sL "https://arxiv.org/pdf/1501.00379v2" -o rash17.pdf

# [Gr13] Gronwall, "Some asymptotic expressions in the theory of numbers", Trans. AMS 14 (1913) 113-122, DOI 10.1090/S0002-9947-1913-1500940-6. Mertens' product (p. 117) and limsup sigma(x)/(x log log x) = e^C, equation (25) on p. 119. The AMS copy sits behind a Cloudflare check, so this is the JSTOR Early Journal Content copy on the Internet Archive.
curl -sL "https://archive.org/download/jstor-1988773/1988773.pdf" -o gronwall13.pdf

# The problem statement itself. Not hash-checked: the page embeds the access date.
curl -sL "https://www.erdosproblems.com/1086" -o erdos1086.html

# The PDFs read for NOTES.md on September 23, 2026 had these hashes. A mismatch means the source changed since then.
shasum -a 256 -c <<'EOF'
7b4b8d25d052d970a4b794e88d15dfd31639f2081c53a054082d1ea33f0073cb  erpu71.pdf
fc1b5966a19a3e0a0eebbc65ae277037f5d3f818e1b1c45f0fe525564a5d32cf  er75f.pdf
0eac0e9f213c7184878464ce09a768a2119c6d23179faaa777faa771de2b25c6  dst09.pdf
80cfac3e1220e5d29e9f18693e6ca33d00117561a687d0547050e209e1953ff7  apsh10.pdf
2a985dc2bc206bf4978466521173f5082b4a199224ac64b962f56052912c7745  rash17.pdf
80694a34edc5055ca7a7d0c96a469e7a2f47daaa780e8189bebc9878515b1071  gronwall13.pdf
EOF

echo "fetched into $(pwd)"
