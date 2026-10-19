#!/bin/sh
# Re-download the primary sources used in NOTES.md for Erdos problem #1162.
#
# None of these files are redistributed in this repository: the Roney-Dougal
# and Tracey preprint is copyright its authors, and the OEIS pages have their
# own terms.  Run this to populate refs/ locally.  src/final_check.py reads
# the PDF to verify quotations and the GroupNames page for section 9.1, and
# skips those checks if a file is absent.
#
# Everything downloaded here is ignored by .gitignore.
set -e
cd "$(dirname "$0")"

# [RoTr25] Roney-Dougal and Tracey, "Subgroups of symmetric groups:
# enumeration and asymptotic properties", arXiv:2503.05416.  This is the
# primary source for: the base-2 convention (p. 1), Theorem 1, Lemma 2.4,
# and Theorems 4-6 including the proof step that NOTES.md section 7 analyses.
#
# Note the https and the -L: http://export.arxiv.org returns a 301.
curl -sL "https://arxiv.org/pdf/2503.05416v1" -o roney-dougal_tracey_2503.05416.pdf

# The arXiv API record, for the submission date and the version count.
curl -sL "https://export.arxiv.org/api/query?id_list=2503.05416" \
     -o arxiv_2503.05416.xml

# Citations of the preprint, to support the (weak) novelty claim in section 8.
curl -sL "https://api.semanticscholar.org/graph/v1/paper/arXiv:2503.05416/citations?fields=title,year,externalIds&limit=100" \
     -o semanticscholar_citations.json

# OEIS.  A005432 is the number of subgroups of S_n counting conjugates as
# distinct -- its comment states Pyber's brackets with c = 2^(1/16) and
# d = 24^(1/6), which is an independent witness that the logarithm is base 2.
# A000638 is conjugacy classes of subgroups; both are the targets that
# src/subgroups.c is validated against.  A000085 (involutions) validates
# a(n,1) and A006116 (Galois numbers) validates galois_number().
for id in A005432 A000638 A000085 A006116; do
  curl -sL "https://oeis.org/search?q=id:$id&fmt=text" -o "oeis_$id.txt"
done

# [GN] GroupNames, transitive groups of degree up to 15.  NOTES.md section
# 9.1 takes the orders of 8T1-8T35 and the entry for 8T22 from it, and
# src/final_check.py parses the saved page to check them.  (LMFDB blocks
# scripted downloads, so it is not used.)
curl -sL "https://people.maths.bris.ac.uk/~matyd/GroupNames/T15.html" \
     -o groupnames_T15.html

echo "refs/ populated:"
ls -la
