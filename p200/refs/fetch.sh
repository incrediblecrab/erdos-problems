#!/bin/sh
# Re-download the primary sources used in NOTES.md.
#
# None of these files are redistributed in this repository: the Monthly article and
# the Green-Tao note are copyright their authors/publishers, and the OEIS and Luhn
# pages have their own terms.  Run this to populate refs/ locally.  verify_records.py
# reads refs/luhn.html and refs/oeis_A005115.txt, so run this before it.
#
# Everything downloaded here is ignored by .gitignore.
set -e
cd "$(dirname "$0")"

# OEIS -- A005115 is the sequence a(k) = least last term of a k-AP of primes.
for id in A005115 A113827 A093364 A133277; do
  curl -sL "https://oeis.org/search?q=id:$id&fmt=text" -o "oeis_$id.txt"
done
curl -sL "https://oeis.org/A005115/b005115.txt" -o b005115.txt

# Norman Luhn's AP-of-primes records -- an OEIS-independent source, and the source
# for the minimality proofs of a(23)..a(26) and for the a(27) upper bound.
curl -sL "https://www.pzktupel.de/PAP/aprecords.php" -o luhn.html

# Granville, "Prime Number Patterns", Amer. Math. Monthly 115 (2008) 279-296.
# Section 2.1 has eq. (2.1) and the Green-Tao tower.  Read it from a render, not the
# text layer: pdftoppm -png -r 400 -f 3 -l 3 granville_PrimePatterns.pdf out
curl -sL "https://dms.umontreal.ca/~andrew/PDF/PrimePatterns.pdf" -o granville_PrimePatterns.pdf

# Green & Tao, "A bound for progressions of length k in the primes" (expository note).
# The abstract states the tower height ("seven"); the final display on p.3 has 7 arrows.
curl -sL "https://people.maths.ox.ac.uk/greenbj/papers/back-of-an-envelope.pdf" -o green_tao_envelope.pdf

# Problem status: the teorth/erdosproblems database, and the Lean formalisation.
curl -sL "https://raw.githubusercontent.com/teorth/erdosproblems/main/README.md" -o erdos_readme.md
curl -sL "https://raw.githubusercontent.com/google-deepmind/formal-conjectures/main/FormalConjectures/ErdosProblems/200.lean" -o ErdosProblem200.lean

echo "done; see refs/README.md for what each file supports"
