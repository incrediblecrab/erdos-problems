#!/bin/sh
# Download the sources that NOTES.md and README.md quote, for Erdős problem 40. None is redistributed here: the site pages belong to erdosproblems.com, the formal-conjectures files are Apache-2.0 (© Google LLC) and fetched from the pinned commit, the three papers are copyright their authors, and the Zenodo note is CC BY 4.0. .gitignore excludes everything this writes. src/final_check.py checks every quotation against these files and reports a quotation whose source is missing as skipped, not passed.
set -e
cd "$(dirname "$0")"

# GitHub's REST API allows few unauthenticated requests an hour. Its four files below are fetched with one retry and, if both attempts fail, a warning instead of an exit: without them src/final_check.py reports the facts that need them as skipped, not passed.
gh_api() {
  curl -sfL "$1" -o "$2" && return 0
  sleep 10
  curl -sfL "$1" -o "$2" && return 0
  rm -f "$2"
  echo "warning: GitHub API request failed twice, probably rate-limited; $2 not fetched" >&2
}

# The formal-conjectures commit is defined once, in lean/lakefile.toml.
FC_PIN=$(sed -n 's/^rev = "\([0-9a-f]\{40\}\)"$/\1/p' ../lean/lakefile.toml)
[ -n "$FC_PIN" ] || { echo "no 40-hex rev in lean/lakefile.toml" >&2; exit 1; }

# erdosproblems.com pages. The quotations were checked against the versions of September 23, 2026 and rechecked against those of November 2, 2026; the site edits its commentary, so a later fetch can differ. Since October 6, 2026 the pages show no problem status (blog_9.html below), so the gate takes statuses from teorth_problems.yaml instead.
for n in 28 39 40 158 1191; do
  curl -sfL "https://www.erdosproblems.com/$n" -o "site_$n.html"
done

# The site owner's blog post of October 6, 2026 announcing, among other changes, that the site no longer displays problem statuses.
curl -sfL "https://www.erdosproblems.com/forum/thread/blog:9" -o blog_9.html

# The LaTeX source of problem 40, against which the statement in README.md is diffed.
curl -sfL "https://www.erdosproblems.com/latex/40" -o site_latex_40.html

# The site's bibliography entries for the problem's two sources, which its page loads from /bibs/: the titles, pages and MathSciNet numbers in NOTES.md §9. Neither source itself was read.
for k in Er95 Er97c; do
  curl -sfL "https://www.erdosproblems.com/bibs/$k" -o "bib_$k.html"
done

# Discussion threads: 40 has no comments; 28 has the remark on whether 0 is in ℕ; 158 has the greedy B2[2] computation that src/greedy_check.py reproduces.
for n in 28 40 158; do
  curl -sfL "https://www.erdosproblems.com/forum/thread/$n" -o "thread_$n.html"
done

# formal-conjectures statements at the commit that lean/lakefile.toml pins, and the file defining sumRep, the representation count they use.
for n in 28 40 158; do
  curl -sfL "https://raw.githubusercontent.com/google-deepmind/formal-conjectures/$FC_PIN/FormalConjectures/ErdosProblems/$n.lean" -o "fc_$n.lean"
done
curl -sfL "https://raw.githubusercontent.com/google-deepmind/formal-conjectures/$FC_PIN/FormalConjecturesForMathlib/Combinatorics/Additive/Convolution.lean" -o fc_Convolution.lean

# formal-conjectures' main branch as read on November 2, 2026, not the pin, for §1.2 and §3: the four files above at that commit, the file for problem 1191, and the commit's lake-manifest.json; GitHub's comparison of the pin with that commit; and GitHub's list of the commits up to that one that touched the files of problems 40, 158 and 1191, with their dates and pull-request numbers. Nothing in lean/ builds against this commit.
FC_HEAD=7bdc4a42f24c0316f81c73b26121794b2bed96a3
for f in 28 40 158 1191; do
  curl -sfL "https://raw.githubusercontent.com/google-deepmind/formal-conjectures/$FC_HEAD/FormalConjectures/ErdosProblems/$f.lean" -o "fchead_$f.lean"
done
curl -sfL "https://raw.githubusercontent.com/google-deepmind/formal-conjectures/$FC_HEAD/FormalConjecturesForMathlib/Combinatorics/Additive/Convolution.lean" -o fchead_Convolution.lean
curl -sfL "https://raw.githubusercontent.com/google-deepmind/formal-conjectures/$FC_HEAD/lake-manifest.json" -o fchead_lake-manifest.json
gh_api "https://api.github.com/repos/google-deepmind/formal-conjectures/compare/$FC_PIN...$FC_HEAD?per_page=1" fchead_compare.json
for f in 40 158 1191; do
  gh_api "https://api.github.com/repos/google-deepmind/formal-conjectures/commits?path=FormalConjectures/ErdosProblems/$f.lean&sha=$FC_HEAD&per_page=100" "fchead_log_$f.json"
done

# The problem statuses, from teorth/erdosproblems, the site's companion repository: status.state of each problem in data/problems.yaml, read from the main branch.
curl -sfL "https://raw.githubusercontent.com/teorth/erdosproblems/main/data/problems.yaml" -o teorth_problems.yaml

# The licences named in the attribution section of README.md: that of teorth/erdosproblems, the site's companion repository, and that of formal-conjectures at the pin.
curl -sfL "https://raw.githubusercontent.com/teorth/erdosproblems/main/LICENSE" -o license_teorth_erdosproblems.txt
curl -sfL "https://raw.githubusercontent.com/google-deepmind/formal-conjectures/$FC_PIN/LICENSE" -o license_formal_conjectures.txt

# J. Pliego, "On the Erdős-Turán Conjecture and the growth of B2[g] sequences", arXiv:2405.04154v1, May 7, 2024: statement (1.4) and its attribution to Erdős and Fuchs, the case g = 1, and Corollary 1.2. The API record shows whether a later version exists.
curl -sfL "https://arxiv.org/pdf/2405.04154v1" -o pliego_2405.04154v1.pdf
curl -sfL "https://export.arxiv.org/api/query?id_list=2405.04154" -o arxiv_2405.04154.xml

# J. Cilleruelo, "Probabilistic constructions of B2[g] sequences", preprint dated June 3, 2008, from the author's page: its account of the Erdős–Rényi method and exponent, and its remark that the case g ≥ 2 is open.
curl -sfL "https://matematicas.uam.es/~franciscojavier.cilleruelo/Preprints/probabilistic%20constructions%20of%20B2g.pdf" -o cilleruelo_2008_B2g.pdf

# OEIS A005282, the Mian-Chowla sequence (the greedy B2[1] set), as the known answer for src/greedy_check.py: the entry, whose %N line gives the definition, and the b-file of terms. OEIS content is CC BY-SA 4.0.
curl -sfL "https://oeis.org/search?q=id:A005282&fmt=text" -o oeis_A005282.txt
curl -sfL "https://oeis.org/A005282/b005282.txt" -o b005282.txt

# K. O'Bryant, "A Complete Annotated Bibliography of Work Related to Sidon Sequences", arXiv:math/0407117v1, July 8, 2004: its title, which shows that the arXiv identifier the Zenodo note gives for a paper of Ruzsa is this bibliography, and section 3.6 on Erdős and Rényi.
curl -sfL "https://arxiv.org/pdf/math/0407117v1" -o obryant_math0407117v1.pdf
# Its LaTeX source, for the MathSciNet reviews it carries: its \mathreview macro prints nothing, so the reviews of Garrison, Zhang and Kiss that §3 and §4.2 use are in the source only. The e-print is a gzipped tar whose Sidon.arXiv.tex is ISO-8859-1; it is stored as UTF-8.
curl -sfL "https://arxiv.org/e-print/math/0407117v1" -o obryant_math0407117v1.tar.gz
tar -xzOf obryant_math0407117v1.tar.gz Sidon.arXiv.tex > obryant_math0407117v1.latin1.tex
iconv -f ISO-8859-1 -t UTF-8 obryant_math0407117v1.latin1.tex > obryant_math0407117v1.tex

# R. Zeraoulia, "Computational Evidence for Erdős Problem #158 via the Greedy B2[2] Construction", Zenodo, February 1, 2026, DOI 10.5281/zenodo.18452185, the note linked from the problem 158 thread: its Table 1, which src/greedy_check.py reproduces. The record lists one file; Zenodo's MD5 for it is checked, so a replaced file is noticed.
curl -sfL -H "Accept: application/json" "https://zenodo.org/api/records/18452185" -o zenodo_18452185.json
curl -sfL "https://zenodo.org/api/records/18452185/files/erdos_problem_158_greedy_computation_note_v7.pdf/content" -o zenodo_18452185_v7.pdf
[ "$(md5 -q zenodo_18452185_v7.pdf 2>/dev/null || md5sum zenodo_18452185_v7.pdf | cut -d' ' -f1)" = 034aec0a94f5388acfd421ec3b18fa76 ] || { echo "zenodo_18452185_v7.pdf does not match the MD5 in its Zenodo record" >&2; exit 1; }

echo "refs/ populated:"
for f in site_*.html blog_9.html bib_*.html thread_*.html fc_*.lean fchead_*.lean fchead_lake-manifest.json fchead_compare.json fchead_log_*.json teorth_problems.yaml license_*.txt pliego_2405.04154v1.pdf arxiv_2405.04154.xml cilleruelo_2008_B2g.pdf oeis_A005282.txt b005282.txt obryant_math0407117v1.pdf obryant_math0407117v1.tex zenodo_18452185.json zenodo_18452185_v7.pdf; do
  if [ -e "$f" ]; then shasum -a 256 "$f"; else echo "missing: $f"; fi
done
