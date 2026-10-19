#!/usr/bin/env python3
"""The gate for Erdős #40: check the committed results, the Lean records, and every number, quotation, excerpt and generated table in the documentation.

Each check is a function of a context that holds the documentation, the committed files and, once refs/fetch.sh has run, the downloaded sources. It returns the problems it found, or raises Skip when an input it needs is absent. Before reporting, the script applies each planted defect to a copy of the context and requires the check it targets to fail; every check has at least one.

The greedy list is rechecked here in a third implementation, sharing no code with src/greedy_b2g.c or src/greedy_check.py: a B2[2] test by sorting all pairwise sums a <= b as uint32, the greedy rule for the first RULE_TERMS terms against a table of counts, and the checkpoints recomputed in exact decimal arithmetic. The Theorem B constants are recomputed with mpmath at the recorded minimising n and at sample points, and the table of k in exact rationals.

In the prose of the documentation every numeral must be accounted for: decimals and integers of three or more digits everywhere, formulas included, and smaller integers outside formulas. A numeral is accounted for when it is compared with a fresh value (ANCHORS); lies inside a text that a source backs (FACTS, QUOTES) or a worded claim that the records back (CLAIMS, whose number word is compared with a count); names a problem that refs/fetch.sh downloads, a section of NOTES.md, or a lemma or step that NOTES.md defines; is an exit status, 0 or 1; is the major version of Lean in lean/lean-toolchain or of the Python running this script; is an identifier that refs/fetch.sh fetches; is an ordered-list marker; or lies in one of the PROVENANCE labels, which are printed as unchecked. Small integers inside formulas are mathematics, which the proofs check and this scan does not. Code, links and generated tables are exempt from the scan; the generated tables must equal a fresh rendering instead.

Usage: python src/final_check.py [--write] [--strict] [--plant NAME] [--list]
  --write       regenerate the generated tables in the documentation, then exit
  --strict      count a skipped check as a failure; use it once refs/ is populated
  --plant NAME  apply one planted defect and run every check: exit 1 if some check fails, as it should, 0 if none does
  --list        list the checks and the plants
Exit 0 means no check failed and every plant whose check ran was caught; with --strict, also that no check and no plant was skipped. A run that passes with skips ends "PASS, incomplete".
"""
import argparse
import copy
import hashlib
import html
import json
import math
import re
import shutil
import sqlite3
import subprocess
import sys
import time
import unicodedata
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

import mpmath
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DOCS = ("README.md", "NOTES.md", "lean/README.md", "lean/Erdos40/README.md", "src/README.md", "results/README.md", "refs/README.md", "scripts/README.md")
SUBFOLDER_WORD_LIMIT = 200
FC_PKG = ROOT / "lean" / ".lake" / "packages" / "formal_conjectures"
FC_FILES = {
    "fc_28.lean": "FormalConjectures/ErdosProblems/28.lean",
    "fc_40.lean": "FormalConjectures/ErdosProblems/40.lean",
    "fc_158.lean": "FormalConjectures/ErdosProblems/158.lean",
    "fc_Convolution.lean": "FormalConjecturesForMathlib/Combinatorics/Additive/Convolution.lean",
}
# The tracker database of the collection this repository came from; not part of the repository, compared only when present.
TRACKER_DB = ROOT.parent / "tracker-and-solver" / "data" / "erdos.db"
RULE_TERMS = 1000
PLANT_CASES = ("control", "pin", "fc_edited", "sorry_listed", "sorry_unlisted", "axiom_unlisted", "kernel_off", "statement_lhs", "statement_rhs", "statement_158", "statement_forall_B2", "statement_hijack")
MONTHS = ("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
EXCLUDED_DIRS = {".git", ".lake", "__pycache__"}


class Skip(Exception):
    pass


# ---------------------------------------------------------------- text helpers

def squash(s):
    return re.sub(r"\s+", "", unicodedata.normalize("NFC", s))


def html_text(raw):
    t = re.sub(r"(?is)<(script|style)\b.*?</\1>", " ", raw)
    return html.unescape(re.sub(r"(?s)<[^>]+>", " ", t))


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def blank(text, pattern, flags=0):
    """Replace each match of pattern by spaces of the same length, keeping offsets."""
    return re.sub(pattern, lambda m: re.sub(r"[^\n]", " ", m.group()), text, flags=flags)


GEN_RE = re.compile(r"<!-- generated: (\w+) -->\n(.*?)<!-- end generated -->", re.S)
FENCE_RE = r"(?ms)^```.*?^```[ \t]*$"


def prose(text):
    """The text with code, comments, generated blocks, headings, link targets and URLs blanked."""
    t = blank(text, GEN_RE.pattern, re.S)
    t = blank(t, FENCE_RE)
    t = blank(t, r"(?s)<!--.*?-->")
    t = blank(t, r"`[^`\n]*`")
    t = blank(t, r"(?m)^#{1,6} .*$")
    t = blank(t, r"\]\([^)\s]*\)")
    t = blank(t, r"https?://[^\s)>\]]+")
    return t


def fmt_int(n):
    return f"{n:,}"


def iso_date_words(stamp):
    """2026-09-23T20:18:38-0400 -> September 23, 2026 at 20:18:38 (UTC-04:00)."""
    m = re.fullmatch(r"(\d{4})-(\d\d)-(\d\d)T(\d\d:\d\d:\d\d)([+-])(\d\d):?(\d\d)", stamp)
    if not m:
        raise ValueError(f"unrecognised timestamp {stamp!r}")
    y, mo, d, hms, sign, oh, om = m.groups()
    return f"{MONTHS[int(mo) - 1]} {int(d)}, {y} at {hms} (UTC{'−' if sign == '-' else '+'}{oh}:{om})"


# ---------------------------------------------------------------- context

class Ctx:
    def __init__(self):
        self.docs = {d: (ROOT / d).read_text(encoding="utf-8") for d in DOCS}
        self.files = {}
        self.json = {}
        self.refs = {}
        self.extra_files = []  # planted paths, for the hygiene check

    def file(self, rel):
        if rel not in self.files:
            self.files[rel] = (ROOT / rel).read_bytes()
        return self.files[rel]

    def text(self, rel):
        return self.file(rel).decode("utf-8")

    def js(self, rel):
        if rel not in self.json:
            self.json[rel] = json.loads(self.text(rel))
        return self.json[rel]

    def ref(self, name, mode="squash"):
        """A source in refs/: its text squashed, raw, or for a PDF in pdftotext -layout mode ("layout") or one page squashed ("page:N")."""
        key = (name, mode)
        if key not in self.refs:
            p = ROOT / "refs" / name
            if not p.exists():
                raise Skip(f"refs/{name} missing; run sh refs/fetch.sh")
            if mode == "raw":
                self.refs[key] = p.read_text(encoding="utf-8", errors="replace")
            elif name.endswith(".pdf"):
                if shutil.which("pdftotext") is None:
                    raise Skip("pdftotext not installed")
                args = ["pdftotext"] + (["-layout"] if mode == "layout" else []) + (["-f", mode[5:], "-l", mode[5:]] if mode.startswith("page:") else []) + [str(p), "-"]
                out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
                self.refs[key] = out if mode == "layout" else squash(out)
            else:
                t = p.read_text(encoding="utf-8")
                if name.endswith(".html"):
                    t = html_text(t)
                elif name.endswith(".json"):
                    t = json.dumps(json.loads(t), ensure_ascii=False)
                self.refs[key] = squash(t)
        return self.refs[key]

    def clone(self):
        c = copy.copy(self)
        c.docs = dict(self.docs)
        c.files = dict(self.files)
        c.json = copy.deepcopy(self.json)
        c.refs = dict(self.refs)
        c.extra_files = list(self.extra_files)
        return c

    # derived data
    def greedy_list(self):
        raw = self.text("results/greedy_b2g_g2.txt")
        if not re.fullmatch(r"(?:[1-9][0-9]*\n)+", raw):
            raise ValueError("results/greedy_b2g_g2.txt is not one positive integer per line")
        return [int(x) for x in raw.split()]

    def fc_pin(self):
        m = re.search(r'^rev = "([0-9a-f]{40})"$', self.text("lean/lakefile.toml"), re.M)
        return m.group(1)

    def manifest(self):
        return self.js("lean/lake-manifest.json")["packages"]

    def fetched_problems(self):
        fetch = self.text("refs/fetch.sh")
        m = re.search(r'^for n in ([0-9 ]+); do\n\s*curl -sfL "https://www\.erdosproblems\.com/\$n"', fetch, re.M)
        return {int(x) for x in m.group(1).split()}


def repo_files(ctx):
    """The files a commit would contain: git's view when this is a work tree, otherwise a walk that mirrors .gitignore."""
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "ls-files", "-co", "--exclude-standard"], capture_output=True, text=True, check=True).stdout
        if (ROOT / ".git").exists():
            files = sorted(set(out.split("\n")) - {""})
            return files + ctx.extra_files
    except (OSError, subprocess.CalledProcessError):
        pass
    files = []
    for p in sorted(ROOT.rglob("*")):
        rel = p.relative_to(ROOT).as_posix()
        parts = rel.split("/")
        if any(x in EXCLUDED_DIRS for x in parts) or p.is_dir() or p.name == ".DS_Store":
            continue
        if parts[0] == "refs" and rel not in ("refs/README.md", "refs/fetch.sh"):
            continue
        if rel == "results/fetch.log" or p.suffix == ".pyc":
            continue
        files.append(rel)
    return files + ctx.extra_files


# ---------------------------------------------------------------- registries

# (document, source in refs/, quotation). Each must appear in the document in straight double quotes and, squashed, in the source.
QUOTES = [
    ("NOTES.md", "site_40.html", "This is open, and cannot be resolved with a finite computation."),
    ("NOTES.md", "site_40.html", r"This is a stronger form of the Erdős-Turán conjecture [28] (since establishing this for any function $g(N)\to \infty$ would imply a positive solution to [28])."),
    ("NOTES.md", "site_28.html", r"If $A\subseteq \mathbb{N}$ is such that $A+A$ contains all but finitely many integers then $\limsup 1_A\ast 1_A(n)=\infty$."),
    ("NOTES.md", "site_28.html", r"Another stronger conjecture would be that the hypothesis $\lvert A\cap [1,N]\rvert \gg N^{1/2}$ for all large $N$ suffices."),
    ("NOTES.md", "thread_28.html", r"I'd rather leave it ambiguous whether $0\in \mathbb{N}$ - sometimes this makes sense, sometimes not."),
    ("NOTES.md", "site_39.html", "Erdős proved that for every infinite Sidon set $A$ we have"),
    ("NOTES.md", "site_158.html", "If we replace $2$ by $1$ then $A$ is a Sidon set, for which Erdős proved this is true."),
    ("NOTES.md", "fc_158.lean", "This is proved in [ESS94]."),
    ("NOTES.md", "site_1191.html", "Erdős proved (see for example [HaRo66]) that if $A$ is an infinite Sidon set then"),
    ("NOTES.md", "site_1191.html", "for some constant $c>0$."),
    ("NOTES.md", "pliego_2405.04154v1.pdf", "The above though was proved when g = 1 by Erdős [18, §2 Theorem 8] by showing that in fact any Sidon sequence A satisfies"),
    ("NOTES.md", "pliego_2405.04154v1.pdf", "would in turn follow from the stronger conjectural statement (see Erdős and Fuchs [11])"),
    ("NOTES.md", "pliego_2405.04154v1.pdf", "for any g ≥ 2 then every B2[g] sequence A ⊂ N has the property that"),
    ("NOTES.md", "pliego_2405.04154v1.pdf", "since asymptotic bases of order 2 cannot satisfy (1.4)"),
    ("NOTES.md", "cilleruelo_2008_B2g.pdf", "it is an old open problem to decide whether"),
    ("NOTES.md", "cilleruelo_2008_B2g.pdf", "cannot hold for B2[g] sequences with g ≥ 2 either."),
    ("NOTES.md", "site_39.html", r"Erdős and Rényi have constructed, for any $\epsilon>0$, a set $A$ such that"),
    ("NOTES.md", "site_39.html", r"for all large $N$ and $1_A\ast 1_A(n)\ll_\epsilon 1$ for all $n$."),
    ("NOTES.md", "obryant_math0407117v1.pdf", "The seminal paper of Erdős & Rényi [13] introducing the probabilistic method to combinatorial number theory"),
    ("NOTES.md", "obryant_math0407117v1.pdf", "P. Erdős and A. Rényi, Additive properties of random sequences of positive integers, Acta Arith. 6 (1960), 83–110."),
    ("NOTES.md", "cilleruelo_2008_B2g.pdf", "applying the Borel-Cantelli lemma to conclude that with probability 1, the sequences satisfy the B2[g] property after removing a finite number of elements."),
    ("NOTES.md", "cilleruelo_2008_B2g.pdf", "The exponent 2 + 1/g improves the previous, 2 + 2/g, obtained by Erdős and Renyi in 1960."),
    ("NOTES.md", "pliego_2405.04154v1.pdf", "For any g ≥ 2 there is a B2[g] sequence A ⊂ N satisfying for large x the estimate"),
    ("NOTES.md", "thread_158.html", "I posted a computational note on greedy $B_2[2]$ sets here:"),
    ("NOTES.md", "thread_158.html", "$a_{2000}=7{,}445{,}662$"),
    ("NOTES.md", "thread_158.html", r"stays around $0.7$ in this range (e.g. $\approx 0.733$ at $N=a_{2000}$)"),
    ("NOTES.md", "thread_158.html", r"This supports the possibility that $\liminf_{N\to\infty}|A\cap[1,N]|/\sqrt{N}>0$, i.e. a counterexample to Problem 158 may exist."),
    ("NOTES.md", "zenodo_18452185_v7.pdf", "decreases from 1.8974 at N = 10 to 0.7330 at N = 7,445,662"),
    ("NOTES.md", "zenodo_18452185_v7.pdf", "This behavior is consistent with the possibility that the liminf equals 0, but it is far from conclusive"),
    ("NOTES.md", "oeis_A005282.txt", "Mian-Chowla sequence (a B_2 sequence)"),
]

# (document or "*", text in the document, [(source, needles, absent)]). A needle prefixed raw: is searched in the unnormalised file, and one prefixed page:N: only in page N of a PDF; a source prefixed repo: is read from the repository. The numbers inside the document text are accounted for by the fact.
FACTS = [
    ("NOTES.md", r"offers a prize of \$500", [("site_40.html", ["$500"], [])]),
    ("NOTES.md", "Its discussion thread has no comments.", [("thread_40.html", ["Comments (0)"], [])]),
    ("NOTES.md", "Its sources are [Er95] and [Er97c]", [("site_40.html", ["#40 : [Er95] [Er97c]"], [])]),
    ("NOTES.md", "The files for problems 28, 40 and 158 carry no `formal_proof` attribute", [("fc_28.lean", [], ["formal_proof"]), ("fc_40.lean", [], ["formal_proof"]), ("fc_158.lean", [], ["formal_proof"])]),
    ("NOTES.md", "the site's page for problem 40 says its statement is formalised and links this file on the main branch", [("site_40.html", ["Formalised statement? Yes", "raw:formal-conjectures/blob/main/FormalConjectures/ErdosProblems/40.lean"], [])]),
    ("NOTES.md", "Thomas Bloom wrote in the problem 28 thread", [("thread_28.html", ["I'd be interested to know. Thomas Bloom —"], [])]),
    ("NOTES.md", "the problem, which the page lists as open", [("site_1191.html", ["OPEN This is open"], [])]),
    ("NOTES.md", "his Conjecture 1.1", [("pliego_2405.04154v1.pdf", ["Conjecture 1.1"], [])]),
    ("NOTES.md", "(1.4)", [("pliego_2405.04154v1.pdf", ["(1.4)"], [])]),
    ("NOTES.md", "Pliego's Corollary 1.2", [("pliego_2405.04154v1.pdf", ["Corollary 1.2. For any g ≥ 2 there is a B2 [g] sequence"], [])]),
    ("NOTES.md", "read from a render of p. 3", [("pliego_2405.04154v1.pdf", ["page:3:Corollary 1.2. For any g ≥ 2 there is a B2 [g] sequence"], [])]),
    ("NOTES.md", "Pliego, p. 2:", [("pliego_2405.04154v1.pdf", ["page:2:The above though was proved when g = 1 by Erdős [18, §2 Theorem 8]"], [])]),
    ("NOTES.md", "I read (1.4) from a render of p. 2.", [("pliego_2405.04154v1.pdf", ["page:2:(1.4)", "page:2:since asymptotic bases of order 2 cannot satisfy (1.4)"], [])]),
    ("NOTES.md", "Cilleruelo, p. 1, writes", [("cilleruelo_2008_B2g.pdf", ["page:1:it is an old open problem to decide whether"], [])]),
    ("NOTES.md", "Erdős's paper cited by Pliego as [18]", [("pliego_2405.04154v1.pdf", ["Erdős [18, §2 Theorem 8]"], [])]),
    ("NOTES.md", "note's Table 1", [("zenodo_18452185_v7.pdf", ["Table 1: Greedy B2 [2] set"], [])]),
    ("NOTES.md", "on February 1, 2026, R. Zeraoulia wrote", [("thread_158.html", ["counterexample to Problem 158 may exist. Zeraoulia Rafik — 23:05 on 01 Feb 2026"], [])]),
    ("NOTES.md", "linking a Zenodo record (DOI 10.5281/zenodo.18452185)", [("thread_158.html", ["raw:https://zenodo.org/records/18452185"], []), ("zenodo_18452185.json", ['"doi": "10.5281/zenodo.18452185"'], [])]),
    ("README.md", "the one computation posted on the closely related problem 158", [("site_158.html", ["Comments (1)"], []), ("thread_158.html", ["Zeraoulia Rafik — 23:05 on 01 Feb 2026"], [])]),
    ("NOTES.md", "P. Erdős and A. Rényi (1960)", [("obryant_math0407117v1.pdf", ["P. Erdős and A. Rényi, Additive properties of random sequences of positive integers, Acta Arith. 6 (1960)"], [])]),
    ("README.md", "maintained by Thomas Bloom", [("site_40.html", ["T. F. Bloom, Erdős Problem #40"], []), ("thread_28.html", ["Thomas Bloom —"], [])]),
    ("NOTES.md", "The site is maintained by Thomas Bloom.", [("site_40.html", ["T. F. Bloom, Erdős Problem #40"], []), ("thread_28.html", ["Thomas Bloom —"], [])]),
    ("README.md", "[teorth/erdosproblems](https://github.com/teorth/erdosproblems) is licensed Apache 2.0", [("license_teorth_erdosproblems.txt", ["Apache License", "Version 2.0, January 2004"], [])]),
    ("README.md", "licensed Apache 2.0, which `lean/` imports", [("license_formal_conjectures.txt", ["Apache License", "Version 2.0, January 2004"], [])]),
    ("NOTES.md", "*On the Erdős-Turán Conjecture and the growth of $B_{2}[g]$ sequences*", [("arxiv_2405.04154.xml", ["<title>On the Erdős-Turán Conjecture and the growth of $B_{2}[g]$ sequences</title>"], [])]),
    ("NOTES.md", "(v1, May 7, 2024; the only version, with no journal reference)", [("arxiv_2405.04154.xml", ["raw:<id>http://arxiv.org/abs/2405.04154v1</id>", "raw:<published>2024-05-07T09:46:01Z</published>", "raw:<updated>2024-05-07T09:46:01Z</updated>"], ["raw:journal_ref", "raw:2405.04154v2"])]),
    ("NOTES.md", "*Probabilistic constructions of B2[g] sequences*, preprint dated June 3, 2008, from the author's page at the Universidad Autónoma de Madrid", [("cilleruelo_2008_B2g.pdf", ["Probabilistic constructions of B2[g] sequences", "June 3, 2008", "Universidad Autónoma de Madrid"], [])]),
    ("NOTES.md", "*A Complete Annotated Bibliography of Work Related to Sidon Sequences*", [("obryant_math0407117v1.pdf", ["A Complete Annotated Bibliography of Work Related to Sidon Sequences"], [])]),
    ("NOTES.md", "(v1, July 8, 2004)", [("obryant_math0407117v1.pdf", ["arXiv:math/0407117v1 [math.NT] 8 Jul 2004"], [])]),
    ("NOTES.md", "*Computational Evidence for Erdős Problem #158 via the Greedy B2[2] Construction*, Zenodo, February 1, 2026", [("zenodo_18452185.json", ["Computational Evidence for Erd\u02ddos Problem #158 via the Greedy B2[2] Construction", '"publication_date": "2026-02-01"'], [])]),
    ("NOTES.md", "the file `erdos_problem_158_greedy_computation_note_v7.pdf`", [("zenodo_18452185.json", ['"key": "erdos_problem_158_greedy_computation_note_v7.pdf"'], [])]),
    ("NOTES.md", "*Some of my favourite problems in number theory, combinatorics, and geometry*, Resenhas (1995), 165-186, MR 1370501", [("bib_Er95.html", ["Some of my favourite problems in number theory, combinatorics, and geometry", "Resenhas (1995), 165-186.", "raw:mr=1370501"], [])]),
    ("NOTES.md", "*Some of my favorite problems and results*, The mathematics of Paul Erdős, I (1997), 47-67, MR 1425174", [("bib_Er97c.html", ["Some of my favorite problems and results", "The mathematics of Paul Erdős, I (1997), 47-67.", "raw:mr=1425174"], [])]),
    ("NOTES.md", "Problem C9 of Guy's collection, [Gu04] on the site, which the problem 39 page cites", [("site_39.html", ["problem C9 of Guy's collection [Gu04]"], [])]),
    ("refs/README.md", "The Zenodo record spells Erdős with U+02DD.", [("zenodo_18452185.json", ["Erd\u02ddos Problem #158"], [])]),
    ("refs/README.md", "The note's reference [2] attributes arXiv:math/0407117, O'Bryant's bibliography, to Ruzsa.", [("zenodo_18452185_v7.pdf", ["[2] I. Z. Ruzsa, Erdős and Sidon sets, arXiv:math/0407117"], []), ("obryant_math0407117v1.pdf", ["arXiv:math/0407117v1", "A Complete Annotated Bibliography of Work Related to Sidon Sequences"], [])]),
]

# (document, text just before the number, name of the fresh value). The number that follows must be a correct rounding of the value to the places it is printed with, or equal it if it is an integer.
ANCHORS = [
    ("README.md", "Carried to the ", "terms"),
    ("README.md", "the ratio falls to ", "last_ratio"),
    ("NOTES.md", "the ratio keeps falling, to ", "last_ratio"),
    ("NOTES.md", "0.515 at the ", "terms"),
]

# The session that did the work: the model, the CLI and the date. No script can check them, so the gate prints them. The date has no lasting record: a rerun of scripts/reproduce.sh overwrites the logs, and the site stamps each page with the date of the fetch in its own time zone, which read 2026-09-24 for a fetch at 20:57 on September 23, US Eastern time.
PROVENANCE = ["Claude Opus 5.5", "Copilot CLI 1.0.88", "September 23, 2026"]


def fresh_values(ctx):
    a = ctx.greedy_list()
    with localcontext() as c:
        c.prec = 50
        last_ratio = Decimal(len(a)) / Decimal(a[-1]).sqrt()
    return {"terms": Decimal(len(a)), "last_ratio": last_ratio}


NUMBER_WORDS = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10, "eleven": 11, "twelve": 12}


def phrase_count(phrase):
    """The one number word in a claim's phrase, as an integer."""
    words = [w for w in re.findall(r"[a-z]+", phrase.lower()) if w in NUMBER_WORDS]
    if len(words) != 1:
        raise ValueError(f"the claim {phrase!r} must contain exactly one number word")
    return NUMBER_WORDS[words[0]]


def phrase_numerals(phrase):
    return [int(x) for x in re.findall(r"(?<![\w.])\d+(?![\w.])", phrase)]


def _decreasing(xs):
    return all(x > y for x, y in zip(xs, xs[1:]))


def _increasing(xs):
    return all(x < y for x, y in zip(xs, xs[1:]))


def _greedy(ctx):
    return ctx.js("results/greedy_b2g.json")


def _ratios(ctx):
    return [c["ratio_at_a_k"] for c in _greedy(ctx)["g2"]["checkpoints"]]


def _exponents(ctx):
    return [c["exponent_from_previous"] for c in _greedy(ctx)["g2"]["checkpoints"][1:]]


def _cases(ctx):
    return {c["name"]: c for c in ctx.js("results/lean_plants.json")["cases"]}


def _modules(ctx):
    return re.search(r'^MODULES="([^"]*)"', ctx.text("lean/verify.sh"), re.M).group(1).split()


def _step_count(ctx):
    """The steps verify.sh echoes, if they are numbered from 0 and NOTES.md lists the same numbers after naming the count."""
    steps = re.findall(r'^echo "== (\d+)\. ', ctx.text("lean/verify.sh"), re.M)
    m = re.search(r"`lean/verify\.sh` runs \w+ steps and stops at the first failure\.\n\n((?:\d+\. .*\n)+)", ctx.docs["NOTES.md"])
    listed = re.findall(r"^(\d+)\. ", m.group(1), re.M) if m else []
    return len(steps) if steps == [str(i) for i in range(len(steps))] == listed else None


def _fc_count(ctx):
    fetch = ctx.text("refs/fetch.sh")
    loop = re.search(r'^for n in ([0-9 ]+); do\n\s*curl -sfL "https://raw\.githubusercontent\.com/google-deepmind/formal-conjectures/', fetch, re.M)
    n = (len(loop.group(1).split()) if loop else 0) + len(re.findall(r"-o fc_\w+\.lean$", fetch, re.M))
    return n if n == len(FC_FILES) else None


def note_table_rows(ctx):
    """The rows (N, count, square root, ratio) of the Zenodo note's Table 1, as printed."""
    layout = ctx.ref("zenodo_18452185_v7.pdf", "layout")
    m = re.search(r"Table 1 reports(.*?)Table 1:", layout, re.S)
    return re.findall(r"^\s*([\d,]+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s*$", m.group(1), re.M) if m else []


def _note_places(ctx):
    places = {len(r.partition(".")[2]) for _, _, _, r in note_table_rows(ctx)}
    return places.pop() if len(places) == 1 else None


def _passes_steps(ctx, phrase, names):
    """Each named case passed the steps numbered in the phrase and stopped at the next one."""
    last = max(phrase_numerals(phrase))
    return all(_cases(ctx)[n]["step"] == last + 1 and _cases(ctx)[n]["ok"] for n in names)


# (file, phrase, kind, test, what the test reads): claims made in words. For kind "count" the phrase holds one number word, and test(ctx) must return that number; for kind "holds", test(ctx, phrase) must be true, reading any numeral from the phrase. The numbers inside a phrase are accounted for by its claim.
CLAIMS = [
    ("README.md", "carried to five times as many elements", "count", lambda c: Fraction(_greedy(c)["g2"]["terms"], _greedy(c)["thread_claim"]["k"]), "the terms computed over the thread's k"),
    ("README.md", "Its ratio $A(N)/\\sqrt N$ keeps falling", "holds", lambda c, p: _decreasing(_ratios(c)), "the checkpoint ratios decrease"),
    ("NOTES.md", "Past the range of the comment the ratio keeps falling", "holds", lambda c, p: _decreasing(_ratios(c)), "the checkpoint ratios decrease"),
    ("README.md", "the local growth exponent of the terms is above $2$ at every checkpoint", "holds", lambda c, p: all(e > phrase_numerals(p)[0] for e in _exponents(c)), "every local exponent exceeds the phrase's number"),
    ("NOTES.md", "rises slowly and is above $2$ at every checkpoint", "holds", lambda c, p: _increasing(_exponents(c)) and all(e > phrase_numerals(p)[0] for e in _exponents(c)), "the local exponents increase and exceed the phrase's number"),
    ("NOTES.md", "within each window the smallest ratio is close to the ratio at the window's right end", "holds", lambda c, p: all(x["min_ratio_in_window"] >= 0.97 * x["ratio_at_a_k"] for x in _greedy(c)["g2"]["checkpoints"]), "each window minimum is at least 0.97 times the ratio at a_k"),
    ("NOTES.md", "when lowered to one second", "count", lambda c: _greedy(c)["checks"]["cpu_limit"]["trip_limit_s"] if _greedy(c)["checks"]["cpu_limit"]["ok"] else None, "cpu_limit.trip_limit_s, if the lowered limit stopped the run"),
    ("NOTES.md", "`lean/verify.sh` runs five steps", "count", _step_count, "the steps verify.sh echoes, numbered from 0 as the list in NOTES.md is"),
    ("NOTES.md", "compares six statements", "count", lambda c: len(parse_verify_log(c)["statements"]), "the statements line of results/lean_verify.log"),
    ("NOTES.md", "every constant the three modules add", "count", lambda c: len(_modules(c)), "MODULES in lean/verify.sh"),
    ("NOTES.md", "covers this project's three modules", "count", lambda c: len(_modules(c)), "MODULES in lean/verify.sh"),
    ("NOTES.md", "for the two cases that alter it", "count", lambda c: sum(1 for x in _cases(c).values() if x["expect_step"] == 0), "the cases expected to stop at step 0"),
    ("NOTES.md", "The four statement cases also run the probe alone", "count", lambda c: sum(1 for x in _cases(c).values() if "probe" in x), "the cases with a probe run"),
    ("NOTES.md", "`sorry_unlisted`, `axiom_unlisted` and `statement_hijack` pass steps 1 and 2", "holds", lambda c, p: _passes_steps(c, p, ("sorry_unlisted", "axiom_unlisted", "statement_hijack")), "each case stopped at the step after the phrase's last"),
    ("NOTES.md", "`kernel_off` passes steps 1 to 3 and is caught only by the kernel replay", "holds", lambda c, p: _passes_steps(c, p, ("kernel_off",)), "the case stopped at the step after the phrase's last"),
    ("NOTES.md", "the four formal-conjectures files these notes rely on", "count", _fc_count, "the formal-conjectures files refs/fetch.sh downloads, if they are the ones FC_FILES names"),
    ("NOTES.md", "the note prints its ratios to four places", "count", _note_places, "the decimal places of the ratios in the note's Table 1"),
    ("scripts/reproduce.sh", "on eleven planted defects", "count", lambda c: len(_cases(c)) - 1, "the cases in results/lean_plants.json besides the control"),
]


# ---------------------------------------------------------------- renderers

def render_theorem_b(ctx, doc):
    d = ctx.js("results/theorem_b_check.json")
    n = fmt_int(d["nmax"]).replace(",", "{,}")
    out = [f"Lemma 1 at every $2\\le n\\le {n}$:", "", "| $\\alpha$ | smallest bound / sum | at $n$ | limit as $n\\to\\infty$ |", "|---|---|---|---|"]
    out += [f"| {r['alpha']} | {r['min_rhs_over_lhs']:.4f} | {fmt_int(r['at_n'])} | {r['asymptotic_ratio']:.4f} |" for r in d["sum_bound"]]
    out += ["", f"Lemma 2 at every $3\\le n\\le {n}$:", "", "| $\\epsilon$ | smallest bound / $\\mu_n$ | at $n$ |", "|---|---|---|"]
    out += [f"| {r['epsilon']} | {r['min_rhs_over_mu']:.4f} | {fmt_int(r['at_n'])} |" for r in d["mean_bound"]]
    out += ["", "The choice of $k$, in exact arithmetic:", "", "| $\\epsilon$ | $k$ | $2\\epsilon k$ | eventual bound $2k-1$ on $r_A$ | density exponent $1/2-\\epsilon$ |", "|---|---|---|---|---|"]
    out += [f"| {r['epsilon']} | {r['k']} | {r['two_eps_k']} | {r['eventual_bound_on_sumRep']} | {r['density_exponent']} |" for r in d["k_table"]]
    held = "Every inequality held" if d["all_hold"] else "NOT every inequality held"
    out += ["", f"{held} with the relative float64 margin of {d['float64_margin']:.0e} that the script requires."]
    return "\n".join(out) + "\n"


def render_greedy_table(ctx, doc):
    g = ctx.js("results/greedy_b2g.json")["g2"]
    out = [f"The first {fmt_int(g['terms'])} terms, with candidates up to {fmt_int(g['vmax'])}; the last is {fmt_int(g['last'])}.", "",
           "| $k$ | $a_k$ | $k/\\sqrt{a_k}$ | $a_{k/2}$ | smallest $A(N)/\\sqrt N$ on $[a_{k/2},a_k]$ | local exponent |", "|---|---|---|---|---|---|"]
    for c in g["checkpoints"]:
        e = f"{c['exponent_from_previous']:.4f}" if "exponent_from_previous" in c else ""
        out.append(f"| {fmt_int(c['k'])} | {fmt_int(c['a_k'])} | {c['ratio_at_a_k']:.4f} | {fmt_int(c['window'][0])} | {c['min_ratio_in_window']:.4f} | {e} |")
    return "\n".join(out) + "\n"


def render_note_table(ctx, doc):
    rows = ctx.js("results/greedy_b2g.json")["note_table"]["rows"]
    out = ["| $N$ | $A(N)$, note | $A(N)$, computed | ratio, note | ratio, computed |", "|---|---|---|---|---|"]
    out += [f"| {fmt_int(r['N'])} | {fmt_int(r['count_claimed'])} | {fmt_int(r['count_computed'])} | {r['ratio_claimed']:.4f} | {r['ratio_computed']:.4f} |" for r in rows]
    return "\n".join(out) + "\n"


def parse_verify_log(ctx):
    log = ctx.text("results/lean_verify.log")
    def one(pattern):
        m = re.search(pattern, log, re.M)
        if not m:
            raise ValueError(f"results/lean_verify.log has no line matching {pattern!r}")
        return m
    return {
        "date": one(r"^# date: (\S+)$").group(1),
        "lean": one(r"^# Lean \(version ([0-9.]+),").group(1),
        "platform": one(r"^# platform: (.+)$").group(1),
        "pass": one(r"^PASS: (\d+) theorems on \{propext, Classical\.choice, Quot\.sound\}; statements match FormalConjectures at ([0-9a-f]{40}); probe clean; (.+) replayed by the kernel\.$"),
        "statements": one(r"^statements: (.*) match FormalConjectures$").group(1).split(", "),
        "probe": one(r"^probe: (\d+) constants checked, (\d+) problems$"),
        "replayed": re.findall(r"^replayed (\S+)$", log, re.M),
        "real": one(r"^\s*([0-9.]+) real\b").group(1),
        "rss": int(one(r"^\s*(\d+)\s+maximum resident set size$").group(1)),
        "exit": int(one(r"^# verify\.sh exit (\d+)$").group(1)),
        "hashes": dict((f, h) for h, f in re.findall(r"^# ([0-9a-f]{64}) (\S+)$", log, re.M)),
    }


def render_lean_summary(ctx, doc):
    v = parse_verify_log(ctx)
    stm = ", ".join(f"`{s}`" for s in v["statements"])
    mods = v["replayed"]
    mod_text = ", ".join(f"`{m}`" for m in mods[:-1]) + f" and `{mods[-1]}`"
    return (f"Run on {iso_date_words(v['date'])} with Lean {v['lean']} on {v['platform']}: `verify.sh` exited {v['exit']}. "
            f"All {v['pass'].group(1)} listed theorems depend on no axiom beyond `propext`, `Classical.choice` and `Quot.sound`; "
            f"the statements of {stm} match formal-conjectures at the pin; the probe checked {v['probe'].group(1)} constants and found {v['probe'].group(2)} problems; "
            f"and the kernel replayed {mod_text}. Wall time {v['real']} s; maximum resident set size {fmt_int(v['rss'])} bytes, as `/usr/bin/time -l` reports it.\n")


def render_plants_table(ctx, doc):
    d = ctx.js("results/lean_plants.json")
    cases = d["cases"]
    ok = sum(1 for c in cases if c["ok"])
    lean = re.search(r"version ([0-9.]+)", d["lean"]).group(1)
    stamp = re.sub(r"([+-]\d\d):(\d\d)$", r"\1\2", d["date"])
    out = [f"Run on {iso_date_words(stamp)} with Lean {lean} on {d['platform']}: {ok} of {len(cases)} cases behaved as expected.", "",
           "| case | planted | expected | observed | probe alone |", "|---|---|---|---|---|"]
    def where(rc, step):
        return f"exit {rc}" + (f" at step {step}" if step is not None else "")
    for c in cases:
        probe = f"exit {c['probe']['exit']}" if "probe" in c else ""
        out.append(f"| `{c['name']}` | {c['planted']} | {where(c['expect_exit'], c['expect_step'])} | {where(c['exit'], c['step'])} | {probe} |")
    return "\n".join(out) + "\n"


def render_lean_env(ctx, doc):
    tool = ctx.text("lean/lean-toolchain").strip()
    out = ["| component | revision | requested as |", "|---|---|---|", f"| Lean | `{tool}` | `lean/lean-toolchain` |"]
    for p in ctx.manifest():
        out.append(f"| {p['name']} | `{p['rev']}` | `{p['inputRev']}` |")
    return "\n".join(out) + "\n"


RENDERERS = {"theorem_b": render_theorem_b, "greedy_table": render_greedy_table, "note_table": render_note_table,
             "lean_summary": render_lean_summary, "plants_table": render_plants_table, "lean_env": render_lean_env}
PLACEMENT = {name: "NOTES.md" for name in RENDERERS}


def regenerate(ctx, doc):
    def sub(m):
        name = m.group(1)
        body = RENDERERS[name](ctx, doc) if name in RENDERERS else m.group(2)
        return f"<!-- generated: {name} -->\n{body}<!-- end generated -->"
    return GEN_RE.sub(sub, ctx.docs[doc])


# ---------------------------------------------------------------- checks

def check_statement_docs(ctx):
    problems = []
    notes = ctx.docs["NOTES.md"]
    m = re.search(r"<!-- statement: site_latex_40\.html -->\n```latex\n(.*?)\n```", notes)
    if not m:
        return ["NOTES.md has no statement block after <!-- statement: site_latex_40.html -->"]
    stmt = m.group(1)
    q = re.search(r"^> (For what functions.*)$", ctx.docs["README.md"], re.M)
    want = " ".join(stmt.replace("\\[", " $").replace("\\]", "$ ").split())
    if not q or " ".join(q.group(1).split()) != want:
        problems.append("README.md's quoted question differs from the statement in NOTES.md with \\[ \\] read as $ $")
    if TRACKER_DB.exists():
        con = sqlite3.connect(f"file:{TRACKER_DB}?mode=ro", uri=True)
        row = con.execute("SELECT statement FROM problems WHERE id = 40").fetchone()
        con.close()
        if not row or row[0].strip() != stmt.strip():
            problems.append("the statement in NOTES.md differs from the tracker database's")
    return problems


def check_statement_site(ctx):
    raw = ctx.ref("site_latex_40.html", "raw")
    m = re.search(r'<div id="content" style="white-space: pre-line;">\s*\n(.*?)\n\s*</div>', raw, re.S)
    if not m:
        return ["refs/site_latex_40.html: no content div"]
    site = m.group(1).strip()
    notes = re.search(r"<!-- statement: site_latex_40\.html -->\n```latex\n(.*?)\n```", ctx.docs["NOTES.md"])
    if not notes or notes.group(1) != site:
        return ["the statement in NOTES.md is not byte for byte the content of refs/site_latex_40.html"]
    return []


def check_excerpts(ctx):
    problems, unchecked = [], 0
    notes = ctx.docs["NOTES.md"]
    for name, rel in FC_FILES.items():
        pkg = FC_PKG / rel
        ref = ROOT / "refs" / name
        if pkg.exists() and ref.exists() and pkg.read_bytes() != ctx.file(f"refs/{name}"):
            problems.append(f"refs/{name} differs from the built formal-conjectures {rel}")
    for m in re.finditer(r"<!-- (fc|src): (\S+) -->\n```lean\n(.*?)```", notes, re.S):
        kind, name, block = m.groups()
        # an excerpt may stop before the end of its last line, as a statement does before `:= by`
        block = block.rstrip("\n")
        if kind == "src":
            if block not in ctx.text(name):
                problems.append(f"a Lean block in NOTES.md is not in {name}")
            continue
        if name not in FC_FILES:
            problems.append(f"unknown formal-conjectures excerpt source {name}")
            continue
        texts = []
        if (FC_PKG / FC_FILES[name]).exists():
            texts.append(("the built package", (FC_PKG / FC_FILES[name]).read_text(encoding="utf-8")))
        if (ROOT / "refs" / name).exists():
            texts.append((f"refs/{name}", ctx.text(f"refs/{name}")))
        if not texts:
            unchecked += 1
            continue
        for where, t in texts:
            if block not in t:
                problems.append(f"a Lean excerpt marked {name} is not in {where}: {block.splitlines()[0][:60]}")
    blocks = re.findall(r"```lean\n", notes)
    marked = re.findall(r"<!-- (?:fc|src): \S+ -->\n```lean\n", notes)
    if len(blocks) != len(marked):
        problems.append(f"{len(blocks) - len(marked)} Lean block(s) in NOTES.md carry no source marker")
    if unchecked and not problems:
        raise Skip(f"{unchecked} formal-conjectures excerpt(s) unchecked: neither refs/fc_*.lean nor the built formal-conjectures package is present")
    return problems


def check_quotes(ctx):
    problems = []
    for doc, src, q in QUOTES:
        if f'"{q}"' not in ctx.docs[doc]:
            problems.append(f"{doc} no longer quotes {q[:60]!r}")
        if squash(q) not in ctx.ref(src):
            problems.append(f"refs/{src} does not contain {q[:60]!r}")
    return problems


def check_quote_coverage(ctx):
    problems = []
    registered = {(d, q) for d, _, q in QUOTES}
    for doc, text in ctx.docs.items():
        t = blank(text, GEN_RE.pattern, re.S)
        t = blank(t, FENCE_RE)
        t = blank(t, r"(?s)<!--.*?-->")
        t = blank(t, r"`[^`\n]*`")
        for m in re.finditer(r'"([^"\n]*)"', t):
            if (doc, text[m.start(1):m.end(1)]) not in registered:
                problems.append(f"{doc}:{text.count(chr(10), 0, m.start()) + 1}: quotation not in the registry: {m.group(1)[:60]!r}")
        if re.search(r"[“”]", t):
            problems.append(f"{doc}: curly double quotes, which the registry does not cover")
    return problems


def fact_source_text(ctx, src, needle):
    """(needle, text to search, needle as searched) for one needle of a fact."""
    mode, n = "squash", needle
    if n.startswith("raw:"):
        mode, n = "raw", n[4:]
    elif n.startswith("page:"):
        page, n = n[5:].split(":", 1)
        mode = f"page:{page}"
    if src.startswith("repo:"):
        t = ctx.text(src[5:])
        return n, (t if mode == "raw" else squash(t)), (n if mode == "raw" else squash(n))
    return n, ctx.ref(src, mode), (n if mode == "raw" else squash(n))


def check_facts(ctx):
    problems = []
    for doc, text, sources in FACTS:
        docs = DOCS if doc == "*" else (doc,)
        if not any(text in ctx.docs[d] for d in docs):
            problems.append(f"{doc}: the registered text {text[:60]!r} is not in the document")
        for src, needles, absent in sources:
            for needle in needles:
                n, t, key = fact_source_text(ctx, src, needle)
                if key not in t:
                    problems.append(f"{src} does not contain {n[:70]!r} (for {text[:40]!r})")
            for needle in absent:
                n, t, key = fact_source_text(ctx, src, needle)
                if key in t:
                    problems.append(f"{src} contains {n[:70]!r}, which {text[:40]!r} says it does not")
    return problems


NUM_RE = re.compile(r"(?<![\w.^{/\\])(?:\d{1,3}(?:,\d{3})+(?!\d)|\d+)(?:\.\d+)*(?:st|nd|rd|th)?(?!\w)")


def number_value(token):
    return Decimal(re.sub(r"(st|nd|rd|th)$", "", token).replace(",", ""))


def agrees(token, value):
    printed = number_value(token)
    places = -printed.as_tuple().exponent
    return abs(Decimal(value) - printed) <= Decimal(5) / Decimal(10) ** (places + 1)


def spans(pattern, text):
    return [(m.start(), m.end()) for m in re.finditer(pattern, text)]


def check_numbers(ctx):
    problems = []
    values = fresh_values(ctx)
    fetched = ctx.fetched_problems()
    fetch = ctx.text("refs/fetch.sh")
    notes = ctx.docs["NOTES.md"]
    top = re.findall(r"^## (\d+)\. ", notes, re.M)
    sections = set(top) | set(re.findall(r"^### (\d+\.\d+) ", notes, re.M))
    defined = set(re.findall(r"^\*(Lemma|Step) (\d+)\.?\*", notes, re.M))
    lean_version = re.search(r":v((\d+)\.\d+\.\d+)$", ctx.text("lean/lean-toolchain").strip())
    contents = re.search(r"^\| § \| contents \| label \|\n\|---\|---\|---\|\n((?:\|.*\|\n)+)", notes, re.M)
    listed = re.findall(r"^\| (\d+) \|", contents.group(1), re.M) if contents else []
    if listed != top:
        problems.append(f"NOTES.md: the contents table lists sections {listed}, but the headings are {top}")
    rows = contents.span(1) if contents else (0, 0)
    for doc, text in ctx.docs.items():
        covered = []
        for fdoc, ftext, _ in FACTS:
            if fdoc in ("*", doc):
                covered += spans(re.escape(ftext), text)
        covered += [x for qdoc, _, q in QUOTES if qdoc == doc for x in spans(re.escape(f'"{q}"'), text)]
        covered += [x for rel, phrase, *_ in CLAIMS if rel == doc for x in spans(re.escape(phrase), text)]
        covered += [x for label in PROVENANCE for x in spans(re.escape(label), text)]
        anchored = set()
        for adoc, before, name in ANCHORS:
            if adoc != doc:
                continue
            m = re.search(re.escape(before) + r"(" + NUM_RE.pattern + r")", text)
            if not m:
                problems.append(f"{doc}: anchor {before!r} not followed by a number")
                continue
            anchored.add(m.start(1))
            if not agrees(m.group(1), values[name]):
                problems.append(f"{doc}: {before}{m.group(1)} but the fresh value is {values[name]:.6f}")
        # quotations are checked against their sources by check_quotes; any other straight-quoted span is caught by check_quote_coverage
        t = blank(prose(text), r'"[^"\n]*"')
        outside_math = blank(blank(blank(t, r"\\\$"), r"(?s)\$\$.*?\$\$"), r"\$[^$\n]+\$")
        for m in NUM_RE.finditer(t):
            tok, pos = m.group(), m.start()
            short = "." not in tok and len(re.sub(r"\D", "", tok)) < 3
            if short and outside_math[pos:m.end()] != tok:
                continue  # a small integer in a formula
            if pos in anchored or any(a <= pos < b for a, b in covered):
                continue
            line_start = text.rfind("\n", 0, pos) + 1
            before = text[max(line_start, pos - 80):pos]
            where = f"{doc}:{text.count(chr(10), 0, pos) + 1}"
            if short and not text[line_start:pos].strip() and text.startswith(". ", m.end()):
                continue  # an ordered-list marker
            if doc == "NOTES.md" and rows[0] <= pos < rows[1] and text[line_start:pos] == "| ":
                continue  # the contents table, compared with the headings above
            if before.endswith("§"):
                if tok not in sections:
                    problems.append(f"{where}: §{tok} names no section of NOTES.md")
                continue
            if re.search(r"(?i)(?:\bproblems?\s+|#)(?:\d+(?:,\s*|\s+and\s+))*$", before):
                if not re.fullmatch(r"\d+", tok) or int(tok) not in fetched:
                    problems.append(f"{where}: problem {tok} is not one that refs/fetch.sh downloads")
                continue
            ref = re.search(r"\b(Lemma|Step)s?\s+(?:\d+(?:,\s*|\s+and\s+|\s+to\s+))*$", before)
            if ref:
                if (ref.group(1), tok) not in defined:
                    problems.append(f"{where}: {ref.group(1)} {tok} is not defined in NOTES.md")
                continue
            if re.search(r"\bexit(?:s|ed)?\s+$", before):
                if tok not in ("0", "1"):
                    problems.append(f"{where}: exit status {tok}; the scripts here exit 0 or 1")
                continue
            if re.search(r"\bLean\s+$", before):
                if not lean_version or tok not in (lean_version.group(1), lean_version.group(2)):
                    problems.append(f"{where}: Lean {tok} is not the release in lean/lean-toolchain")
                continue
            if re.search(r"\bPython\s+$", before):
                if tok != str(sys.version_info.major):
                    problems.append(f"{where}: Python {tok}, but this script runs under Python {sys.version_info.major}")
                continue
            if re.search(r"(?:arXiv:|doi:|DOI )$", before):
                full = re.match(r"[\w./]+", text[pos:]).group().rstrip(".")
                if full not in fetch:
                    problems.append(f"{where}: identifier {full} is not one refs/fetch.sh fetches")
                continue
            problems.append(f"{where}: number {tok!r} is not accounted for: ...{text[max(0, pos - 40):m.end() + 10]!r}")
    return problems


def check_hex_versions(ctx):
    problems = []
    pin = ctx.fc_pin()
    revs = {pin} | {p["rev"] for p in ctx.manifest()}
    inputs = {p["inputRev"] for p in ctx.manifest()} | {ctx.text("lean/lean-toolchain").strip().split(":")[-1]}
    recorded = set()
    for rel in ("results/greedy_b2g.json", "results/theorem_b_check.json", "results/lean_plants.json"):
        recorded |= set(re.findall(r"[0-9a-f]{64}", ctx.text(rel)))
    recorded |= set(parse_verify_log(ctx)["hashes"].values())
    fetch = ctx.text("refs/fetch.sh")
    for doc, text in ctx.docs.items():
        for m in re.finditer(r"(?<![0-9a-zA-Z])[0-9a-f]{12,}(?![0-9a-zA-Z])", text):
            h, line = m.group(), text.count("\n", 0, m.start()) + 1
            if not re.search(r"[a-f]", h) or not re.search(r"\d", h):
                continue
            ok = (len(h) == 40 and h in revs) or (len(h) == 64 and h in recorded) or (len(h) == 32 and h in fetch) or (len(h) < 40 and any(r.startswith(h) for r in revs))
            if not ok:
                problems.append(f"{doc}:{line}: hex {h[:16]}… is not the pin, a manifest revision or a recorded hash")
        for m in re.finditer(r"(?<![\w.])v(\d+\.\d+\.\d+)(?![\w.])", text):
            if f"v{m.group(1)}" not in inputs:
                problems.append(f"{doc}:{text.count(chr(10), 0, m.start()) + 1}: version v{m.group(1)} is neither the Lean release nor a manifest input")
        if pin not in text and doc == "NOTES.md":
            problems.append("NOTES.md does not name the formal-conjectures pin")
    return problems


def check_generated(ctx):
    problems = []
    for doc, text in ctx.docs.items():
        names = GEN_RE.findall(text)
        for name, body in names:
            if name not in RENDERERS:
                problems.append(f"{doc}: unknown generated block {name}")
            elif PLACEMENT[name] != doc:
                problems.append(f"{doc}: generated block {name} belongs in {PLACEMENT[name]}")
            elif body != RENDERERS[name](ctx, doc):
                problems.append(f"{doc}: generated block {name} is stale; run src/final_check.py --write")
        if text.count("<!-- generated:") != len(names):
            problems.append(f"{doc}: a generated block is malformed")
    for name, doc in PLACEMENT.items():
        n = sum(1 for x, _ in GEN_RE.findall(ctx.docs[doc]) if x == name)
        if n != 1:
            problems.append(f"{doc} has {n} generated blocks named {name}, expected 1")
    return problems


def check_greedy_hash(ctx):
    problems = []
    g = ctx.js("results/greedy_b2g.json")
    if sha256(ctx.file("results/greedy_b2g_g2.txt")) != g["g2"]["list_sha256"]:
        problems.append("results/greedy_b2g_g2.txt does not match list_sha256")
    if sha256(ctx.file("src/greedy_b2g.c")) != g["source_sha256"]:
        problems.append("src/greedy_b2g.c changed since results/greedy_b2g.json was written")
    if sha256(ctx.file("src/greedy_check.py")) != g["script_sha256"]:
        problems.append("src/greedy_check.py changed since results/greedy_b2g.json was written")
    a = ctx.greedy_list()
    if len(a) != g["g2"]["terms"] or a[-1] != g["g2"]["last"] or a[-1] > g["g2"]["vmax"]:
        problems.append("the list's length, last term or bound disagrees with results/greedy_b2g.json")
    return problems


def check_greedy_b2(ctx):
    a = np.array(ctx.greedy_list(), dtype=np.uint64)
    if a[0] != 1 or not np.all(np.diff(a.astype(np.int64)) > 0):
        return ["the list does not start at 1 and increase strictly"]
    if 2 * int(a[-1]) >= 2 ** 32:
        return ["pair sums would overflow uint32"]
    a32 = a.astype(np.uint32)
    n = len(a32)
    sums = np.empty(n * (n + 1) // 2, dtype=np.uint32)
    pos = 0
    for i in range(n):
        sums[pos:pos + n - i] = a32[i] + a32[i:]
        pos += n - i
    sums.sort(kind="stable")
    triple = sums[2:] == sums[:-2]
    if triple.any():
        n = int(sums[2:][triple][0])
        return [f"{n} has at least three representations a + b with a <= b"]
    return []


def check_greedy_rule(ctx, g=2):
    a = ctx.greedy_list()[:RULE_TERMS]
    R = np.zeros(2 * a[-1] + 2, dtype=np.int32)
    S = np.array([], dtype=np.int64)
    prev = 0
    for k, nxt in enumerate(a):
        cand = np.arange(prev + 1, nxt + 1, dtype=np.int64)
        bad = R[2 * cand] >= g
        if len(S):
            bad |= (R[cand[:, None] + S[None, :]] >= g).any(axis=1)
        if bad[-1]:
            return [f"term {k + 1}, {nxt}, breaks the B2[{g}] condition when added"]
        if not bad[:-1].all():
            m = int(cand[:-1][~bad[:-1]][0])
            return [f"term {k + 1} is {nxt}, but {m} was acceptable first"]
        np.add.at(R, nxt + S, 1)
        R[2 * nxt] += 1
        S = np.append(S, nxt)
        prev = nxt
    return []


def check_greedy_numbers(ctx):
    problems = []
    g = ctx.js("results/greedy_b2g.json")
    a = ctx.greedy_list()
    with localcontext() as c:
        c.prec = 50
        def ratio(k, n):
            return Decimal(k) / Decimal(n).sqrt()
        prev = None
        for cp in g["g2"]["checkpoints"]:
            k = cp["k"]
            lo = k // 2
            dips = min(ratio(j, a[j] - 1) for j in range(lo, k))
            checks = [("a_k", cp["a_k"] == a[k - 1]), ("window", cp["window"] == [a[lo - 1], a[k - 1]]),
                      ("ratio_at_a_k", agrees(f"{cp['ratio_at_a_k']:.6f}", ratio(k, a[k - 1]))), ("min_ratio_in_window", agrees(f"{cp['min_ratio_in_window']:.6f}", dips))]
            if prev is not None:
                e = (Decimal(a[k - 1]) / Decimal(a[prev - 1])).ln() / (Decimal(k) / Decimal(prev)).ln()
                checks.append(("exponent_from_previous", agrees(f"{cp['exponent_from_previous']:.4f}", e)))
            problems += [f"checkpoint k={k}: {name} disagrees with the list" for name, ok in checks if not ok]
            prev = k
        t = g["thread_claim"]
        r = ratio(t["k"], a[t["k"] - 1])
        if t["a_k_computed"] != a[t["k"] - 1] or not agrees(f"{t['ratio_computed']:.6f}", r):
            problems.append("thread_claim: the computed values disagree with the list")
        if t["a_k_claimed"] != a[t["k"] - 1] or not agrees(str(t["ratio_claimed"]), r) or not t["reproduced"]:
            problems.append("thread_claim: the thread's numbers are not reproduced by the list")
        for row in g["note_table"]["rows"]:
            n = row["N"]
            if n > a[-1]:
                problems.append(f"note row N={n} lies beyond the list")
                continue
            cnt = sum(1 for x in a if x <= n)
            r = ratio(cnt, n)
            if row["count_computed"] != cnt or not agrees(f"{row['ratio_computed']:.6f}", r):
                problems.append(f"note row N={n}: the computed values disagree with the list")
            if row["count_claimed"] != cnt or not agrees(f"{row['ratio_claimed']:.4f}", r):
                problems.append(f"note row N={n}: the note's numbers are not reproduced by the list")
    return problems


def check_claimed_numbers_source(ctx):
    problems = []
    g = ctx.js("results/greedy_b2g.json")
    printed = note_table_rows(ctx)
    if not printed:
        return ["the Zenodo note's Table 1 was not found"]
    rows = [(int(n.replace(",", "")), int(c), Decimal(r)) for n, c, _, r in printed]
    want = [(r["N"], r["count_claimed"], Decimal(str(r["ratio_claimed"]))) for r in g["note_table"]["rows"]]
    if rows != want:
        problems.append(f"the note's Table 1 reads {rows}, but results/greedy_b2g.json records {want}")
    thread = ctx.ref("thread_158.html")
    t = g["thread_claim"]
    a_text = f"{t['a_k_claimed']:,}".replace(",", "{,}")
    if squash(f"$a_{{{t['k']}}}={a_text}$") not in thread or squash(f"\\approx {t['ratio_claimed']}$ at $N=a_{{{t['k']}}}$") not in thread:
        problems.append("the thread's a_k or ratio differs from results/greedy_b2g.json")
    return problems


def check_theorem_b(ctx):
    problems = []
    d = ctx.js("results/theorem_b_check.json")
    if sha256(ctx.file("src/theorem_b_check.py")) != d["script_sha256"]:
        problems.append("src/theorem_b_check.py changed since results/theorem_b_check.json was written")
    if not d["all_hold"] or d["float64_margin"] != 1e-10:
        problems.append("theorem_b_check.json does not record every inequality holding with margin 1e-10")
    mpmath.mp.dps = 30
    def S(n, alpha):
        return mpmath.fsum((a * (n - a)) ** (-alpha) for a in range(1, n))
    def mu(n, eps):
        return mpmath.fsum((a * (n - a)) ** (-(mpmath.mpf(1) / 2 + eps)) for a in range(1, (n + 1) // 2))
    for r in d["sum_bound"]:
        al = mpmath.mpf(str(r["alpha"]))
        bound = lambda n: 2 ** (2 * al) * mpmath.mpf(n) ** (1 - 2 * al) / (1 - al)
        limit = 2 ** (2 * al) / ((1 - al) * mpmath.beta(1 - al, 1 - al))
        if abs(limit - r["asymptotic_ratio"]) > 1e-12 * limit:
            problems.append(f"alpha={r['alpha']}: asymptotic ratio {r['asymptotic_ratio']} but mpmath gives {mpmath.nstr(limit, 15)}")
        at = bound(r["at_n"]) / S(r["at_n"], al)
        if abs(at - r["min_rhs_over_lhs"]) > 1e-9 * at:
            problems.append(f"alpha={r['alpha']}: ratio at n={r['at_n']} is {mpmath.nstr(at, 15)}, recorded {r['min_rhs_over_lhs']}")
        for n in (2, 3, 10, 101, 1000, d["nmax"] // 2 + 1):
            if bound(n) / S(n, al) < mpmath.mpf(r["min_rhs_over_lhs"]) * (1 - mpmath.mpf("1e-9")):
                problems.append(f"alpha={r['alpha']}: the ratio at n={n} is below the recorded minimum")
        if not (r["holds"] and r["min_rhs_over_lhs"] > 1):
            problems.append(f"alpha={r['alpha']}: Lemma 1 not recorded as holding")
    for r in d["mean_bound"]:
        e = mpmath.mpf(Fraction(r["epsilon"]).numerator) / Fraction(r["epsilon"]).denominator
        bound = lambda n: 2 ** (2 * e) * mpmath.mpf(n) ** (-2 * e) / (mpmath.mpf(1) / 2 - e)
        at = bound(r["at_n"]) / mu(r["at_n"], e)
        if abs(at - r["min_rhs_over_mu"]) > 1e-9 * at:
            problems.append(f"epsilon={r['epsilon']}: ratio at n={r['at_n']} is {mpmath.nstr(at, 15)}, recorded {r['min_rhs_over_mu']}")
        for n in (3, 4, 11, 100, 1001, d["nmax"] // 2):
            if bound(n) / mu(n, e) < mpmath.mpf(r["min_rhs_over_mu"]) * (1 - mpmath.mpf("1e-9")):
                problems.append(f"epsilon={r['epsilon']}: the ratio at n={n} is below the recorded minimum")
        if not (r["holds"] and r["min_rhs_over_mu"] > 1):
            problems.append(f"epsilon={r['epsilon']}: Lemma 2 not recorded as holding")
    for r in d["k_table"]:
        eps = Fraction(r["epsilon"])
        k = math.floor(Fraction(1) / (2 * eps)) + 1
        want = {"k": k, "two_eps_k": str(2 * eps * k), "summable": 2 * eps * k > 1, "eventual_bound_on_sumRep": 2 * k - 1, "density_exponent": str(Fraction(1, 2) - eps)}
        if {x: r[x] for x in want} != want or not want["summable"]:
            problems.append(f"k table row epsilon={r['epsilon']} is not {want}")
    return problems


def lean_hashes(ctx):
    lean = ROOT / "lean"
    files = sorted(p.relative_to(lean).as_posix() for p in list(lean.glob("*.lean")) + list(lean.glob("Erdos40/*.lean")))
    files += ["verify.sh", "lakefile.toml", "lake-manifest.json", "lean-toolchain"]
    return {f: sha256(ctx.file(f"lean/{f}")) for f in files}


def local_path_problems(name, text):
    home = str(Path.home())
    marks = ["/" + "Users" + "/", "/" + "home" + "/", "/private/var/"]
    found = [m for m in marks + [home] if m in text]
    return [f"{name} contains a local path ({found[0]})"] if found else []


def check_lean_plants(ctx):
    problems = []
    d = ctx.js("results/lean_plants.json")
    want = dict(lean_hashes(ctx), **{"plants.py": sha256(ctx.file("lean/plants.py"))})
    if d["sources_sha256"] != want:
        changed = sorted(k for k in set(want) | set(d["sources_sha256"]) if want.get(k) != d["sources_sha256"].get(k))
        problems.append(f"lean_plants.json was recorded for other sources: {', '.join(changed)}")
    if d["formal_conjectures"] != ctx.fc_pin():
        problems.append("lean_plants.json records another formal-conjectures pin")
    names = [c["name"] for c in d["cases"]]
    if tuple(names) != PLANT_CASES:
        problems.append(f"the plant cases are {names}")
    log = ctx.text("results/lean_plants.log")
    heads = list(re.finditer(r"^### (\S+): (.*)$", log, re.M))
    sections = {}
    for h, nxt in zip(heads, heads[1:] + [None]):
        body = log[h.end():nxt.start() if nxt else len(log)]
        kind = "probe" if h.group(2).startswith("Probe.lean run directly") else "verify"
        sections[(h.group(1), kind)] = (h.group(2), body)
    all_ok = True
    for c in d["cases"]:
        n = c["name"]
        head, body = sections.get((n, "verify"), (None, ""))
        if head != c["planted"]:
            problems.append(f"{n}: no log section headed with what was planted")
            all_ok = False
            continue
        m = re.match(r"\n### verify\.sh exit (\d+) after ", body)
        steps = re.findall(r"^== (\d+)\. ", body, re.M)
        step = int(steps[-1]) if steps and c["expect_exit"] != 0 else None
        found = c["expect_message"] in body
        ok = bool(m) and int(m.group(1)) == c["exit"] == c["expect_exit"] and found == c["message_found"] and found and (c["expect_step"] is None or step == c["step"] == c["expect_step"])
        if c["expect_exit"] == 0:
            ok = ok and c["step"] is None and "PASS:" in body
        if "probe" in c:
            p = c["probe"]
            phead, pbody = sections.get((n, "probe"), ("", ""))
            pm = re.match(r"Probe\.lean run directly, exit (\d+) after ", phead)
            pfound = p["expect_message"] in pbody
            ok = ok and bool(pm) and int(pm.group(1)) == p["exit"] == 1 and pfound and p["message_found"]
        if ok != c["ok"] or not ok:
            problems.append(f"{n}: the log does not support ok={c['ok']}")
        all_ok = all_ok and ok
    if d["all_ok"] != all_ok or not all_ok:
        problems.append("lean_plants.json all_ok is not supported by its cases")
    problems += local_path_problems("results/lean_plants.log", log) + local_path_problems("results/lean_plants.json", ctx.text("results/lean_plants.json"))
    return problems


def check_lean_verify(ctx):
    problems = []
    v = parse_verify_log(ctx)
    verify = ctx.text("lean/verify.sh")
    expected = int(re.search(r"^EXPECTED_AXIOM_LINES=(\d+)$", verify, re.M).group(1))
    modules = re.search(r'^MODULES="([^"]*)"', verify, re.M).group(1).split()
    if v["exit"] != 0:
        problems.append(f"verify.sh exited {v['exit']}")
    if int(v["pass"].group(1)) != expected or v["pass"].group(2) != ctx.fc_pin() or v["pass"].group(3).split() != modules:
        problems.append("the PASS line disagrees with verify.sh or the pin")
    if f"{expected} axiom reports clean; statement check passed" not in ctx.text("results/lean_verify.log"):
        problems.append("the log lacks the clean axiom report count")
    if v["probe"].group(2) != "0":
        problems.append("the probe reported problems")
    if v["replayed"] != modules:
        problems.append(f"the kernel replayed {v['replayed']}, not {modules}")
    if v["hashes"] != lean_hashes(ctx):
        changed = sorted(k for k in set(v["hashes"]) | set(lean_hashes(ctx)) if v["hashes"].get(k) != lean_hashes(ctx).get(k))
        problems.append(f"results/lean_verify.log was recorded for other sources: {', '.join(changed)}")
    return problems + local_path_problems("results/lean_verify.log", ctx.text("results/lean_verify.log"))


def check_claims(ctx):
    problems, skipped = [], []
    for rel, phrase, kind, test, why in CLAIMS:
        text = ctx.docs[rel] if rel in ctx.docs else ctx.text(rel)
        if phrase not in text:
            problems.append(f"{rel}: the claim {phrase!r} is no longer there; update the registry")
            continue
        try:
            if kind == "count":
                want, got = phrase_count(phrase), test(ctx)
                if got != want:
                    problems.append(f"{rel}: {phrase!r} says {want}, but {why} gives {got}")
            elif not test(ctx, phrase):
                problems.append(f"{rel}: {phrase!r} is not supported: {why}")
        except Skip as e:
            skipped.append(f"{phrase[:40]!r}: {e}")
    if skipped and not problems:
        raise Skip("; ".join(skipped))
    return problems


SECRET_RES = [r"ghp_[A-Za-z0-9]{36}", r"github_pat_[A-Za-z0-9_]{20,}", r"sk-[A-Za-z0-9_-]{20,}", r"AKIA[0-9A-Z]{16}", r"-----BEGIN [A-Z ]*PRIVATE KEY-----", r"xox[baprs]-[A-Za-z0-9-]{10,}"]


def check_hygiene(ctx):
    problems = []
    files = repo_files(ctx)
    dirs = {str(Path(f).parent.as_posix()) for f in files} - {"."}
    for d in sorted(dirs | {"."}):
        rel = "README.md" if d == "." else f"{d}/README.md"
        if rel not in files:
            problems.append(f"{d}/ has no README.md")
        elif d != ".":
            text = ctx.docs[rel] if rel in ctx.docs else ctx.text(rel)
            if len(text.split()) >= SUBFOLDER_WORD_LIMIT:
                problems.append(f"{rel} has {len(text.split())} words, not under {SUBFOLDER_WORD_LIMIT}")
    for rel in files:
        if rel not in ctx.docs and rel.endswith("README.md") and rel not in [f"{d}/README.md" for d in dirs] + ["README.md"]:
            problems.append(f"{rel} is not a folder README")
    lic = ctx.text("LICENSE")
    if not lic.startswith("MIT License") or "Permission is hereby granted, free of charge" not in lic:
        problems.append("LICENSE is not the MIT licence")
    ignore = ctx.text(".gitignore").split("\n")
    for rule in ("refs/*", "!refs/README.md", "!refs/fetch.sh", "lean/.lake/", "results/fetch.log", "__pycache__/"):
        if rule not in ignore:
            problems.append(f".gitignore lacks {rule}")
    for rel in files:
        if rel.startswith("refs/") and rel not in ("refs/README.md", "refs/fetch.sh"):
            problems.append(f"{rel} would be committed; refs/ content is not redistributed")
        if rel.lower().endswith(".pdf"):
            problems.append(f"{rel}: a PDF would be committed")
        p = ROOT / rel
        data = ctx.files.get(rel) or (p.read_bytes() if p.exists() else b"")
        if b"\0" in data[:4096]:
            continue
        text = data.decode("utf-8", errors="replace")
        for pat in SECRET_RES:
            if re.search(pat, text):
                problems.append(f"{rel} contains something shaped like a secret ({pat[:12]}…)")
        if rel != "src/final_check.py":
            problems += local_path_problems(rel, text)
    return problems


CHECKS = [
    ("statement_docs", check_statement_docs), ("statement_site", check_statement_site), ("excerpts", check_excerpts),
    ("quotes", check_quotes), ("quote_coverage", check_quote_coverage), ("facts", check_facts), ("numbers", check_numbers),
    ("hex_versions", check_hex_versions), ("generated", check_generated), ("claims", check_claims),
    ("greedy_hash", check_greedy_hash), ("greedy_b2", check_greedy_b2), ("greedy_rule", check_greedy_rule),
    ("greedy_numbers", check_greedy_numbers), ("claimed_numbers_source", check_claimed_numbers_source), ("theorem_b", check_theorem_b),
    ("lean_plants", check_lean_plants), ("lean_verify", check_lean_verify), ("hygiene", check_hygiene),
]


# ---------------------------------------------------------------- plants

def _doc_sub(ctx, doc, old, new):
    if old not in ctx.docs[doc]:
        raise RuntimeError(f"plant target {old!r} not in {doc}")
    ctx.docs[doc] = ctx.docs[doc].replace(old, new, 1)


def _gen_sub(ctx, doc, name, old, new):
    m = next(m for m in GEN_RE.finditer(ctx.docs[doc]) if m.group(1) == name)
    body = m.group(2)
    if old not in body:
        raise RuntimeError(f"plant target {old!r} not in generated block {name}")
    ctx.docs[doc] = ctx.docs[doc][:m.start(2)] + body.replace(old, new, 1) + ctx.docs[doc][m.end(2):]


def _set_list(ctx, a):
    ctx.files["results/greedy_b2g_g2.txt"] = "".join(f"{x}\n" for x in a).encode()


def plant_number(ctx):
    _doc_sub(ctx, "README.md", "the ratio falls to 0.515", "the ratio falls to 0.516")


def plant_unregistered(ctx):
    _doc_sub(ctx, "NOTES.md", "An exploratory run for more terms", "An exploratory run for 20,000 terms")


def plant_short_number(ctx):
    _doc_sub(ctx, "NOTES.md", "Its discussion thread has no comments.", "Its discussion thread has no comments. Its page cites 3 sources.")


def plant_section(ctx):
    _doc_sub(ctx, "README.md", "`NOTES.md` §9, so no claim of novelty is made", "`NOTES.md` §11, so no claim of novelty is made")


def plant_problem_ref(ctx):
    _doc_sub(ctx, "README.md", "Erdős–Turán conjecture, problem 28.", "Erdős–Turán conjecture, problem 29.")


def plant_lemma_ref(ctx):
    _doc_sub(ctx, "NOTES.md", "By Lemmas 2 and 3,", "By Lemmas 2 and 4,")


def plant_generated(ctx):
    _gen_sub(ctx, "NOTES.md", "greedy_table", "| 0.5154 |", "| 0.5155 |")


def plant_quote(ctx):
    global QUOTES
    q = "This is open, and cannot be resolved with a finite computation."
    bad = "This is solved, and cannot be resolved with a finite computation."
    _doc_sub(ctx, "NOTES.md", f'"{q}"', f'"{bad}"')
    QUOTES = [(d, s, bad if x == q else x) for d, s, x in QUOTES]


def plant_unquoted(ctx):
    _doc_sub(ctx, "NOTES.md", "Its discussion thread has no comments.", 'Its discussion thread has no comments. The page also says "this was checked by computer".')


def plant_statement(ctx):
    _doc_sub(ctx, "README.md", r"\rvert \gg \frac{N^{1/2}}{g(N)}$ implies", r"\rvert \ll \frac{N^{1/2}}{g(N)}$ implies")


def plant_statement_site(ctx):
    ctx.ref("site_latex_40.html", "raw")
    key = ("site_latex_40.html", "raw")
    ctx.refs[key] = ctx.refs[key].replace(r"\gg \frac{N^{1/2}}{g(N)}", r"\gg \frac{N^{1/3}}{g(N)}", 1)


def plant_excerpt(ctx):
    _doc_sub(ctx, "NOTES.md", "(fun N : ℕ ↦ √N / g N) =O[atTop]", "(fun N : ℕ ↦ √N / g N) =o[atTop]")


def plant_fact(ctx):
    ctx.ref("site_40.html")
    key = ("site_40.html", "squash")
    ctx.refs[key] = ctx.refs[key].replace("$500", "$1000")


def plant_hex(ctx):
    pin = ctx.fc_pin()
    _doc_sub(ctx, "NOTES.md", f"At the pinned commit `{pin}`", f"At the pinned commit `{pin[:-1]}{'0' if pin[-1] != '0' else '1'}`")


def plant_claim(ctx):
    """The document and the registry changed together: the number word must still be compared with the count."""
    global CLAIMS
    _doc_sub(ctx, "NOTES.md", "`lean/verify.sh` runs five steps", "`lean/verify.sh` runs six steps")
    CLAIMS = [(r, p.replace("runs five steps", "runs six steps"), k, t, w) for r, p, k, t, w in CLAIMS]


def plant_claim_numeral(ctx):
    global CLAIMS
    old = "`kernel_off` passes steps 1 to 3"
    _doc_sub(ctx, "NOTES.md", old, "`kernel_off` passes steps 1 to 4")
    CLAIMS = [(r, p.replace(old, "`kernel_off` passes steps 1 to 4"), k, t, w) for r, p, k, t, w in CLAIMS]


def plant_page(ctx):
    """The document and the registry moved a quotation to the wrong page of its source."""
    global FACTS
    _doc_sub(ctx, "NOTES.md", "Cilleruelo, p. 1, writes", "Cilleruelo, p. 2, writes")
    FACTS = [(d, "Cilleruelo, p. 2, writes", [(s, [n.replace("page:1:", "page:2:") for n in ns], a) for s, ns, a in srcs]) if t == "Cilleruelo, p. 1, writes" else (d, t, srcs) for d, t, srcs in FACTS]


def plant_greedy_hash(ctx):
    g = ctx.js("results/greedy_b2g.json")
    h = g["g2"]["list_sha256"]
    g["g2"]["list_sha256"] = h[:-1] + ("0" if h[-1] != "0" else "1")


def plant_greedy_b2(ctx):
    a = ctx.greedy_list()
    s = set(a)
    reps = {}
    for i in range(60):
        for j in range(i, 60):
            reps[a[i] + a[j]] = reps.get(a[i] + a[j], 0) + 1
    n = next(n for n, r in sorted(reps.items()) if r == 2 and n > 2 * a[30])
    x = next(n - y for y in a if y < n and n - y not in s and n - y > 0 and n - y != y)
    _set_list(ctx, sorted(a + [x]))


def plant_greedy_rule(ctx):
    a = ctx.greedy_list()
    _set_list(ctx, a[:250] + a[251:])


def plant_greedy_numbers(ctx):
    cp = next(c for c in ctx.js("results/greedy_b2g.json")["g2"]["checkpoints"] if c["k"] == 5000)
    cp["ratio_at_a_k"] = round(cp["ratio_at_a_k"] + 0.0001, 6)


def plant_note_source(ctx):
    row = next(r for r in ctx.js("results/greedy_b2g.json")["note_table"]["rows"] if r["N"] == 100000)
    row["count_claimed"] += 1


def plant_theorem_b(ctx):
    row = ctx.js("results/theorem_b_check.json")["k_table"][0]
    row["k"] -= 1


def plant_theorem_b_limit(ctx):
    row = ctx.js("results/theorem_b_check.json")["sum_bound"][2]
    row["asymptotic_ratio"] *= 1 + 1e-9


def plant_lean_plants(ctx):
    c = next(c for c in ctx.js("results/lean_plants.json")["cases"] if c["name"] == "kernel_off")
    c["exit"] = 0


def plant_lean_verify(ctx):
    t = ctx.text("results/lean_verify.log")
    ctx.files["results/lean_verify.log"] = t.replace("replayed Erdos40.Equivalence\n", "", 1).encode()


def plant_hygiene(ctx):
    ctx.docs["src/README.md"] = ctx.docs["src/README.md"] + "\n" + " ".join(["padding"] * SUBFOLDER_WORD_LIMIT) + "\n"


def plant_hygiene_readme(ctx):
    ctx.extra_files.append("data/orphan.json")
    ctx.files["data/orphan.json"] = b"{}\n"


PLANTS = {
    "number": ("numbers", plant_number), "unregistered": ("numbers", plant_unregistered), "short_number": ("numbers", plant_short_number),
    "section": ("numbers", plant_section), "problem_ref": ("numbers", plant_problem_ref), "lemma_ref": ("numbers", plant_lemma_ref), "generated": ("generated", plant_generated),
    "quote": ("quotes", plant_quote), "unquoted": ("quote_coverage", plant_unquoted), "statement": ("statement_docs", plant_statement),
    "statement_site": ("statement_site", plant_statement_site), "excerpt": ("excerpts", plant_excerpt), "fact": ("facts", plant_fact),
    "hex": ("hex_versions", plant_hex), "claim": ("claims", plant_claim), "claim_numeral": ("claims", plant_claim_numeral), "page": ("facts", plant_page),
    "greedy_hash": ("greedy_hash", plant_greedy_hash),
    "greedy_b2": ("greedy_b2", plant_greedy_b2), "greedy_rule": ("greedy_rule", plant_greedy_rule), "greedy_numbers": ("greedy_numbers", plant_greedy_numbers),
    "note_source": ("claimed_numbers_source", plant_note_source), "theorem_b": ("theorem_b", plant_theorem_b), "theorem_b_limit": ("theorem_b", plant_theorem_b_limit),
    "lean_plants": ("lean_plants", plant_lean_plants), "lean_verify": ("lean_verify", plant_lean_verify),
    "hygiene": ("hygiene", plant_hygiene), "hygiene_readme": ("hygiene", plant_hygiene_readme),
}


# ---------------------------------------------------------------- driver

def evaluate(fn, ctx):
    """('pass'|'fail'|'skip'|'error', problems or reason)."""
    try:
        problems = fn(ctx)
    except Skip as e:
        return "skip", [str(e)]
    except Exception as e:  # a check that crashes on its input has found a problem with it
        return "error", [f"{type(e).__name__}: {e}"]
    return ("fail" if problems else "pass"), problems


def with_plant(name, ctx):
    """Apply a plant to a clone and return (clone, restore). A plant that cannot be applied restores the globals and raises."""
    global QUOTES, FACTS, CLAIMS
    saved = (QUOTES, FACTS, CLAIMS)
    def restore():
        global QUOTES, FACTS, CLAIMS
        QUOTES, FACTS, CLAIMS = saved
    c = ctx.clone()
    try:
        PLANTS[name][1](c)
    except BaseException:
        restore()
        raise
    return c, restore


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--plant", choices=sorted(PLANTS))
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for name, _ in CHECKS:
            print(f"check {name}: plants {', '.join(p for p, (t, _) in PLANTS.items() if t == name)}")
        return 0
    ctx = Ctx()
    if args.write:
        for doc in DOCS:
            new = regenerate(ctx, doc)
            if new != ctx.docs[doc]:
                (ROOT / doc).write_text(new, encoding="utf-8")
                print(f"wrote the generated blocks of {doc}")
        return 0

    if args.plant:
        try:
            c, restore = with_plant(args.plant, ctx)
        except Skip as e:
            print(f"plant {args.plant}: cannot be applied here ({e})")
            return 0
        failed = []
        for name, fn in CHECKS:
            status, problems = evaluate(fn, c)
            if status in ("fail", "error"):
                failed.append(name)
                print(f"{name}: {status}: {problems[0]}")
        restore()
        target = PLANTS[args.plant][0]
        if target in failed:
            print(f"plant {args.plant}: caught by {target}")
            return 1
        print(f"plant {args.plant}: NOT caught by {target}" + (f" (other checks failed: {', '.join(failed)})" if failed else ""))
        return 0

    counts = {"pass": 0, "fail": 0, "skip": 0, "error": 0}
    statuses = {}
    t0 = time.time()
    for name, fn in CHECKS:
        t = time.time()
        status, problems = evaluate(fn, ctx)
        counts[status] += 1
        statuses[name] = status
        print(f"{status.upper():5} {name} ({time.time() - t:.1f} s)")
        for p in problems[:20]:
            print(f"      {p}")
        if len(problems) > 20:
            print(f"      … and {len(problems) - 20} more")
    print(f"not checked by any script, printed as written: {'; '.join(PROVENANCE)}")
    if TRACKER_DB.exists():
        print("statement_docs also compared the tracker database")
    else:
        verb = {"pass": "passed", "fail": "failed", "skip": "was skipped", "error": "failed"}[statuses["statement_site"]]
        print(f"statement_docs: no tracker database here, so the only comparison of the statement with a source is statement_site, which {verb}")

    caught, missed, skipped = 0, [], []
    targets = {t for t, _ in PLANTS.values()}
    for name, _ in CHECKS:
        if name not in targets:
            missed.append(f"{name} (no plant)")
    for pname, (target, _) in PLANTS.items():
        try:
            c, restore = with_plant(pname, ctx)
        except Skip as e:
            skipped.append(f"{pname} (cannot be applied here: {e})")
            continue
        except Exception as e:
            missed.append(f"{pname} (the plant could not be applied: {type(e).__name__}: {e})")
            continue
        try:
            status, problems = evaluate(dict(CHECKS)[target], c)
        finally:
            restore()
        if status in ("fail", "error"):
            caught += 1
        elif status == "skip":
            skipped.append(f"{pname} (target {target} skipped)")
        else:
            missed.append(f"{pname} (target {target} passed)")
    for m in skipped:
        print(f"plant skipped: {m}")
    for m in missed:
        print(f"plant not caught: {m}")
    ok = counts["fail"] == 0 and counts["error"] == 0 and not missed and (not args.strict or (counts["skip"] == 0 and not skipped))
    print(f"{counts['pass']} checks passed, {counts['fail'] + counts['error']} failed, {counts['skip']} skipped; {caught} of {len(PLANTS)} plants caught, {len(skipped)} skipped; {time.time() - t0:.0f} s; {('PASS, incomplete' if counts['skip'] or skipped else 'PASS') if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
