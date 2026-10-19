#!/usr/bin/env python3
"""Verify that every VERBATIM block in FCShim.lean is copied character-for-character
from the current formal-conjectures sources.

Exit 0 if every block matches upstream, 1 otherwise.
"""
import os
import re
import sys
import urllib.request

RAW = "https://raw.githubusercontent.com/google-deepmind/formal-conjectures/main/"
SHIM = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Erdos200", "FCShim.lean")

BLOCK = re.compile(
    r"^-- BEGIN VERBATIM: (?P<src>\S+)\n(?P<body>.*?)^-- END VERBATIM$",
    re.MULTILINE | re.DOTALL,
)


def fetch(path: str) -> str:
    with urllib.request.urlopen(RAW + path) as fh:
        return fh.read().decode()


def main() -> int:
    shim = open(SHIM).read()
    blocks = list(BLOCK.finditer(shim))
    if not blocks:
        print("FAIL: no VERBATIM blocks found in", SHIM)
        return 1

    cache: dict[str, str] = {}
    ok = True
    for m in blocks:
        src, body = m.group("src"), m.group("body")
        if src not in cache:
            cache[src] = fetch(src)
        upstream = cache[src]
        nlines = len(body.rstrip("\n").split("\n"))
        if body.rstrip("\n") in upstream:
            print(f"  ok       {nlines:>2} lines <- {src}")
        else:
            ok = False
            print(f"  MISMATCH {nlines:>2} lines <- {src}")
            for line in body.rstrip("\n").split("\n"):
                mark = " " if line in upstream else ">"
                print(f"    {mark} {line}")

    print(("PASS: " if ok else "FAIL: ") + f"{len(blocks)} verbatim block(s) checked")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
