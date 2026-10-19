/-
Copyright (c) 2026 Max Marquardt.
Released under the MIT licence; see LICENSE at the repository root.

Erdős problem #200, upper bound, built against `PrimeNumberTheoremAnd`.

The sibling project `../lean` proves the bound with constant `1` under an explicit
hypothesis `PNT`, because Mathlib has no Prime Number Theorem.  This project discharges
that hypothesis using `PrimeNumberTheoremAnd.chebyshev_asymptotic`, giving

  `ErdosProblem200.erdos_200_variants_upper_unconditional`

with no hypothesis: the statement `formal-conjectures` leaves as `sorry` at
`Erdos200.erdos_200.variants.upper`.

Nothing here addresses the open question, which is whether the bound can be improved to
`o(log N)`.
-/
import Erdos200.Structure
import Erdos200.Bound
import Erdos200.FCShim
import Erdos200.Asymptotic
import Erdos200.Unconditional
import Erdos200.Audit
