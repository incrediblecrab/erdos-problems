import Lean

/-!
# Axiom and statement probe, run from a separate process

`verify.sh` runs this with `lake env lean --run Probe.lean`. The probe imports only `Lean`, so this code is elaborated without the modules it audits and loads them at runtime with `importModules`, which does not run their initializers. No instance, macro or elaborator they define can change how this code runs. `Erdos40/Audit.lean` makes the same statement checks during the build, but it runs inside the environment it audits.

Axioms: for every constant that `Erdos40.Basic`, `Erdos40.Equivalence` or `Erdos40.Problem158` adds, auxiliary ones included, it collects the axioms the constant depends on and flags any outside `propext`, `Classical.choice`, `Quot.sound`.

Statements: it compares elaborated types, as `Expr`s, with the declarations in `FormalConjectures`, making the same comparisons as `Erdos40/Audit.lean`. `strong` below is `Erdos40.Erdos40For (fun _ : ℕ ↦ (1 : ℝ))` written out by hand, in the form Lean elaborates it to.

Exit 0 means nothing was flagged.
-/

open Lean

def targets : Array Name := #[`Erdos40.Basic, `Erdos40.Equivalence, `Erdos40.Problem158]

def standardAxioms : Array Name := #[`propext, `Classical.choice, `Quot.sound]

def required : Array Name := #[`ErdosProblem40.exists_iff, `ErdosProblem40.answerSet_nonempty_iff, `ErdosProblem40.strong_implies_erdos_28, `ErdosProblem40.uniform_dip, `ErdosProblem40.Erdos40For.of_isBigO, `ErdosProblem40.not_erdos40For_sqrt, `ErdosProblem40.not_erdos40For_of_sqrt_isBigO, `ErdosProblem40.univ_hypothesis_and_conclusion, `ErdosProblem40.sumRep_le_of_B2, `ErdosProblem40.B2_of_sumRep_le, `ErdosProblem40.count_sq_le, `ErdosProblem40.liminf_eq_zero_of_bound, `ErdosProblem40.strong_iff_forall_B2, `ErdosProblem40.weaker_iff_forall_B2, `ErdosProblem40.strong_implies_erdos_158]

/-- `(1 : ℝ)` as Mathlib elaborates it. -/
def realOne : Expr :=
  mkApp3 (mkConst `OfNat.ofNat [Level.zero]) (mkConst `Real) (mkRawNatLit 1)
    (mkApp2 (mkConst `One.toOfNat1 [Level.zero]) (mkConst `Real) (mkConst `Real.instOne))

/-- `Erdos40.Erdos40For (fun _ : ℕ ↦ (1 : ℝ))`. -/
def strong : Expr := mkApp (mkConst `Erdos40.Erdos40For) (.lam `x (mkConst `Nat) realOne .default)

/-- The problems found in the statements, as messages; empty means all checks passed. -/
def statementProblems (env : Environment) : Array String := Id.run do
  let mut out := #[]
  let typeOf (n : Name) : Option Expr := (env.find? n).map (·.type)
  let moduleOf (n : Name) : Option Name := (env.getModuleIdxFor? n).bind fun i ↦ env.header.moduleNames[i.toNat]?
  for (n, m) in [(`Erdos40.Erdos40For, `FormalConjectures.ErdosProblems.«40»), (`Erdos40.erdos_40, `FormalConjectures.ErdosProblems.«40»), (`Erdos40.erdos_40.variants.weaker, `FormalConjectures.ErdosProblems.«40»), (`Erdos28.erdos_28, `FormalConjectures.ErdosProblems.«28»), (`Erdos158.B2, `FormalConjectures.ErdosProblems.«158»), (`Erdos158.erdos_158, `FormalConjectures.ErdosProblems.«158»)] do
    unless moduleOf n == some m do
      out := out.push s!"{n} is not from {m}"
  -- exists_iff : (the right side of erdos_40.variants.weaker) ↔ strong
  match (typeOf `Erdos40.erdos_40.variants.weaker).bind Expr.iff?, (typeOf `ErdosProblem40.exists_iff).bind Expr.iff? with
  | some (_, fcRhs), some (lhs, rhs) =>
    unless lhs == fcRhs do out := out.push "exists_iff: left side is not the right side of erdos_40.variants.weaker"
    unless rhs == strong do out := out.push "exists_iff: right side is not Erdos40For (fun _ ↦ 1)"
  | _, _ => out := out.push "exists_iff or erdos_40.variants.weaker is missing or not an ↔"
  -- answerSet_nonempty_iff : (the set in erdos_40).Nonempty ↔ strong
  match (typeOf `Erdos40.erdos_40).bind Expr.eq?, (typeOf `ErdosProblem40.answerSet_nonempty_iff).bind Expr.iff? with
  | some (_, fcSet, _), some (lhs, rhs) =>
    unless lhs.isAppOfArity `Set.Nonempty 2 && lhs.appArg! == fcSet do out := out.push "answerSet_nonempty_iff: left side is not the set of erdos_40, nonempty"
    unless rhs == strong do out := out.push "answerSet_nonempty_iff: right side is not Erdos40For (fun _ ↦ 1)"
  | _, _ => out := out.push "answerSet_nonempty_iff or erdos_40 is missing or has another shape"
  -- strong_implies_erdos_28 : strong → (the type of erdos_28)
  match typeOf `ErdosProblem40.strong_implies_erdos_28, typeOf `Erdos28.erdos_28 with
  | some (.forallE _ hyp body _), some e28 =>
    unless hyp == strong do out := out.push "strong_implies_erdos_28: hypothesis is not Erdos40For (fun _ ↦ 1)"
    unless !body.hasLooseBVars && body == e28 do out := out.push "strong_implies_erdos_28: conclusion is not the type of erdos_28"
  | _, _ => out := out.push "strong_implies_erdos_28 or erdos_28 is missing or has another shape"
  -- strong_implies_erdos_158 : strong → (the right side of erdos_158)
  let fc158 := (typeOf `Erdos158.erdos_158).bind Expr.iff? |>.map (·.2)
  match typeOf `ErdosProblem40.strong_implies_erdos_158, fc158 with
  | some (.forallE _ hyp body _), some rhs158 =>
    unless hyp == strong do out := out.push "strong_implies_erdos_158: hypothesis is not Erdos40For (fun _ ↦ 1)"
    unless !body.hasLooseBVars && body == rhs158 do out := out.push "strong_implies_erdos_158: conclusion is not the right side of erdos_158"
  | _, _ => out := out.push "strong_implies_erdos_158 or erdos_158 is missing or has another shape"
  -- strong_iff_forall_B2 : strong ↔ ∀ g A, B2 g A → L A; at g := 2 the right side is erdos_158's `∀ A, A.Infinite → B2 2 A → L A` without `A.Infinite`
  let forallB2 := (typeOf `ErdosProblem40.strong_iff_forall_B2).bind Expr.iff?
  match forallB2, fc158 with
  | some (lhs, rhs@(.forallE _ (.const `Nat []) gBody _)), some (.forallE _ tA' (.forallE _ _ fcA _) _) =>
    unless lhs == strong do out := out.push "strong_iff_forall_B2: left side is not Erdos40For (fun _ ↦ 1)"
    match gBody.instantiate1 (mkNatLit 2) with
    | .forallE _ tA oursA _ =>
      unless tA == tA' && !fcA.hasLooseBVar 0 && oursA == fcA.lowerLooseBVars 1 1 do out := out.push "strong_iff_forall_B2: right side at g = 2 is not the right side of erdos_158 without A.Infinite"
    | _ => out := out.push "strong_iff_forall_B2: right side has no ∀ A"
    -- weaker_iff_forall_B2 : (the right side of erdos_40.variants.weaker) ↔ (the right side of strong_iff_forall_B2)
    match (typeOf `Erdos40.erdos_40.variants.weaker).bind Expr.iff?, (typeOf `ErdosProblem40.weaker_iff_forall_B2).bind Expr.iff? with
    | some (_, fcRhs), some (lhs6, rhs6) =>
      unless lhs6 == fcRhs && rhs6 == rhs do out := out.push "weaker_iff_forall_B2: a side is not the right side of erdos_40.variants.weaker or of strong_iff_forall_B2"
    | _, _ => out := out.push "weaker_iff_forall_B2 is missing or not an ↔"
  | _, _ => out := out.push "strong_iff_forall_B2 or erdos_158 is missing or has another shape"
  return out

unsafe def main : IO UInt32 := do
  let env ← importModules (targets.map fun m ↦ { module := m }) Options.empty 0
  let mut problems : Nat := 0
  let mut checked : Nat := 0
  for m in targets do
    let some idx := env.getModuleIdx? m
      | throw (IO.userError s!"module {m} is not in the loaded environment")
    for (name, _) in env.constants.toList do
      if env.getModuleIdxFor? name == some idx then
        let (axs, _) ← ((Lean.collectAxioms name : CoreM (Array Name))).toIO
          { fileName := "probe", fileMap := default } { env := env }
        checked := checked + 1
        let extra := axs.filter fun a ↦ !standardAxioms.contains a
        unless extra.isEmpty do
          problems := problems + 1
          IO.println s!"NONSTANDARD {name}: {extra.toList}"
        if required.contains name then
          IO.println s!"{name}: {axs.toList}"
  for r in required do
    unless env.contains r do
      problems := problems + 1
      IO.println s!"MISSING {r}"
  let stmt := statementProblems env
  for s in stmt do
    IO.println s!"STATEMENT {s}"
  problems := problems + stmt.size
  if stmt.isEmpty then
    IO.println "statements: exists_iff, answerSet_nonempty_iff, strong_implies_erdos_28, strong_implies_erdos_158, strong_iff_forall_B2, weaker_iff_forall_B2 match FormalConjectures"
  IO.println s!"probe: {checked} constants checked, {problems} problems"
  return if problems == 0 then 0 else 1
