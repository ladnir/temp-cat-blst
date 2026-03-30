# Lazy Renorm Strategy

This note captures the current methodology for placing lazy repair and
representation-conversion operations in AVX field circuits.

The original question was where to place `LUFp` `renorm` and `normalize`, but
the more general and more accurate framing is:

- this is a representation-placement problem over a circuit DAG.

## Goal

Use lazy representations to make local combine regions cheaper than strict
representations, without paying more in repair and conversion overhead than we
save.

The key principle is:

- defer repair and conversion work as long as the circuit remains within a
  documented bound and no helper precondition is violated.

## General Framing: Representation Placement

Each intermediate in the circuit should be thought of as carrying a
representation state, not just a value.

At minimum, the representation state should include:

- width:
  - narrow
  - wide
- low-half slack class
- high-half drift class
- legality for the next helper:
  - add/sub legal
  - multiply legal
  - normalize legal

Examples of representation dimensions:

- `(narrow, exact)`
- `(wide, LNeg0_Pos1_HNeg0_Pos1)`
- `(wide, LNeg1_Pos2_HNeg2_Pos4)`

The placement problem is then:

- where to stay in the current representation;
- where to widen;
- where to narrow;
- where to `renorm`;
- where to `normalize`.

This is a DAG optimization problem. Fanout matters because a single earlier
repair or conversion may be cheaper than repeating that cost on every outgoing
edge.

## Representation Model

Model a widened value as:

- `X = lo + R * hi`
- `R = 2^384`

Track bounds separately:

- low-half slack in units of `R`
- high-half slack in units of `p`

Recommended naming:

- `LNegA_PosB`: `-A*R <= lo < B*R`
- `HNegA_PosB`: `-A*p <= hi < B*p`

Examples:

- narrow `Fp` embedded in widened form:
  - `LNeg0_Pos1_HNeg0_Pos0`
- strict widened `UFp`:
  - `LNeg0_Pos1_HNeg0_Pos1`

If future lazy narrow types are introduced, the same methodology should be
applied there too. The representation state should stay uniform even when the
width changes.

## Why Low Slack Matters

If low limbs are forced canonical after every conceptual add/sub, then even a
single low carry can be conservatively rounded into an extra full `+p` in the
high range. That is too pessimistic.

Allowing low-half slack means:

- carries/borrows stay local for longer;
- high bounds do not inflate prematurely;
- one final bounded normalize can often replace many intermediate repairs.

## Mixed Wide And Narrow Inputs

Treat narrow as widened with exact zero high half:

- narrow embedded: `LNeg0_Pos1_HNeg0_Pos0`

That means mixed operations are tracked the same way as wide-wide operations,
except one operand starts with exact `hi = 0`.

Important consequence:

- wide + narrow does not automatically require a full extra `+p` in high;
- wide - narrow does not automatically require a full extra `-p` in high.

The low carry/borrow can remain deferred until `renorm` or `normalize`.

## Default Decision Rule

If a circuit:

- has a single output,
- is a linear add/sub chain,
- has no helper precondition issues,
- and will not overflow the allowed bound,

then the default strategy is:

- do not `renorm` in the middle;
- do one bounded final `normalize`.

This is the baseline lazy strategy.

## When Early Renorm Is Worth It

Early `renorm` becomes attractive when one intermediate feeds multiple
downstream consumers.

Example:

```text
z = w - t1 - t2
u = z + a
v = z - b
```

Here `z` fans out to two branches.

If `z` is renormed once before the split:

- both downstream branches start from a tighter bound;
- one renorm may replace two more expensive final repairs.

So the first question is not "how many ops deep is this chain?" but rather:

- where does the circuit fan out?

This same logic applies to width conversion:

- if one wide intermediate fans out to many consumers, narrowing it once may be
  cheaper than narrowing each branch separately;
- if one narrow intermediate fans out to many wide combine regions, widening it
  once may likewise be cheaper than repeating the promotion.

## Internal Split Helpers

Even a circuit with one final output can contain internal fanout.

Example:

```text
(x0, x1) -> (x0 - x1, x0 + x1)
```

This is effectively a split point. If the input is very wide, both outputs
inherit that width. A selective `renorm` before such a helper can be cheaper
than letting both branches carry inflated bounds.

This is why helpers like `mul_by_u_plus_1` must be considered carefully.

The same reasoning will apply to any future lazy narrow helper that turns one
value into two independently used outputs.

## Multiply As A Representation Boundary

If we do not introduce a wider-than-wide representation, multiply is a hard
boundary.

That means:

- lazy wide add/sub regions can usually span only the combine space between
  multiplies;
- before feeding a value into multiply again, it must be brought into a
  multiply-legal representation;
- future lazy narrow types may become valuable specifically because they can
  reduce the cost of that boundary.

This is why the full problem is at least two-dimensional:

- width
- drift

and in practice more than two-dimensional once legality and fanout are included.

## Repair Strategy

`renorm` and `normalize` should be treated as members of the same repair family.
The difference is the target bound class:

- `renorm`
  - repairs into a cheaper lazy target class;
- `normalize`
  - repairs into a stricter final consumer-legal class.

Neither should be thought of as "single-`p` only" by default. If the source
bound is wide enough, both may need scaled-`p` repair.

Final normalize should be fixed-pass SIMD, not repeated single-`p` correction.

If the source or final high range is bounded by a power-of-two multiple of `p`,
use a scaled-`p` repair ladder.

Example:

- for a class bounded by roughly `HNeg8_Pos8`:
  - conditional `+8p`
  - conditional `-8p`
  - conditional `-4p`
  - conditional `-2p`
  - conditional `-p`

This scales with `log2(bound)` and is compatible with vector masks.

The same idea can apply to `renorm`.

Example:

- if `renorm` needs to move from a class like `HNeg8_Pos8` to something tighter
  such as `HNeg0_Pos2`, then a scaled-`p` ladder may also be the right fixed-pass
  vector strategy there.

## Repair Implementation Options

At this stage, the methodology should make the repair options explicit without
mandating a single implementation style yet. Different bound classes may prefer
different kernels.

### Option 1: One-Sided Ladder

Use this when the source class is already one-sided or effectively one-sided.

Examples:

- `HNeg0_Pos8 -> HNeg0_Pos2`
- `HNeg4_Pos1 -> HNeg0_Pos1`

Shape:

- only upper repair, or only lower repair;
- scaled-`p` ladder on one side only.

This is the simplest case and should be preferred whenever the bound allows it.

### Option 2: Biased One-Sided Ladder

For a two-sided source class, unconditionally bias by a fixed multiple of `p`
to shift the bound entirely to one side, then run a one-sided ladder.

Example:

- `HNeg2_Pos4`
- unconditional `+2p`
- transformed class: `HNeg0_Pos6`
- upper-side ladder only

Advantages:

- simpler than a true two-sided repair;
- no sign-side split in the main ladder;
- often good when both negative and positive drift are material.

Disadvantages:

- the bias cost is paid on every lane;
- the one-sided bound becomes wider;
- may lose if the source class is already strongly skewed to one side.

### Option 3: Ternary Symmetric Ladder

At each scale `k`, classify lanes into:

- negative
- in-range
- positive

and apply one correction from `{+k*p, 0, -k*p}`.

This is a three-way per-scale repair step.

Advantages:

- often lower latency than maintaining separate negative and positive ladders;
- especially attractive for balanced two-sided classes such as:
  - `HNeg4_Pos4`
  - `HNeg8_Pos8`
- does not pay unconditional bias on lanes already near the target.

Disadvantages:

- needs a three-way lane classification;
- may be a little more complex in masks/registers than a one-sided ladder.

### Option 4: Independent Two-Sided Ladders

Conceptually:

- upper-side repair ladder
- lower-side repair ladder

This is the simplest way to think about two-sided repair algebraically, but it
is not necessarily the preferred optimized implementation.

It remains useful as:

- a correctness baseline;
- a reference model for source/target class reasoning;
- a fallback implementation if a ternary kernel is awkward for a given class.

But for latency-oriented optimized kernels, ternary repair will often be a
better default starting point than fully independent upper and lower ladders.

## Choosing A Repair Kernel

The methodology does not yet mandate one repair style globally. Instead:

- one-sided bounds should use one-sided ladders;
- balanced two-sided bounds should consider ternary repair first;
- two-sided bounds with large offset may consider biased one-sided repair;
- independent two-sided ladders remain a reference option, not the default
  optimized target.

The final decision should be made by benchmarked implementation, not by theory
alone.

## Compile-Time Selection Plan

Later, after benchmarking, repair kernels should be selected by compile-time
logic based on the documented source and target bound classes.

That implies:

- each source/target transform may have more than one implementation candidate;
- the code should keep the bound contract fixed even if the chosen kernel changes;
- benchmark results should drive which implementation is selected for each class.

The bound contract should be stable. The repair implementation behind that
contract may evolve.

## Practical Methodology

For any candidate circuit:

1. Identify exact representation state for each input:
   - width
   - low slack
   - high drift
   - helper legality
2. Propagate that state through each operation.
3. Mark every point where one value fans out into multiple consumers.
4. Mark every point where width changes are possible:
   - narrow -> wide
   - wide -> narrow
5. Mark every helper whose precondition requires a tighter bound or a different
   legality class.
6. Delay `renorm` unless fanout, overflow risk, or a helper precondition forces it.
7. Delay width conversion unless fanout or helper legality makes an earlier
   shared conversion cheaper.
8. Prefer one final bounded `normalize` over many intermediate `renorm`s when
   the final bounded class still admits an efficient fixed-pass repair.
9. When repair is needed, choose the cheapest repair ladder for the source and
   target classes; this may be scaled-`p` for either `renorm` or `normalize`.

## Heuristic Summary

- Linear single-output chain: defer repair.
- Fanout point: consider renorm before the split.
- Split helper: treat as fanout.
- Width conversion is also a placement decision.
- Shared conversion before fanout can beat repeated per-branch conversion.
- Overflow risk: renorm before the bound is violated.
- Final normalize should absorb as much remaining bounded drift as practical.

## What This Means For Fp6

Current working hypothesis:

- `Fp6 mul`
  - likely best with zero intermediate `renorm`;
  - preserve the baseline widened schedule;
  - rely on one bounded final normalize if the final class allows it.

- `Fp6 sqr`
  - likely needs at most one selective `renorm`;
  - the most likely location is the long `c2` subtraction chain or just before a
    split-style helper if bounds widen too much.

If lazy narrow representations are added later, they should be evaluated with
the same topology-aware methodology. The question should not be "is lazy narrow
good?" but rather:

- where in the DAG does narrow-vs-wide conversion minimize total repair cost?

## Circuit Cost View

One useful way to think about the problem is:

- each node has a representation state;
- each edge/helper has a cost and a legality contract;
- each repair or conversion has a cost;
- fanout multiplies downstream conversion cost.

So the real optimization target is:

- minimize total cost over the circuit DAG subject to legality and bound
  constraints.

`renorm` placement is one part of that optimization. Width placement is another.

## Rule For Future Work

When evaluating a placement decision, ask:

- is this value on a single linear path, or is it about to fan out?
- does the next helper force a tighter drift bound?
- does the next helper force a width change?
- can one shared conversion or `renorm` replace several downstream ones?

If it is only on a single linear path and bounds remain legal, the default is:

- do not repair or convert yet.
