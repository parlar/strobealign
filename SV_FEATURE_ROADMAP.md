# SV Feature Roadmap: From Spec to Implementation

## Context

FUTURE_SAM_OUTPUT.md defines a comprehensive set of SV-aware SAM tags. This document
describes which features to implement first, why they matter, and how they fit into a
concrete use case: CYP2D6 pharmacogenomics star allele calling (cyrius-rs project).

CYP2D6 is one of the hardest loci in the human genome — a 4.4kb gene with a 97%
identical pseudogene (CYP2D7) 13kb downstream, connected by an 8kb spacer region.
Structural variants include tandem duplications (*2x2), gene deletions (*5), and
hybrid gene conversions (*36, *13, *68). Current callers fail on complex arrangements
because the linear reference can't represent tandem junctions, and reads spanning these
junctions get soft-clipped or mismapped.

This is the ideal test case for SV-aware alignment: a small region (60kb), well-studied
structural variants, ground truth available, and existing callers that fail on specific
samples due to missing structural signals.

## Priority 1: Supplementary Alignments (SA tag)

**What:** When a read aligns to two distinct locations (chimeric/split read), emit a
supplementary alignment record (FLAG 0x800) and the `SA:Z` tag on both primary and
supplementary records.

**Why this is the highest priority:**
- SA tags are the standard way to represent split reads. Every downstream SV caller
  (DELLY, Manta, GRIDSS, our sv_caller) consumes SA tags as primary evidence.
- Without SA, split reads appear as soft-clipped reads with no indication of where the
  clipped portion maps. The sv_caller has a clip_realign module to compensate, but it's
  slower and less accurate than having the aligner do it.
- For CYP2D6: reads spanning a *36 hybrid junction (D7 exon9 → D6 intron4) would
  produce SA tags showing the two mapping locations, directly proving the hybrid event.

**How:**
- After primary alignment, check if large soft clips (>30bp) can be aligned elsewhere.
  Use the existing NAM/chain machinery on the clipped portion.
- If the clipped portion aligns with sufficient score, emit a supplementary record with
  FLAG 0x800 and add SA:Z to both records.
- Format: `SA:Z:rname,pos,strand,CIGAR,MAPQ,NM;`
- The `tp:A` tag should distinguish: `P` = primary, `S` = supplementary.

**Implementation complexity:** Medium. The chaining/alignment infrastructure exists.
The main work is detecting good split points and managing the two-record output.

## Priority 2: Pair Type Classification (YT tag)

**What:** Classify each read pair into one of: `CP` (concordant proper), `DO`
(discordant orientation), `DI` (discordant insert size), `IC` (inter-chromosomal),
`UU` (both unmapped).

**Why:**
- Discordant pairs are the second most important SV signal after split reads.
- RF (reverse-forward) orientation = tandem duplication. This is the classic signature
  for CYP2D6 *2x2, *36+*10, etc.
- DI (large insert) = deletion. Would directly detect *5 gene deletions.
- The information is already available (mate position, orientation, insert size) — it
  just needs to be classified and tagged.

**How:**
- After paired-end alignment, compare observed orientation and insert size to expected.
- Use the already-computed insert size distribution (mu, sigma) from insertsize.rs.
- Classify: same-chrom FR within 3 sigma = CP, same-chrom RF = DO (dup), same-chrom
  FR but |insert| > mu + 3*sigma = DI, different chrom = IC.
- Emit `YT:Z:XX` tag.
- Also emit `YD:f` (insert size Z-score) for quantitative discordance.

**Implementation complexity:** Low. All inputs are already computed. Just needs
classification logic and tag emission.

## Priority 3: SV Type Tags (YS, YB, YM)

**What:** For reads with supplementary alignments, classify the SV type and compute
breakpoint coordinates with microhomology.

**Why:**
- `YS:Z:DUP` vs `YS:Z:DEL` tells the SV caller what structural event the read
  supports, without the caller having to re-derive it from the alignment geometry.
- `YB:Z` breakpoint intervals give precise positions with uncertainty from
  microhomology.
- `YM:i` microhomology length indicates the mechanism (NAHR for CYP2D6 events).

**How:**
- From primary + supplementary alignment coordinates:
  - Same strand, gap on reference, no gap on query → DEL
  - Same strand, overlap on reference → DUP
  - Opposite strands → INV
  - Different chromosomes → TRA
- Breakpoint = junction between primary and supplementary aligned segments.
- Microhomology = bidirectional scan from breakpoint for matching bases.

**Implementation complexity:** Medium. Requires SA tag (Priority 1) to be implemented
first.

## Priority 4: Soft-Clip Realignment (XK tag)

**What:** When a read has large soft clips that can't form a clean supplementary
alignment, attempt local realignment of the clipped bases and report results in `XK:Z`.

**Why:**
- Not all split reads produce clean supplementary alignments. In repetitive regions
  like CYP2D6's REP6, soft clips may hit low-complexity sequence where the chaining
  fails.
- XK provides a fallback: "these 58 clipped bases best align to chr22:42143719" even
  when the alignment isn't strong enough for a full supplementary record.
- For CYP2D6: the breakpoint reads at position 42130053 in sample NA18565 have
  38-107bp soft clips with repetitive content. SA tags might not fire, but XK could
  still capture the realignment signal.

**Implementation complexity:** Medium-high. Needs a local realignment strategy for
short, possibly repetitive sequences.

## Priority 5: Two-Pass SV Hotspot Realignment

**What:** First pass collects SV signals (soft-clip clusters, discordant pairs). Second
pass re-examines reads at these hotspots with relaxed alignment parameters, potentially
producing better alignments. Original CIGAR stored in `OC:Z`.

**Why:**
- Standard alignment parameters are tuned for speed on normal reads. At SV breakpoints,
  reads may need more aggressive extension or different gap penalties to align correctly.
- For CYP2D6: the D6/D7 homology region produces many low-MAPQ reads. A second pass
  focused on known breakpoint positions could rescue some of these.

**Implementation complexity:** High. Requires a first-pass signal collection framework
and a targeted re-alignment pass. This is the most complex feature.

## Priority 6: Repeat & Mappability Tags (XP, XW)

**What:** Classify the mapping region as UNQ/TAN/SEG/INV/DIS and report local
mappability.

**Why:**
- `XP:Z:SEG` would directly flag CYP2D6/D7 as a segmental duplication, telling
  downstream tools to expect ambiguous mappings.
- `XW:f` mappability estimate helps CNV callers adjust depth signals for region-specific
  biases.

**Implementation complexity:** Low-medium. XW can be derived from seed occurrence
counts (already computed). XP requires external annotation or runtime detection.

## How the System Should Work End-to-End

The goal is a pipeline where strobealign provides rich per-read SV signals that
downstream callers consume without redundant re-analysis:

```
Input: FASTQ + Reference
           |
    Strobealign (with SV features)
           |
    BAM with: SA, YT, YS, YB, YM, XK, XP, XW tags
           |
     +-----+-----+
     |             |
  General SV    cyrius-rs
  tools           |
     |          CYP2D6 star
  SV VCF       allele calls
```

For CYP2D6 calling (~/dev/cyrius-rs/):
- SA tags identify reads crossing hybrid junctions (D6↔D7)
- YT tags detect tandem dup signatures (RF pairs in the spacer region)
- Per-read SV evidence feeds into a structural diplotype scoring framework
- Replaces the current heuristic if/else pattern matching with probabilistic scoring

**Interim: algorithms from sv_tools brought into cyrius-rs and strobealign:**

Before strobealign's SV features are fully implemented, useful algorithms from
~/dev/sv_tools/ (experimental) can be copied directly into cyrius-rs or strobealign.
The principle is to bring the code in, not depend on sv_tools as a crate:

- **LocalIndex (sv_align)** — q-gram + Myers bit-parallel alignment, ~400 lines, zero
  deps. Useful for both strobealign (powering SA/XK tag generation) and cyrius-rs
  (realigning clips against known CYP2D6 structural junctions).
- **ClipCluster (sv_break)** — clip clustering logic for identifying breakpoints from
  soft-clip consensus. Algorithm adaptable into either project.
- **Per-exon depth ratios (sv_exon)** — interval depth computation for detecting
  partial duplications. Useful in cyrius-rs for CN gradient detection.
- **CBS segmentation (sv_cnv)** — changepoint detection for depth signals.

Once strobealign natively emits SA/YT/YS tags, the interim clip realignment in
cyrius-rs becomes a fallback for reads where the aligner can't produce clean
supplementary alignments (e.g., repetitive breakpoint sequences in REP6).

## Test Case: CYP2D6 Region (chr22:42100000-42160000)

Three specific samples from 1000 Genomes demonstrate why these features matter:

**NA18545** (expected *5/*36x2+*10x2): D6 CN=4, spacer CN=4, 180 problem reads.
Currently miscalled as *36+*10/*36+*10. SA tags on reads at the tandem junction would
confirm the multi-copy arrangement. YT:DI on pairs spanning the *5 deletion would
confirm the deletion.

**NA18565** (expected *1/*36+*10): D6 CN=3.08, spacer CN=2.7, 222 problem reads,
1 breakpoint contig at 42130053. Currently miscalled as *1/*2 because CN rounds to 2.
SA tags on the breakpoint reads (currently 38-107bp soft clips) would provide direct
structural evidence for the tandem arrangement, independent of depth.

**NA19143** (expected *1/*2x2): D6 CN=2.45, 122 problem reads. Exon-level depth shows
partial duplication (CN gradient from 1.7 to 3.0 across the gene). SA tags from reads
spanning the partial dup breakpoint would pinpoint where the duplication starts.

Full-region BAMs for these samples are at /tmp/fullwgs/*_full_region.bam. Detailed
analysis is in ~/dev/cyrius-rs/ANALYSIS_3_WRONG_SAMPLES.md.

## Implementation Order Summary

| Priority | Feature | Tags | Complexity | Impact |
|----------|---------|------|------------|--------|
| 1 | Supplementary alignments | SA, tp | Medium | Critical — enables all downstream SV calling |
| 2 | Pair type classification | YT, YD | Low | High — tandem dup detection via RF pairs |
| 3 | SV type + breakpoints | YS, YB, YM | Medium | High — pre-computed SV evidence per read |
| 4 | Clip realignment | XK | Medium-high | Medium — fallback for ambiguous splits |
| 5 | Two-pass realignment | OC | High | Medium — rescues reads at known hotspots |
| 6 | Repeat/mappability | XP, XW | Low-medium | Low — annotation, not direct evidence |

Priorities 1-2 would unlock the most value for the least implementation effort.
Priority 1 (SA tags) is the single most impactful feature — it's the foundation that
the sv_caller, cyrius-rs, and every other SV tool depends on.
