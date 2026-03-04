# Aligner SAM Output Specification

Complete reference for all SAM fields, custom tags, and SV-related signals emitted by the aligner.

---

## Table of Contents

- [Standard SAM Columns](#standard-sam-columns)
- [FLAG Bits](#flag-bits)
- [CIGAR Operations](#cigar-operations)
- [Core Alignment Tags](#core-alignment-tags)
- [Seeding & Chaining Tags](#seeding--chaining-tags)
- [Paired-End Tags](#paired-end-tags)
- [Supplementary & Split-Read Tags](#supplementary--split-read-tags)
- [SV Detection Tags](#sv-detection-tags)
- [Repeat & Genomic Context Tags](#repeat--genomic-context-tags)
- [Haplotype Phasing Tags](#haplotype-phasing-tags)
- [Filtering & Pipeline Tags](#filtering--pipeline-tags)
- [Unmapped Records](#unmapped-records)
- [CLI Flags That Control Tag Emission](#cli-flags-that-control-tag-emission)
- [Tag Summary by Use Case](#tag-summary-by-use-case)

---

## Standard SAM Columns

| Col | Field   | Description |
|-----|---------|-------------|
| 1   | QNAME   | Read name from FASTQ |
| 2   | FLAG    | Bitwise flags (see below) |
| 3   | RNAME   | Reference sequence name (chromosome), `*` if unmapped |
| 4   | POS     | 1-based leftmost mapping position, `0` if unmapped |
| 5   | MAPQ    | Mapping quality (0-60 scale) |
| 6   | CIGAR   | Alignment operations (see below), `*` if unmapped |
| 7   | RNEXT   | Mate reference name (`=` if same chromosome, `*` if unpaired/unmapped) |
| 8   | PNEXT   | 1-based mate position, `0` if unpaired/unmapped |
| 9   | TLEN    | Signed template length (positive for leftmost mate, negative for rightmost) |
| 10  | SEQ     | Read sequence (reverse-complemented if on minus strand) |
| 11  | QUAL    | Base qualities (Phred+33 encoded) |

---

## FLAG Bits

| Bit    | Hex    | Name                | Meaning |
|--------|--------|---------------------|---------|
| 1      | 0x1    | PAIRED              | Read is part of a pair |
| 2      | 0x2    | PROPER_PAIR         | Both mates mapped in expected orientation and insert size |
| 4      | 0x4    | UNMAPPED            | This read is unmapped |
| 8      | 0x8    | MATE_UNMAPPED       | Mate is unmapped |
| 16     | 0x10   | REVERSE             | Read aligned to reverse complement strand |
| 32     | 0x20   | MATE_REVERSE        | Mate aligned to reverse complement strand |
| 64     | 0x40   | FIRST_IN_PAIR       | This is the R1 mate |
| 128    | 0x80   | SECOND_IN_PAIR      | This is the R2 mate |
| 256    | 0x100  | SECONDARY           | Secondary alignment (not primary) |
| 2048   | 0x800  | SUPPLEMENTARY       | Supplementary alignment (chimeric/split read) |

---

## CIGAR Operations

| Op | Meaning | Notes |
|----|---------|-------|
| M  | Alignment match (match or mismatch) | Default mode |
| =  | Sequence match | With `--cigar-eqx` (not yet exposed as CLI flag) |
| X  | Sequence mismatch | With `--cigar-eqx` |
| I  | Insertion to reference | Bases present in read but not reference |
| D  | Deletion from reference | Bases present in reference but not read |
| S  | Soft clip | Unaligned bases at read ends (kept in SEQ) |

Indels are left-aligned per SAM spec convention.

---

## Core Alignment Tags

Always present on mapped reads.

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `NM` | i    | `NM:i:3` | Edit distance: total mismatches + inserted bases + deleted bases |
| `MD` | Z    | `MD:Z:50A49` | Mismatch/deletion string per SAM spec (for SNP/indel calling without reference) |
| `AS` | i    | `AS:i:145` | Alignment score from the extension stage |
| `ms` | i    | `ms:i:145` | minimap2-compatible match score (same value as AS) |
| `de` | f    | `de:f:0.0200` | Per-base divergence = 1.0 - (matches / aligned_length) |
| `rl` | i    | `rl:i:150` | Original read length in bp |
| `nn` | i    | `nn:i:0` | Count of N bases in the aligned reference span |

---

## Seeding & Chaining Tags

Present on all mapped reads. Useful for understanding why a read mapped where it did.

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `Xs` | i    | `Xs:i:42` | Total seeds collected for this read across all seed groups |
| `Xc` | i    | `Xc:i:35` | Primary chain score (DP chaining score before extension) |
| `cm` | i    | `cm:i:18` | Number of matching seeds in the primary chain (minimap2-compatible) |
| `s1` | i    | `s1:i:35` | Best chain score (minimap2-compatible) |
| `s2` | i    | `s2:i:12` | Second-best chain score, 0 if only one chain (minimap2-compatible) |
| `tp` | A    | `tp:A:P` | Alignment type: `P` = primary, `S` = supplementary, `I` = supplementary inverted |
| `Xp` | i    | `Xp:i:0` | Pipeline path: `0` = standard alignment, `5` = high-confidence fast path (skipped secondary extensions) |

---

## Paired-End Tags

Present only in paired-end mode.

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `MC` | Z    | `MC:Z:150M` | Mate's CIGAR string |
| `MQ` | i    | `MQ:i:60` | Mate's mapping quality |
| `YD` | f    | `YD:f:1.50` | Insert size Z-score: \|observed - mean\| / stddev. Values > 3 suggest discordance. |
| `YJ` | f    | `YJ:f:0.9523` | Joint posterior probability of this pair configuration (among top-k pair candidates) |
| `YT` | Z    | `YT:Z:CP` | Pair type classification (see table below) |
| `Xr` | i    | `Xr:i:0` | Mate rescue status: `0` = both mapped independently, `1` = normal pair, `2` = mate was rescued via local realignment |

### YT Pair Type Values

| Value | Meaning | SV Relevance |
|-------|---------|--------------|
| `CP`  | Concordant proper pair | Normal — expected orientation and insert size |
| `DO`  | Discordant orientation | Suggests inversion or translocation |
| `DI`  | Discordant insert size | Suggests deletion (large insert) or duplication (small insert) |
| `IC`  | Inter-chromosomal | Suggests translocation |
| `UU`  | Both unmapped | No mapping |

---

## Supplementary & Split-Read Tags

Present on reads with chimeric/split alignments (primary + supplementary records).

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `SA` | Z    | `SA:Z:chr2,1200,+,95M55S,40,2;` | Supplementary alignments. Semicolon-delimited, each entry: `rname,pos,strand,CIGAR,MAPQ,NM`. Trailing semicolon per SAM spec. Present on both primary and supplementary records. |
| `XA` | Z    | `XA:Z:chr1,+500,100M,2;chr3,-800,90M5I5M,4` | Alternative alignment locations from top-k evaluation. Format: `chr,strand+pos,CIGAR,NM` per entry. SE only. |
| `XS` | i    | `XS:i:120` | Score of the best secondary/supplementary chain |

### SA Tag Format Detail

```
SA:Z:chr,pos,strand,CIGAR,MAPQ,NM;chr,pos,strand,CIGAR,MAPQ,NM;...;
```

- **chr**: Reference name
- **pos**: 1-based position
- **strand**: `+` or `-`
- **CIGAR**: Full CIGAR (includes soft clips)
- **MAPQ**: Mapping quality of this supplementary segment
- **NM**: Edit distance of this segment
- Terminated with trailing `;`

---

## SV Detection Tags

These tags encode structural variant evidence directly in the SAM record. Present on reads with supplementary alignments (split reads).

### Breakpoint Tags

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `YS` | Z    | `YS:Z:DEL` | SV type classification. On the primary record: comma-separated list of all SV types from supplementaries. On supplementary records: that segment's SV type. |
| `YB` | Z    | `YB:Z:chr1:97-102\|chr1:196-203` | Breakpoint intervals (0-based, inclusive). Left and right breakpoints separated by `\|`. Interval width reflects positional uncertainty from microhomology. |
| `YC` | Z    | `YC:Z:5,7` | Confidence interval widths in bp for left and right breakpoints. CI = microhomology length + local repeat bonus. |
| `YM` | i    | `YM:i:5` | Microhomology length in bp at the breakpoint junction. Bidirectional scan (forward + backward). `0` if no microhomology. |
| `YI` | Z    | `YI:Z:ACGTACGT` | Inserted sequence at the breakpoint junction. Only present when uncovered read bases exist between primary and supplementary segments. Handles RC orientation. |
| `YA` | Z    | `YA:Z:100,200;101,201` | Alternative breakpoint placements within the microhomology zone. Each entry is `left_pos,right_pos`. Only present when microhomology > 0. |

### YS SV Type Values

| Value | Meaning | Evidence Pattern |
|-------|---------|-----------------|
| `DEL` | Deletion | Split read with gap on reference, same strand, same chromosome |
| `INV` | Inversion | Split read segments on opposite strands, same chromosome |
| `TRA` | Translocation | Split read segments on different chromosomes |
| `DUP` | Duplication | Split read with overlapping reference coordinates, same strand |

### SV Confidence & Filtering Tags

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `XN` | i    | `XN:i:1` | Germline SV flag. `1` = breakpoint matches a known germline SV in the panel-of-normals BED. `0` = novel. Only present when `--panel-of-normals` is provided. |
| `XE` | f    | `XE:f:0.85` | Split chain entropy. Low entropy = confident single split point. High entropy = ambiguous or complex rearrangement. SE only. |
| `XD` | i    | `XD:i:50` | Score delta between primary and best supplementary chain. Large delta = confident primary, small delta = ambiguous. |
| `XR` | f    | `XR:f:2.500` | Ratio of primary to second-best chain score. Higher = more confident primary mapping. |
| `XF` | f    | `XF:f:0.05` | Maximum artifact probability among supplementary alignments. High values suggest the split may be artifactual. |
| `XU` | Z    | `XU:Z:REF:0.10,SV:0.85,RPT:0.05` | Uncertainty classification posteriors. Three classes: REF (reference-consistent), SV (structural variant), RPT (repeat/mapping artifact). Sum to ~1.0. SE only. |

---

## Repeat & Genomic Context Tags

Provide genomic context for each alignment. Useful for filtering and interpretation.

### Copy Number & Mappability

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `XG` | f    | `XG:f:0.45` | GC fraction of the aligned reference region. Useful for GC-bias correction in CNV calling. |
| `XW` | f    | `XW:f:0.85` | Local mappability estimate: `min(1.0 / avg_seed_occurrence, 1.0)`. Low values indicate repetitive regions. |

### Repeat Classification

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `XP` | Z    | `XP:Z:TAN` | Repeat class of the mapping region. Only present when 3+ near-equal chains exist. SE only. |
| `Xe` | Z    | `Xe:Z:LINE` | Mobile element class if the read maps to a known ME. Only present when detected. SE only. |

#### XP Repeat Class Values

| Value | Meaning |
|-------|---------|
| `UNQ` | Unique (non-repetitive) |
| `TAN` | Tandem repeat |
| `SEG` | Segmental duplication |
| `INV` | Inverted repeat |
| `DIS` | Dispersed repeat |

#### Xe Mobile Element Values

`LINE`, `SINE`, `DNA`, `LTR`, `SVA`, `Alu`

### Tandem Repeat Annotation

Present when a tandem repeat BED annotation is provided via `--tr-bed`.

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `TR` | Z    | `TR:Z:CAG` | Repeat unit sequence |
| `TN` | f    | `TN:f:10.5` | Number of repeat unit copies (fractional) |
| `TL` | i    | `TL:i:3` | Repeat unit length in bp |

---

## Haplotype Phasing Tags

Present when `--phase-vcf` is provided with a phased VCF file.

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `HP` | i    | `HP:i:1` | Haplotype assignment: `1` or `2`. Only present when the read overlaps heterozygous sites with sufficient evidence. |
| `PS` | i    | `PS:i:12345` | Phase set ID. Reads sharing a PS value are phased relative to each other. Derived from VCF phase blocks or positional clustering. |

---

## Filtering & Pipeline Tags

### Top-K Evaluation (SE only)

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `YK` | i    | `YK:i:5` | Number of chains extended during top-k evaluation (controlled by `--top-k`). More chains = better MAPQ calibration. |
| `YP` | f    | `YP:f:0.9500` | Boltzmann posterior probability of the primary alignment among all extended chains. Values near 1.0 = high confidence; low values = ambiguous mapping. |

### Soft-Clip Realignment

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `XK` | Z    | `XK:Z:...` | Soft-clip realignment results. Encodes where large soft-clipped segments can be realigned. Present when clips are successfully remapped. |

### Two-Pass Realignment

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `OC` | Z    | `OC:Z:100M50S` | Original CIGAR before two-pass SV-hotspot realignment. Only present on records that were realigned in Pass 2 (with `--two-pass`). |

### Read Group

| Tag  | Type | Example | Description |
|------|------|---------|-------------|
| `RG` | Z    | `RG:Z:sample1` | Read group ID. Present when `--read-group` / `-R` is specified. The ID is extracted from the `@RG` header line. |

---

## Unmapped Records

Unmapped reads (FLAG & 0x4) have minimal fields:

- **RNAME**: `*`
- **POS**: `0`
- **MAPQ**: `0`
- **CIGAR**: `*`
- **SEQ/QUAL**: Present (original read sequence and qualities)
- **Tags**: Only `de` and `rl` are emitted; all alignment-derived tags are absent.

In paired-end mode, if only one mate is unmapped, it inherits the mapped mate's RNAME and POS (per SAM spec), and FLAG includes both 0x4 (unmapped) and 0x8 (mate unmapped) as appropriate.

---

## CLI Flags That Control Tag Emission

| Flag | Tags Affected | Effect |
|------|---------------|--------|
| `--top-k N` | `YK`, `YP`, `XA` | Controls how many chains are extended. Higher values produce `XA` alternatives and better `YP` posteriors. |
| `--two-pass` | `OC`, `YS`, `YB`, `YC`, `YM`, `YI`, `YA` | Enables two-pass SV-hotspot realignment. First pass collects SV signals; second pass realigns reads at hotspots. Produces `OC` on realigned records. |
| `--panel-of-normals <bed>` | `XN` | Enables germline SV filtering. Reads at known germline breakpoints get `XN:i:1`. |
| `--phase-vcf <vcf>` | `HP`, `PS` | Enables haplotype phasing from a phased VCF. |
| `--tr-bed <bed>` | `TR`, `TN`, `TL` | Enables tandem repeat annotation. |
| `--sv-signal <bed>` | (BED file output) | Writes SV signal clusters to a BED file (not a SAM tag). |
| `--scaffold-bed <bed>` | (affects two-pass hotspots) | Loads long-read SV breakpoint hints and merges them into two-pass hotspot regions. |
| `--read-group` / `-R` | `RG` | Adds read group tag + `@RG` header line. |
| `--debug-align [file]` | (JSON debug output) | Writes per-read JSON trace with seeds, chains, extension details, MAPQ breakdown. Not in SAM. |
| `--stats <json>` | (JSON stats output) | Writes summary statistics (MAPQ histogram, mapping rate, identity distribution). Not in SAM. |
| `--clip-overlap` | `OC` | Clips overlapping PE bases at 3' end of forward mate; stores original CIGAR in `OC`. |

---

## Tag Summary by Use Case

### For SV Callers

The primary tags for building an SV caller on top of this aligner:

| Signal | Tags | How to Use |
|--------|------|------------|
| **Split reads** | `SA`, `YS`, `YB`, `YC` | `SA` identifies chimeric alignments. `YS` classifies the SV type. `YB` gives breakpoint intervals with `YC` confidence widths. |
| **Breakpoint precision** | `YM`, `YI`, `YA` | `YM` reports microhomology (ambiguity). `YI` captures inserted sequence at junctions. `YA` lists alternative placements within the microhomology zone. |
| **Discordant pairs** | `YT`, `YD` | `YT` classifies pair discordance type (`DO`=inversion, `DI`=size, `IC`=translocation). `YD` Z-score quantifies insert size deviation. |
| **Confidence** | `XU`, `XE`, `XD`, `XR`, `XF` | `XU` posteriors classify REF/SV/RPT. `XE` entropy, `XD` score delta, `XR` ratio, and `XF` artifact probability provide layered confidence. |
| **Germline filter** | `XN` | Skip `XN:i:1` records to filter known germline SVs. |
| **Repeat context** | `XP`, `Xe`, `TR`/`TN`/`TL` | Annotate whether SVs overlap repeats, mobile elements, or tandem repeat expansions. |

### For CNV Callers

| Signal | Tags | How to Use |
|--------|------|------------|
| **Read depth** | `XG`, `XW` | `XG` for GC-bias correction. `XW` for mappability masking. |
| **Insert size** | `YD`, `YT` | `YD` Z-score for insert size anomalies indicating CNVs. |

### For Variant Callers (SNV/Indel)

| Signal | Tags | How to Use |
|--------|------|------------|
| **Alignment quality** | `NM`, `MD`, `AS`, `de` | Standard tags for variant quality scoring. `de` divergence for filtering. |
| **Mapping confidence** | `MAPQ`, `YP`, `s2` | MAPQ for filtering. `YP` posterior and `s2` second-best score for fine-grained confidence. |
| **Haplotype** | `HP`, `PS` | Phase variants using read-level haplotype assignments. |
| **Tandem repeats** | `TR`, `TN`, `TL` | Flag variants in tandem repeat regions for special handling. |

### For Quality Control

| Signal | Tags | How to Use |
|--------|------|------------|
| **Mapping stats** | `de`, `rl`, `nn`, `MAPQ` | Per-read divergence, length, N-content, mapping quality. |
| **Pipeline diagnostics** | `Xs`, `Xc`, `Xp`, `cm`, `Xr` | Seed count, chain score, pipeline path, rescue status. |
| **Pair metrics** | `YD`, `YT`, `YJ` | Insert size Z-score, pair type, pair posterior. |

---

## Header Lines

The aligner emits standard SAM header lines:

```
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:248956422
...
@RG	ID:sample1	SM:sample1	...    (if --read-group provided)
@PG	ID:aligner	PN:aligner	VN:x.y.z	CL:...
```
