# Figure 3 + S3 — CpG methylation, Fiber-seq accessibility, 45S/RNU2

Long-read methylation/accessibility. Streaming stages ran on the servers; final panels are local.

| Panel | Script | Where |
|---|---|---|
| 3A / S3A,B array methylation vs position | `40_methylation_full215*.py` | Mac (+ server stream) |
| 3B / S3C per-copy methylation distribution | `40_...`, `41_methylation_overview*.py` | Mac |
| 3C,D / S3F,G within-copy positional by class | `40_methylation_full215.py` (authoritative) | Mac |
| 3E dosage compensation (set-point) | `40_...`, `36_lowmeth_copies_donor_stability.py` | Mac |
| 3F,G,H whole-copy meth by variant count | `60_variant_hypo_proportion*.py`, `24b–g` | Mac |
| S3D,E array-edge / border methylation | `62_/64_/65_/66_*` | Mac |
| **3E,F / S3H,I Fiber-seq m6A accessibility** | `fiberseq_5S/scripts/build_spanning_figure.py`, `within_copy_profile_pub.py` | Mac (+ trr237 stream) |
| **S3K,L 45S NOR edge methylation** | `45S_methylation/scripts/{00–04}` | Mac (halves: immuno2 39 + trr237 109) |
| S3M RNU2 own-assembly methylation | `rDNA_dosage_control/09.../RNU2/compute_rnu2_ownasm.py` | Mac |
| S3N edge-biased gene-body variation | `nascent_edge_analysis/edge_analysis.py`, `edge_expression.py` | Mac |
| DB `copy_methylation` build | `37_/38_*`, `60_multicontig_methylation.py` | Mac + immuno2 |

Notes: use `03_methylation_full215` (scripts 40/59_alu/60), NOT the superseded `03_methylation`
(scripts 16/19). The 45S cohort needs BOTH server halves. Tools: modkit, minimap2, fibertools-rs.
