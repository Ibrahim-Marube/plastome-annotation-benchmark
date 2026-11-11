# plastome-annotation-benchmark

How well does reference-based plastome annotation work, and what limits it?
Built while learning genome annotation (MSc Bioinformatics, University of Potsdam);
all runs were done on the university HPC cluster with a Snakemake workflow.

**Question.** PGA (Plastid Genome Annotator) annotates a plastome by transferring genes
from a related reference. Does accuracy depend on how closely related the reference is,
and what else goes wrong?

**Short answer.** On this target the choice of reference species mattered little
(about 1 gene in 109). What mattered most was whether the reference's *gene names*
follow the convention PGA expects: with RefSeq references as downloaded, PGA silently
skipped every rRNA and almost every tRNA (66-70% of genes found); after renaming the
genes in the reference, 93.6-94.5% were found.

## Data

| role | species | accession | note |
|---|---|---|---|
| target (answer key) | *Medicago truncatula* | NC_003119.8 | 124,033 bp; RefSeq has 76 CDS, 29 tRNA, 4 rRNA = 109 features. The sequence is the input, the RefSeq annotation is the truth. |
| reference | *Lotus japonicus* | NC_002694.1 | same family (Fabaceae) |
| reference | *Arabidopsis thaliana* | NC_000932.1 | another eudicot |
| reference | *Amborella trichopoda* | NC_005086.1 | outside the other angiosperms |

## Method

- PGA at commit `0683a0e`, BLAST+ 2.13.0 (module `bio/BLAST+/2.13.0-gompi-2022a`),
  Snakemake 9.27.1. Each reference is run twice: **as downloaded** and **normalized**
  (gene names rewritten, see below).
- A feature is identified by type (CDS, tRNA, rRNA), strand and 3' end. *Found* = same
  strand and 3' end; *exact* = also the same 5' end. Sensitivity = found / truth,
  precision = found / predicted. A second, position-tolerant measure (**+-5 nt**) also
  counts one-to-one matches whose start and stop are each within 5 bases, so a boundary
  slip is not scored as one miss plus one false positive.

## Results

| reference | names | CDS found | tRNA found | tRNA found (+-5 nt) | rRNA found | all: sensitivity | all: sensitivity (+-5 nt) | all: precision | all: exact start |
|---|---|---|---|---|---|---|---|---|---|
| Lotus japonicus | as downloaded | 75/76 | 1/29 | 1/29 | 0/4 | 69.7% | 69.7% | 91.6% | 68/109 |
| Lotus japonicus | normalized | 75/76 | 24/29 | 27/29 | 4/4 | 94.5% | 97.2% | 96.3% | 91/109 |
| Arabidopsis thaliana | as downloaded | 75/76 | 1/29 | 1/29 | 0/4 | 69.7% | 69.7% | 89.4% | 66/109 |
| Arabidopsis thaliana | normalized | 75/76 | 24/29 | 27/29 | 3/4 | 93.6% | 96.3% | 92.7% | 88/109 |
| Amborella trichopoda | as downloaded | 72/76 | 0/29 | 0/29 | 0/4 | 66.1% | 66.1% | 96.0% | 62/109 |
| Amborella trichopoda | normalized | 73/76 | 25/29 | 28/29 | 4/4 | 93.6% | 96.3% | 93.6% | 88/109 |

(`results/comparison.md`; per-run score tables, missed/extra gene lists and the PGA
output are in `results/`.)

## Findings

**1. PGA silently skips genes whose names do not match its conventions (confirmed).**
The PGA source only takes the rRNA branch for names like `rrn16`, while RefSeq writes
`rrn16S` or `rrn16s`; its no-intron tRNA branch requires an anticodon in the name
(`trnQ-UUG`), but RefSeq often has `trnQ`, and Amborella writes `tRNA-His (GUG)`. Two
names, `trnR` and `trnA`, are special-cased in the source, which is why exactly one tRNA
(`trnR`) was annotated with the Lotus reference. `workflow/scripts/normalize_reference.py`
fixes the names: it strips the trailing S/s from rRNA names and, for bare tRNA names,
reads the anticodon from the spliced tRNA sequence and checks it against `/product`
(anticodon -> codon -> amino acid). For Lotus all 35 tRNAs verified; for Arabidopsis 34
of 35 did, and the exception (an unusual `trnS`) is flagged in
`results/reference_reports/`, not forced. Amborella already gives the anticodon in the
name, so those names were only reformatted.

**2. Relatedness had a small effect once names were fixed.** Sensitivity: 94.5% (Lotus),
93.6% (Arabidopsis), 93.6% (Amborella). The closest relative gave the best precision
(96.3% vs 92.7% and 93.6%) and the most exact start sites (91 vs 88 and 88). This is one
target genome, so I did not test whether these differences are meaningful.

**3. Most remaining tRNA "errors" are boundary slips (confirmed for Lotus).** `trnC-GCA`,
`trnF-GAA` and `trnN-GUU` are both "missed" and "extra" because PGA's ends differ from
RefSeq's by 1-3 bases (e.g. `trnC` 90806..90879 in RefSeq, 90806..90876 from PGA). With
the +-5 nt tolerance 27-28 of 29 tRNAs are found.

**4. Gene-level differences that remain (named, mostly not yet explained).**
- `ycf68` is missed with all three references (it lies within the span of the `trnI-GAU` gene).
- `ycf4` is predicted with all three references but absent from the RefSeq record
  (which also has no `rpl22` or `rps16`). Whether it is a real gene RefSeq omits is not
  checked here.
- Arabidopsis and Amborella additionally predict `rpl22` (not in RefSeq) and an `rps12`
  that disagrees with RefSeq's `rps12` (a trans-spliced gene; likely a structure
  difference, not verified).
- Amborella misses `petL` and `psaI`, both very short genes.
- Arabidopsis' RefSeq record has no 5S rRNA gene, so 5S cannot be transferred (3/4 rRNAs).

**5. A bug in PGA made about 1 run in 8 silently produce an unannotated file (confirmed).**
PGA builds temp-file names from 14 random characters; the character list contains `$`.
The names go straight into shell commands, so `$$` becomes a process ID and `$abc`
expands to nothing, and BLAST is asked for a file name different from the one PGA wrote.
The run still "succeeds". Stock PGA failed 4 of 30 times in my test; a copy with `$`
removed from that one line failed 0 of 40 times. The workflow patches a *copy*
(`workflow/scripts/patch_pga.sh`, which checks that exactly one line changed) and runs PGA
through `run_pga.sh`, which accepts a result only if it has at least 30 CDS and otherwise
retries. The scorer also refuses to score an unannotated file.

**What did not matter (exploratory runs, not in the workflow).** `-f linear` and
`-i 100000` (the IR-length threshold, relevant because *Medicago* lacks an inverted repeat)
gave results identical to the default; the BLAST+ 2.13.0 module reproduced my first runs.

## Limitations

- One target genome, one annotation tool, three references: no statistics, and I cannot
  tell whether the (small) reference effect in finding 2 would hold for other targets.
- RefSeq is the answer key but is not perfect truth (`ycf4`, `rpl22`).
- The +-5 nt tolerance is a choice; strict numbers are reported alongside.
- tRNA anticodons are inferred from sequence with a heuristic; one case is unresolved.
- The PGA fix is applied to a local copy and not submitted upstream.
- Reproduction was checked on one cluster only: a fresh clone with the committed
  `results/` removed re-downloaded the data and gave identical score files and
  comparison table (about 1 minute, 2 cores).

## Reproduce (on the cluster)

~~~bash
git clone https://github.com/Ibrahim-Marube/plastome-annotation-benchmark.git
cd plastome-annotation-benchmark
module load lang/Miniconda3/23.9.0-0 && eval "$(conda shell.bash hook)"
conda activate <env with snakemake 9>
snakemake -s workflow/Snakefile --cores 2 --use-envmodules
cat results/comparison.md
~~~

The workflow downloads the four GenBank records from NCBI, clones PGA at the pinned
commit, patches a copy, normalizes the references, runs PGA (6 runs: 3 references x
as downloaded / normalized, with retries), scores them and builds the table.

## Layout

- `config.yaml`: accessions, PGA commit, retry settings
- `workflow/Snakefile`: the pipeline
- `workflow/scripts/`: `compare_gb.py` (scoring), `normalize_reference.py`,
  `patch_pga.sh`, `run_pga.sh`, `gb_to_fasta.py`, `make_table.py`
- `results/`: `comparison.md`, `score/`, `pga/`, `reference_reports/`

PGA: Qu et al. 2019, *Plant Methods* 15:50.
