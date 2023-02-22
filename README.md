# getPIP : find plant-inducible promoter (PIP) boxes and the closest downstream gene.

PIP boxes were first described by [Fenselau et al](https://doi.org/10.1094/mpmi-8-0845) in 1995 in *Xanthomonas campestris* pv. *vesicatoria* as conserved sequences upstream of hrpB genes, and were later found to be activated by the HrpX regulator ([Koebnik et al 2006](https://doi.org/10.1128/JB.00795-06)). PIP boxes are involved in the activation of five of the six structural type 3 secretion system operons and multiple type 3 effectors. A consensus sequence of PIP boxes was defined as: TTCGB-N<sub>15</sub>-TTCGB, where B is any nucleotide except A. "Imperfect" PIP boxes can also be found with the sequence: TTCGB-N<sub>8</sub>-TTCGB.

Investigating the presence of such PIP boxes upstream of genes can be helpful in identifying genes involved in virulence. getPIP will look for PIP boxes and tell you which is the closest downstream gene.

# Dependencies
This script is written in Python 3 and requires the [BioPython](https://biopython.org/wiki/Download)  package to be installed.

# Usage

```
getPIP.py -f <fasta_file> -a <annotation_file> [-o <output_csv>] [--verbose]
-f, --fasta: Fasta file path
-a, --annotation: GFF Annotation file path
-o, --output: Output CSV file path (optional)
-v, --verbose: Enable verbose output (optional)
-h, --help: Show this help message
```

The script takes as input a genome sequence in FASTA format and a genome annotation file in GFF format. The output is a CSV file containing information on each motif match: the motif sequence, position in the genome, closest gene, gene location, strand, and distance to the gene. If you do not specify an output name, a file called <your fasta file name>.csv will be created in your current directory by default.
