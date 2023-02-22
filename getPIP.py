from Bio import SeqIO
import re, csv, os
import sys, getopt

def print_help():
    print("Usage: getPIP.py -f <fasta_file> -a <annotation_file> [-o <output_csv>] [--verbose]")
    print("-f, --fasta: Fasta file path")
    print("-a, --annotation: GFF Annotation file path")
    print("-o, --output: Output CSV file path (optional)")
    print("-v, --verbose: Enable verbose output (optional)")
    print("-h, --help: Show this help message")

fasta_file = None
annotation_file = None
VERBOSE = False
output_file = None

# read command line options
if not any(i in sys.argv for i in ('-a','--annotation','-f','--fasta','-h','--help')) : #if no input argument specified, print help and exit
    print_help()
    sys.exit()

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:a:o:v", ["help", "fasta=", "annotation=", "output=", "verbose"])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)

# handle command line options
    for opt, arg in opts:
        if opt == '-h' or opt == '--help':
            print_help()
            sys.exit()
        if opt in ("-f", "--fasta"):
            fasta_file = arg
        elif opt in ("-a", "--annotation"):
            annotation_file = arg
        elif opt in ("-o", "--output"):
            output_file = arg
        elif opt in ("-v", "--verbose"):
            VERBOSE = True

# check required options
if fasta_file is None:
    print("Fasta file is required.")
    sys.exit(2)
if annotation_file is None:
    print("Annotation file is required.")
    sys.exit(2)
if output_file is None:
    output_file = os.path.splitext(os.path.basename(fasta_file))[0] + ".csv"

#get genome sequence and reverse genome sequence as strings
genome_sequence = ""
reverse_genome_sequence = ""
for record in SeqIO.parse(fasta_file, "fasta"):
    genome_sequence += str(record.seq)
    reverse_genome_sequence += str(record.seq.reverse_complement())
    
# Read in the genome annotation file and store the annotations as a list of dictionaries
annotations = []
with open(annotation_file) as f:
    for line in f:
        if not line.startswith("#"):
            fields = line.strip().split("\t")
            if len(fields) > 5:
                seqid = fields[0]
                source = fields[1]
                feature_type = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                score = fields[5]
                strand = fields[6]
                phase = fields[7]
                attributes = fields[8]
                
                # Parse attributes to get gene ID
                gene_id = None
                if feature_type == "gene":
                    gene_id_match = re.search("ID=([^;\n]+)", attributes)
                    if gene_id_match is not None:
                        gene_id = gene_id_match.group(1)
                
                annotation = {
                    "seqid": seqid,
                    "source": source,
                    "type": feature_type,
                    "start": start,
                    "end": end,
                    "score": score,
                    "strand": strand,
                    "phase": phase,
                    "name": attributes,
                    "gene_id": gene_id
                }
            annotations.append(annotation)

# Search for PIP motifs in the genome sequence
PIP_motifs = ["TTCG[^A]...............TTCG[^A]", "TTCG[^A]........TTCG[^A]"]
all_motif_matches = []
for motif in PIP_motifs:
    # Forward strand matches
    motif_matches = [(match.start(), match.group(), "+") for match in re.finditer(motif, genome_sequence)]
    all_motif_matches.extend(motif_matches)
    
    # Reverse strand matches
    motif_matches = [(match.start() - len(motif), match.group(), "-") for match in re.finditer(motif, reverse_genome_sequence)]
    all_motif_matches.extend(motif_matches)
    
# Write output to CSV file
with open(output_file, 'w', newline='') as csvfile:
    fieldnames = ['Motif Sequence', 'Motif Position', 'Closest Gene', 'Gene Location', 'Strand', 'Distance to Gene']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    # Find the closest downstream gene for each motif match and write to CSV
    for motif_match in all_motif_matches:
        closest_gene = None
        min_distance = float("inf")
        for annotation in annotations:
            if annotation["type"] == "gene":
                if motif_match[2] == "+" and annotation["strand"] == "+":
                    gene_distance = annotation["start"] - motif_match[0]
                elif annotation["strand"] == "-":
                    gene_distance = motif_match[0] - annotation["end"]
                if gene_distance > 0 and gene_distance < min_distance:
                    closest_gene = annotation
                    min_distance = gene_distance
        if closest_gene is not None and min_distance < 1000:
            writer.writerow({'Motif Sequence': motif_match[1][:5] + '-N' + str(len(motif_match[1])-10) + '-' + motif_match[1][-5:], 'Motif Position': motif_match[0], 'Closest Gene': closest_gene['type'] + ' ' + closest_gene['gene_id'], 'Gene Location': closest_gene['seqid'] + ':' + str(closest_gene['start']) + '-' + str(closest_gene['end']), 'Strand': closest_gene['strand'], 'Distance to Gene': min_distance})
            if VERBOSE is True:
                print(f"Motif {motif_match[1][:5]}-N{len(motif_match[1])-10}-{motif_match[1][-5:]} found at position {motif_match[0]}, the closest gene is {closest_gene['gene_id']} on {closest_gene['seqid']} strand {closest_gene['strand']}, {min_distance:,} bp downstream of the motif")
    print(f"Finished! Results written in {output_file}")
