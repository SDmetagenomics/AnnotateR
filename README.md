# AnnotateR
Takes and annotates raw metagenome bins for input into AnnoVisR

This script requires the following programs to be in your path:
  1. Prodigal
  2. HMMER
  3. CheckM (Optional; https://github.com/Ecogenomics/CheckM/wiki/Installation)
  4. pplacer (If running CheckM; https://github.com/matsen/pplacer/releases/tag/v1.1.alpha17)

Input Requirements:

  - Input files must be in fasta format and can either be genome DNA scaffolds or protein sequences
  - DNA scaffold files must be named [my_genome_name].fasta
  - Protein files must be named [my_proteins_name].faa
  
