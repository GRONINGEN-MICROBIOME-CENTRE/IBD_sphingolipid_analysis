#!/usr/bin/env bash
# DIAMOND 2.1.8; UniRef90 database prepared separately.
diamond blastp \
  --query input/enzyme_queries.fasta --db reference/uniref90.dmnd \
  --out output/diamond_hits.tsv --outfmt 6 --evalue 1e-5 \
  --max-target-seqs 500 --id 70 --threads 8

# BLAST+ 2.14.1; apply per MAG protein set.
makeblastdb -in input/mag_proteins.faa -dbtype prot -out reference/mag_proteins
blastp -query input/enzyme_queries.fasta -db reference/mag_proteins \
  -out output/blastp_hits.tsv -evalue 1e-5 -num_threads 2 \
  -max_target_seqs 10 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
# Keep highest-bit-score retained hit per MAG target protein downstream.
