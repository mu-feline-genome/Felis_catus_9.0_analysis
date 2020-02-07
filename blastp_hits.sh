#!/bin/bash


CWD=/home/buckley/Documents/projects/project_files/find_orthologs/Felis_catus-Homo_sapian
SP1name=Felis_catus
SP2name=Homo_sapiens
SP1file=/home/buckley/Documents/projects/reference_files/fasta/Felis_catus/ensembl_97/Felis_catus.Felis_catus_9.0.pep.all.fa
SP2file=/home/buckley/Documents/projects/reference_files/fasta/Homo_sapians/ensembl_97/Homo_sapiens.GRCh38.pep.all.fa


mkdir -p $CWD

makeblastdb -dbtype prot -out $CWD/$SP1name.db -in $SP1file &
makeblastdb -dbtype prot -out $CWD/$SP2name.db -in $SP2file &

wait

blastp -query $SP1file -db $CWD/$SP2name.db -out $CWD/${SP1name}-${SP2name}.out -outfmt 6 -culling_limit 3 -num_threads 4 &

blastp -query $SP2file -db $CWD/$SP1name.db -out $CWD/${SP2name}-${SP1name}.out -outfmt 6 -culling_limit 3 -num_threads 4 &



