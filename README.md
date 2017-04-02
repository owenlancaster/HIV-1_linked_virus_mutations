## Determine HIV-1 mutation frequency (overall and linked mutations) from raw single reads from MinION
---
#### Need to convert MinION reads from fast5 to fasta format for the sample (this example uses [poretools](https://poretools.readthedocs.io/en/latest/)):
```
poretools fasta DIRECTORY_FOR_FAST5 > SAMPLE_NAME.fasta
```
---
#### Run linked mutation python script to count/calculate the total number/frequency of mutations for each virus particle for the sample:

##### Example to run linked mutation script on fasta file, BLASTx will be initially run on this be first of all run. Specifying -l enables calculating linked mutations on each BLAST hit (as well as overall single mutation counts) 

```
python HIV-1_linked_mutations_blastx.py -f reads.fasta -e 0.01 -l
```

##### N.B. this step requires NCBI BLAST+ to be installed, and uses the "all" database (included in this respository) which is a HIV-1 RT/PR/INT pre-formatted protein BLAST database created from the Stanford HIV-1 subtype B protein sequences. The output format from the search in BLAST xml format (outfmt 5) which is then used in the script

##### N.B. if you already have the BLAST xml output file and do not want BLAST to be re-run on the fasta file then use the -x option instead of -f and specify the BLAST xml file instead of the fasta file
---
To get help run:

```
python HIV-1_linked_mutations_blastx.py -h
```
---
##### TODO: example other options for the script e.g. specifying dictionaries of mutations