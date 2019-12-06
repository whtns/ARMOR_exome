# to access a snakemake pipeline:
1. enter project directory
2. cd to ARMOR directory
3. `git pull origin master` to be sure you have most update snakemake pipeline
4. `conda activate ARMOR`
5. `snakemake -n --use-conda` for a 'dry run'
6. `snakemake -j 6 --use-conda` for full run

* try running all of this inside a tmux sesssion 

# to view all of the rules for a snakemake pipeline:

`snakemake --rulegraph | dot -Tpdf > rulegraph.pdf`
