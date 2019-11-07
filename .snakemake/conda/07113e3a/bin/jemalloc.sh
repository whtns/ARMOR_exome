#!/bin/sh

prefix=/dataVolume/storage/single_cell_projects/sc_RB_devel/FACS_20190501_sunlee_H_sapiens_proj/ARMOR/.snakemake/conda/07113e3a
exec_prefix=/dataVolume/storage/single_cell_projects/sc_RB_devel/FACS_20190501_sunlee_H_sapiens_proj/ARMOR/.snakemake/conda/07113e3a
libdir=${exec_prefix}/lib

LD_PRELOAD=${libdir}/libjemalloc.so.2
export LD_PRELOAD
exec "$@"
