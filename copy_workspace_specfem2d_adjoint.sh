#!/bin/bash

indir=$1
outdir=$2

mkdir -p $outdir/DATA
mkdir -p $outdir/OUTPUT_FILES
mkdir -p $outdir/SEM
cp $indir/interfaces_membrane.dat $outdir
cp $indir/run_this_example.sh $outdir
cp $indir/run_this_example_mpi.sh $outdir
cp $indir/submit_janus_run.sh $outdir
cp $indir/DATA/* $outdir/DATA
cp $indir/SEM/* $outdir/SEM
cp $indir/OUTPUT_FILES/absorb* $outdir/OUTPUT_FILES
cp $indir/OUTPUT_FILES/lastframe* $outdir/OUTPUT_FILES
