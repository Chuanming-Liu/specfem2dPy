#!/bin/bash

indir=$1
outdir=$2

mkdir -p $outdir/DATA
cp $indir/interfaces_membrane.dat $outdir
cp $indir/run_this_example.sh $outdir
cp $indir/run_this_example_mpi.sh $outdir
cp $indir/submit_janus_run.sh $outdir
cp $indir/DATA/* $outdir/DATA
