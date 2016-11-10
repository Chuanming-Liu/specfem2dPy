#!/bin/bash

module load gcc
rm -rf xadj_seismogram_t xadj_seismogram_amp 
gfortran adj_seismogram_t.f90 -o xadj_seismogram_t
gfortran adj_seismogram_amp.f90 -o xadj_seismogram_amp
