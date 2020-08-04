#!/bin/bash

polyftsdir=/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release

IDIR=`pwd`

source activate mp

python ./runner_general.py -PolyFTS ${polyftsdir} -stat DGC -p DIS LAM -nt 1 1 -sp $IDIR/SEEDS/ -ds 0.1 -cl AB-Bottlebrush -dt __dt__ -nref 1 -ss 0.001 -fs 1.0 0.1 -L0 5.0 5.0 -npw 512 512 -ftol 1e-7 -stol 1e-8 -bsp 1 2 -sas 1 2 -sac 1.0 1.0 -space 1 1 -nbbmin __nbbmin__ -nbbmax __nbbmax__ -dnbb __dnbb__ -nscmin __nscmin__ __nscmin__ -nscmax __nscmax__ __nscmax__ -ndsc __dnsc__ __dnsc__ -chimin __chimin__ -chimax __chimax__ -dchi __dchi__ -dirs nbb nsc f chi -dirnl None 0 0 0 -runpoly False -subpoly False -chain False

