#!/bin/sh

# Get the cwd
pwd=`pwd`;
cwd=`basename ${pwd}`;

# Initialise the random seed corrector.
r=1;

# Submit to LSF
# Get all the command line arguments as the input filenames.
for f in $@;
do
  bsub -o ${cwd}_LSF.out -x -q 16gb runSeqEnrichmentExperiment -p seqEnrichSim.cpf -i ${f} -r ${r} ;
  r=`expr $r + 1`
done

