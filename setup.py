#!/usr/bin/env python

from distutils.core import setup

setup(name='SequenceEnrichmentSimulator',
  version='0.1',
  description='A framework to conduct Sequencing Enrichment simulations and related bioinformatics approaches.',
  long_description='A computational framework to simulate various staeps in sequencing enrichment and mining of NGS reads for sequences of interest. A set of pure python modules are representing the various precedures, a run script is provided to conduct the experiments and a control parameters file to set the experimental parameters. We intent to use the framewrok as a tool for the general plant (and not only) community to enchance the studies of sequence enrichment approaches.',
  requires=['Bio'],
#  packages=['seqEnrichSim'],
  py_modules=['seqEnrichSim'],
  scripts=['runSeqEnrichmentExperiment', 'conductSimulationExperiments.sh'],
  author='Costas Bouyioukos',
  author_email='costas.bouyioukos@tsl.ac.uk',
  url='http://sites.google.com/site/cbouyio/',
  license='GNU General Public License (GPL) GPLv3 (or newer)',
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Environment :: Console',
    'Intended Audience :: Education',
    'Intended Audience :: Science/Research',
    'License :: Free for non-commercial use',
    'License :: OSI Approved :: GNU General Public License (GPL)',
    'Operating System :: POSIX',
    'Operating System :: POSIX :: Linux',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2.7',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering ::Bio-Informatics',
    'Topic :: Scientific/Engineering ::Computational Biology',
    'Topic :: Scientific/Engineering ::Regulatory Gene Networks',
  ],
)

