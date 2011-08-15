# Makefile to install the modules and cripts and generate the documentation.

PYTHONMODULES = seqEnrichSim.py

html : ${PYTHONMODULES}
	rm -rf html
	epydoc -v --html ${PYTHONMODULES}

install : ${PYTHONMODULES}
	python setup.py install --home=${HOME}

clean :
	rm -rf velvetAssembly* html build *.pyc *SeqHomol.fasta *IlluminaReads.fasta* assemblyStats.txt bestAssemblyContigs.fasta* *.nsq *.nhr *.nin

.PHONY : install clean
