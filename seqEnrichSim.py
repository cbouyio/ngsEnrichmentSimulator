#!/usr/bin/env python


"""
Python module to conduct second generation sequencing simulators.

@author: Costas Bouyioukos
@organization: The Sainsbury Laboratory
@since: Novemeber 2011
@copyright: The program is coming as it is. You have the right to redistribute,
transform and change the source code presuming the apropriate reference and
the lisence is kept free.
@license: GNU GPL3 or newer.
@contact: U{Costas Bouyioukos<mailto:k.bouyioukos@uea.ac.uk>}
@version: 0.0.1
"""


import re
import random
import string
import subprocess
import shlex
import shutil
#import copy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import FastxMetrics


class AbstactReference(object) :
  """Base to represent a reference sequence.

  Abstract class to represnt a genome, a chromosome or parts of genomic
  sequences, a transcriptome.
  """

  def __init__(self) :
    """Empty constructor, raise an exception, each subclass should be contructed
    individualy.
    """
    raise StandardError, 'cannot instantiate AbstactReference class'



class GenomeSequence(AbstactReference) :
  """Class to represent a genome derived collection of sequences.

  @ivar refGenerator: An iterator wich yields the next fasta record in each
  invocation.
  @type refGenerator: C{class 'Generator'}
  """

  def __init__(self, refFile, refFileType) :
    """Constructor takes a fatsa file with the reference Genome sequence.

    Constructs a Biopyhon file parser object.
    @param refFile: An open ready to read Python file object or the name of the
    file  containing the reference genomic sequences.
    @type refFile: C{file} or C{string}
    @param refFileType: A field to designate the type of the reference sequence
    file. Only the "fasta"
    file format is supported at the time.
    @type refFileType: C{string}
    """
    if refFileType == "fasta":
      self.refGenerator = SeqIO.parse(refFile, refFileType, IUPAC.unambiguous_dna)
    else :
      raise StandardError, 'Only fasta format is supporting by the module at the moment...'



class TranscirptomeSequence(AbstactReference) :
  """Class to represent a trnascriptome derived collection of sequences.
  """

  def __init__(self, refFasta) :
    """Constructor takes a fatsa file with the reference Transcriptome sequence

    """

  def next(self) :
    """Return the next fasta sequence

    """



class Predator(object) :
  """Superclass to represent sequences of interest.

  These are the sequences that will be generating the probes (baits) on the
  capture array for sequencing enrichment.
  E.g. genes of interest, protein families of interest.
  """

  def __init__(self, targetFastaFile) :
    """Constructor.
    """



class Experiment(object) :
  """A class to represent an individual experiment.

  """

  def __init__(self) :
    """Constructor.

    """



class CaptureArray(object) :
  """Class to represent a capture array.

  """

  def __init__(self) :
    """Constructor.

    """



class ParametersParser(object) :
  """Class to represent a parameters file parser.

  Parse the control paremeters file and constructs the appropriate parameter
  objects.
  The parameters names are appearing as class variables. Any change to the
  paramters file format should be reflected here too.
  """

  # Class variables to specify the control parameters file attributes.
  magic = '#sequenceEnrichmentParameters'
  parametersIdentifiers = ["experimentType", "LibraryParamaters"]
  experimentTypeValues = ["test", "partialKnowledge"]

  def __init__(self, parametersFile) :
    """Constructor

    The constructor checks for the magic line in the parameters file.
    @param parametersFile: Text file containing the control parameters of the
    experiment. Either a filehandler or a string with the filename.
    @type parametersFile: C{file} or C{str}
    """
    if isinstance(parametersFile, str) :
      fh = open(parametersFile, 'r')
    elif isinstance(parametersFile, file) :
      fh = parametersFile
    else :
      raise StandardError, 'Please specify a filename or an open filehandler as a parameters file argument.'
    self.fh = fh


  def parse(self) :
    """Parse the parameters file and populate the instance variables of the
    relevant parameters objects.
    """

    def parse_next_parameter_block_header(paramBlock) :
      """Nested function to parse the next parameter block.

      @param paramBlock: A line from the control paramters file.
      @type paramBlock: C{'str'}
      """
      line = self.fh.readline().strip()
      if not re.match(paramBlock, line) :
        raise StandardError, 'Expected "%s", got "%s"' % (paramBlock, line)
      if paramBlock == 'randomSeed' :
        return int(line.split(':')[1])
      elif paramBlock == '#experimentType' :
        return self.parse_experiment_type_parameters()
      elif paramBlock == '#LibraryParamaters' :
        return self.parse_library_parameters()
      elif paramBlock == '#NGSparameters' :
        return self.parse_ngs_parameters()
      elif paramBlock == '#SequenceHomologyParameters' :
        return self.parse_seq_homol_params()
      elif paramBlock == '#AssemblyParameters' :
        return self.parse_assembly_params()
      else :
        raise StandardError, 'parameter block header "%s" is not recognised as a seqEnrichSim parameter block' % paramBlock

    # The implementation of the parser.
    line = self.fh.readline().strip()
    if line == '' :
      raise StandardError, 'expected magic but got EOF'
    if line != self.magic :
      raise StandardError, 'File is not starting with magic line "%s", it is not a valid sequence enrichment simulator control parameter file' % self.magic
    rndSeed = parse_next_parameter_block_header('randomSeed')
    expP    = parse_next_parameter_block_header('#experimentType')
    libP    = parse_next_parameter_block_header('#LibraryParamaters')
    ngsP    = parse_next_parameter_block_header('#NGSparameters')
    seqHomP = parse_next_parameter_block_header('#SequenceHomologyParameters')
    assP    = parse_next_parameter_block_header('#AssemblyParameters')
    # Append the random seed to Parameters subclasses.
    libP.append(rndSeed)
    ngsP.append(rndSeed)
    return Parameters(rndSeed, ExperimentParameters(expP), LibraryParameters(libP), NGSParameters(ngsP), SeqHomologyParameters(seqHomP), AssemblyParameters(assP))


  def parse_experiment_type_parameters(self) :
    """Parse the experiment type parameters section of the parametrs file.

    """
    experimetTypeParamsList = []
    line = self.fh.readline().strip()
    if line not in self.experimentTypeValues :
      raise StandardError, 'Exteriment type not supported. Please specify one of the %s' % self.experimentTypeValues
    experimetTypeParamsList.append(line)
    return experimetTypeParamsList


  def parse_library_parameters(self) :
    """Parse the library control parameters section of the parameters file.

    """
    libraryParamtersList = []
    # Go over the block line by line
    lt = self.parse_next_name_value('libraryType')
    libraryParamtersList.append(str(lt))
    ins = self.parse_next_name_value('insertSize')
    libraryParamtersList.append(int(ins))
    sd = self.parse_next_name_value('standardDeviation')
    libraryParamtersList.append(int(sd))
    cv = self.parse_next_name_value('coverage')
    libraryParamtersList.append(int(cv))
    return libraryParamtersList


  def parse_ngs_parameters(self) :
    """Parse the NGS parameters section of the control parameters file.

    """
    ngsParameters = []
    sp = self.parse_next_name_value('sequencingPlatform')
    ngsParameters.append(str(sp))
    rl = self.parse_next_name_value('readLength')
    ngsParameters.append(int(rl))
    pe = self.parse_next_name_value('pairedEnd')
    if pe == 'True' :
      ngsParameters.append(True)
    elif pe == 'False' :
      ngsParameters.append(False)
    else :
      raise StandardError, 'Value "%s" not supported. Specify one of "True" or "False" for "pairedEnd" field.' % pe
    em = self.parse_next_name_value('errorModel')
    if em == 'True' :
      ngsParameters.append(True)
    elif em == 'False' :
      ngsParameters.append(False)
    else :
      raise StandardError, 'Value "%s" not supported. Specify one of "True" or "False" for "erroModel" field.' % em
    return ngsParameters


  def parse_seq_homol_params(self) :
    """Parse the Sequence Homology Parameters section of the control parameters
    file.

    """
    seqhParams = []
    hps = self.parse_next_name_value('hmmProfileFiles')
    hpsl = []
    for hp in hps.split(';') :
      hpsl.append(hp)
    seqhParams.append(hpsl)
    evs = self.parse_next_name_value('hmmEvalues')
    evsl = []
    for ev in evs.split(';') :
      evsl.append(ev)
    seqhParams.append(evsl)
    blastDB = self.parse_next_name_value('BLASTdatabase')
    if blastDB == 'None' :
      blastDB = None
    seqhParams.append(blastDB)
    si = self.parse_next_name_value('seqIdentity')
    seqhParams.append(int(si))
    sl = self.parse_next_name_value('seqLengthAligned')
    seqhParams.append(int(sl))
    return seqhParams


  def parse_assembly_params(self) :
    """Parse the Assembly Parameters block of the control parameters file.

    """
    assParams = []
    ass = self.parse_next_name_value('assembler')
    assParams.append(ass)
    kmers = self.parse_next_name_value('kmerSize')
    kml = []
    for km in kmers.split(';') :
      kml.append(int(km))
    assParams.append(kml)
    assRef = self.parse_next_name_value('referenceAssembly')
    assParams.append(assRef)
    return assParams


  def parse_next_name_value(self, name) :
    """Check a name-colon-value pair line for the correct existence of name,
    return a list with the name value pair.

    @param name: Specify the name of the name:value pair.
    @type name: C{'str'}
    @rtype: C{'str'}
    """
    line = self.fh.readline().strip()
    if not re.match(name, line) :
      raise StandardError, 'Expected "%s" got "%s"' % (name, line.split[0])
    return line.split(':')[1].strip()



class Parameters(object) :
  """Superclass of the paramters objects.

  Keeps the random number generator object of the experiment.
  Implements a str method to print out parameter name:value pairs.
  """

  def __init__(self, rndSeed, expParams, libraryParams, ngsParams, seqHomolParams, assParams) :
    """Constructor.

    The class implements a print parameters method.
    """
    self.rndSeed            = rndSeed
    self.expParameters      = expParams
    self.libraryParameters  = libraryParams
    self.ngsParameters      = ngsParams
    self.seqHomolParameters = seqHomolParams
    self.assemblyParameters = assParams


  def __str__(self) :
    """Method to print the name:value parameter pairs.

    """
    pass



class ExperimentParameters(Parameters) :
  """Class to represent the experiment type parameters.

  """

  def __init__(self, expParamsList) :
    """Constructor
    """
    self.experimentType = expParamsList[0]



class LibraryParameters(Parameters) :
  """Class to represent the library generation parameters.

  The expected library parameters are represented by a list as a class variable.
  The order of the list specifies the order the parameters should be found in
  the parameters file record.
  """

  libraryParameterNames = ["libraryType", "insertSize", "standardDeviation", "coverage"]

  def __init__(self, libraryParametersList) :
    """Constructor.

    @param libraryParametersList: A list of the parameter values from the
    parameter file. The paremeter file parses shouls take care for the
    consistency of this list.
    @type libraryParametersList: C{list}
    """
    # Some sanity checks (the + 1 is required as the random seed is also passed
    # to all the suclasses of parameters.
    if len(self.libraryParameterNames) + 1 != len(libraryParametersList) :
      raise StandardError, "Library parameters list not the same size to library parameters names list."
    self.libraryType       = libraryParametersList[0]
    self.insertSize        = libraryParametersList[1]
    self.standardDeviation = libraryParametersList[2]
    self.coverage          = libraryParametersList[3]
    self.rndSeed           = libraryParametersList[4]



class NGSParameters(Parameters) :
  """Class to represent the next generation sequencing experiment parameters.
  """

  NGSParametersList = ["", "", "", ""]

  def __init__(self, ngsParamsList) :
    """Constructor.

    @param ngsParamsList: A list with the NGS params files as they have been
    parsed from the ParametersParser class.
    @type ngsParamsList: C{'list'}
    """
    # A simple sanity check.
    if len(ngsParamsList) != len(self.NGSParametersList) + 1 :
      raise StandardError, "NGS parameters list not the same size to NGS parameters names list."
    self.sequencingPlatform = ngsParamsList[0]
    self.readLength         = ngsParamsList[1]
    self.PE                 = ngsParamsList[2]
    self.errorModel         = ngsParamsList[3]
    self.rndSeed            = ngsParamsList[4]



class SeqHomologyParameters(Parameters) :
  """Class to contain the sequence homology parameters.

  """

  def __init__(self, seqHomParamList) :
    """The constructor.

    """
    self.HMMProfiles    = seqHomParamList[0]
    self.HMMEvalues     = seqHomParamList[1]
    self.blastDatabase  = seqHomParamList[2]
    self.seqIdentity    = seqHomParamList[3]
    self.seqAlignLength = seqHomParamList[4]



class AssemblyParameters(Parameters) :
  """Container class of the assembly program and reference control parameters.

  """

  def __init__(self, assParamList) :
    """The constructor...

    """
    self.assemblyProgram       = assParamList[0]
    self.assemblyKmerSizeList  = assParamList[1]
    self.assemblyReferenceFile = assParamList[2]



class LibraryFactory(object) :
  """Class to represent and specify the generation of a sequencing library.

  @ivar sequencingLibraryType: Text designating the sequencing library type
  @type sequencingLibraryType: C{string}
  @ivar libraryInsertSize: The insert size of the library.
  @type libraryInsertSize: C{int}
  @ivar libraryInsertSD: The standard deviation of the fragment sizes of the
  library
  @type libraryInsertSD: C{float}
  @ivar libraryCoverage: The deapth of the library coverage.
  @type libraryCoverage: C{int}
  """

  def __init__(self, libraryParameters, rndSeedCorrector = 0) :
    """Constructor.

    Takes the LibraryParameters object as a parameter.
    @param libraryParameters: A C{LibraryParameters} instance to describe the
    parameters of the library generation.
    @type libraryParameters: C{class LibraryParameters}
    @param rndSeedCorrector: A number to be added to the seed of the random
    number generator to avoid repetition of the same randfom number stream when
    multiple sequence fiels are provided.
    @type rndSeedCorrector: C{int}
    """
    if not isinstance(libraryParameters, LibraryParameters) :
      raise StandardError, 'The LibrayFactory object requires a valid LibraryParameters instance to b instantiated.'
    self.sequencingLibraryType = libraryParameters.libraryType
    self.libraryInsertSize     = libraryParameters.insertSize
    self.libraryInsertSD       = libraryParameters.standardDeviation
    self.libraryCoverage       = libraryParameters.coverage
    self.rndSeed               = libraryParameters.rndSeed
    self.rndSeedCorrector      = rndSeedCorrector


  def __call__(self, refFastaFile) :
    """Caller to generate a sequencing library.

    Factory method that generates the clones of a library.
    Returns a Bio.SeqIO.SeqRecord with all the library clones.
    @param refFastaFile: A fasta file containing the reference genome. Either a
    multi-fasta files or single fasta entry files are supported. The file can
    either come as a filehandler or as a filename string.
    @type refFastaFile: C{'file'} or C{'string'}
    @rtype : C{list} of C{class 'Bio.SeqRecord'}s
    """
    # Check that fh is a filehandler and not just a filename, the write method
    # of SeqIO works perfectly well with filehanlers but not with filenames.
    if not isinstance(refFastaFile, (file, str,)) :
      raise StandardError, 'Reference fasta file should be either a filehandler or a string of the fasta file name.'
    # Set the random number generator (add the 1000end multiplicant of the
    # rndSeedCorrector to avoid seed overlaping).
    rng = random.Random(self.rndSeed + 1000 * self.rndSeedCorrector)
    # Initialise the list of library clones.
    libraryClones = []
    # Iterate over the parser object.
    for sequence in SeqIO.parse(refFastaFile, 'fasta') :
      for c in xrange(self.libraryCoverage) :
        previousSlicePoint = 0
        while previousSlicePoint < len(sequence) :
          # While loop to generate only positive length of fragment sizes.
          while True :
            fragmentSize = int(round(rng.gauss(self.libraryInsertSize, self.libraryInsertSD)))
            if fragmentSize > 0 :
              break;
          nextSlicePoint = previousSlicePoint + fragmentSize
          fragment = sequence[previousSlicePoint:nextSlicePoint]
          fragment.description = ''
          fragmntEnd = str(nextSlicePoint)
          if nextSlicePoint >= len(sequence) :
            fragmntEnd = str(len(sequence))
          fragment.id = fragment.id + '_fragmnt_' + str(previousSlicePoint) + '_' + fragmntEnd + '_fragmntLen_' + str(len(fragment.seq))
          libraryClones.append(fragment)
          previousSlicePoint = nextSlicePoint
    return libraryClones



class NGSFactory(object) :
  """Class to represent and specify the run of a Next Generation Sequencing.

  @ivar readLength: The length of the NGS reads. In case of Illumina is the
  exact lengthin case of the 454 is the mode read length.
  @type readLength: C{'int'}
  @ivar pairedEnd: Indicates whether we get paired end sequences from the
  sequencing simulator.
  @type pairedEnd: C{'bool'}
  @ivar errorModel: Indicates whether the appropriate error model (either a
  model representing the Illumina errors or a model representing the 454 errors)
  will be incorporated.
  @type errorModel: C{'bool'}
  @ivar rng: The random number generator object.
  @type rng: C{class 'random.Random'}
  """

  def __init__(self, ngsParameters) :
    """Constructor.

    @param ngsParameters: An C{'NGSParameters'} instance specifying the
    sequencing parameters.
    @type ngsParameters: C{class 'NGSParameters'}
    """
    self.platform   = ngsParameters.sequencingPlatform
    self.readLength = ngsParameters.readLength
    self.pairedEnd  = ngsParameters.PE
    self.errorModel = ngsParameters.errorModel
    self.rng        = random.Random(ngsParameters.rndSeed)
#    self.cloneList  = libraryClones


  def __call__(self, libraryClones) :
    """Generate next generation sequencing reads.

    """
    if self.platform == 'Illumina' :
      seqReads = self.generate_illumina_reads(libraryClones)
    elif self.platform == '454' :
      seqReads = self.generate_454_reads(libraryClones)
    else :
      raise StandardError, 'Sequencing platform "%s" not identified.' % self.platform
    return seqReads


  def generate_illumina_reads(self, cloneList) :
    """Generate simulated Illumina reads.

    @rtype: C{list} of C{class 'Bio.SeqRecord'}
    """
    readsList = []
    for clone in cloneList :
      read1 = clone[:self.readLength]
      rn = self.rng.random()
      if rn >= 0.5 :
        read1.seq = read1.seq.reverse_complement()
      read1.description = ''
      read1.id = clone.id + '_read_1'
      if self.pairedEnd :
        read2 = clone[(len(clone) - self.readLength):]
        if rn < 0.5 :
          read2.seq = read2.seq.reverse_complement()
        read2.description = ''
        read2.id = clone.id + '_read_2'
        read1.id = clone.id + '_read_1'
        readsList.append(read2)
      readsList.append(read1)
    if self.errorModel :
      self.introduce_seq_errors(readsList)
    return readsList


  def generate_454_reads(self, cloneList) :
    """Generate simulated Roche 454 reads.

    """
    pass


  def introduce_seq_errors(self, readsList) :
    """Introduce sequencing errors to reads generated by the NGS simulator.

    Currently implements only an Illumina sequencing error model.
    @param readsList: A list of NGS reads.
    @type readsList: C{'list'} of C{class 'Bio.SeqRecord'}
    @rtype: None
    """

    # A block of nested functions follows.
    def generate_error_distribution(mn, mx, length) :
      """Return the error probability distribution.

      This is the actual implementation of the Dohm 2008 NAR paper error model.
      The distribution of error follows a 4th polynomial relation relative to
      the position of the nucleotide in the read.
      """
      cubesDist = []
      noCubes   = 0
      for i in xrange(0, length) :
        noCubes = noCubes + i**4
        cubesDist.append(noCubes)
      diff = mx - mn
      incr = diff / float(cubesDist[-1])
      probabilityDist = []
      for cubes in cubesDist :
        probabilityDist.append(mn + cubes * incr)
      return probabilityDist

    def substitute_nucleotide(nuc) :
      """Return another nucleotide, or N.

      The implementation at the moment returns any of the rest three
      nucleotides with an equal probability. However a more accurate model
      should take into account a nucleotide substitution matrix, derived either
      from empirical data, or from Illumina/454 biases studies for
      substitutions.
      """
      nuc = string.upper(nuc)
      nucleotides = ["N", "A", "G", "T", "C"]
      nucleotides.remove(nuc)
      rn = self.rng.random()
      if rn <=0.001 :
        return nucleotides[0]
      elif rn > 0.001 and rn <= 0.334 :
        return nucleotides[1]
      elif rn > 0.334 and rn <= 0.667 :
        return nucleotides[2]
      elif rn > 0.667 and rn <= 1 :
        return nucleotides[3]

    def impose_illumina_errors(readList) :
      """Implement an Illumina error model based on empirical data.

      The implementation and the empirical error rates are based on the
      Dohm2008 NAR paper.
      """
      # Empiricaly derived error rates.
      totalErrors     = 0
      totalErrorReads = 0
      errorMin        = 0.001
      errorMax        = 0.01
      error1          = 2.0 * errorMin # The error rate of the first nucleotide
                                       #is double the minimum one.
      # Generate the error rate distribution.
      errorModelDist = generate_error_distribution(errorMin, errorMax, self.readLength - 1)
      # Prepend the error rate for the first position.
      errorModelDist.insert(0, error1)
#      avgError = sum(errorModelDist) / float(len(errorModelDist))
      for read in readList :
        readSeq = read.seq
        # Convert the sequence to mutable so that you can alter it.
        readSeq = readSeq.tomutable()
        # The nucleotide error rates follow a distribution described by a 4th
        # order polynomial. (for details see the implememntation of the
        # generate_error_distribution method)
        errorRead = 0
        for i in xrange(len(readSeq)) :
          rn = self.rng.random()
          if rn <= errorModelDist[i] :
            errorRead = 1
            readSeq[i] = substitute_nucleotide(read[i])
            totalErrors = totalErrors + 1
        if errorRead :
          totalErrorReads = totalErrorReads + 1
        # Substitute the error imposed sequence and convert it back to
        # imutable.
        read.seq = readSeq.toseq()
        # Some sanity check prints.
#        print 'mutated : %s' % readSeq
#      print errorModelDist
#      print avgError
#      print 'total reads      : %i' % len(readsList)
#      print 'total error reads: %i' % totalErrorReads
#      print 'total errors     : %i' % totalErrors

    def impose_454_errors(readList) :
      """Implement a 454 sequencing error model.

      Not implemented so far.
      """
      return readList

    # Here is the actual implementation of the function introduce_seq_errors.
    if self.platform == 'Illumina' :
      impose_illumina_errors(readsList)
    elif self.platform == '454' :
      impose_454_errors(readsList)
    else :
      raise StandardError, 'Sequencing platform "%s" not identified.' % self.platform
    # --- End of nested functions ---



class SequenceHomolgyFactory(object) :
  """Class to represent the execution of a sequence homology search procedure.

  """

  def __init__(self, seqHomolParams, ngsParams) :
    """The constructor, collects the parameters.

    """
    self.hmmProfiles = seqHomolParams.HMMProfiles
    self.hmmEvalues  = seqHomolParams.HMMEvalues
    self.blastDb     = seqHomolParams.blastDatabase
    self.seqIdent    = seqHomolParams.seqIdentity
    self.alignLen    = seqHomolParams.seqAlignLength
    self.PEflag      = ngsParams.PE


  def __call__(self, aaFile, readsFile) :
    """The caller actualy conducts the search and returns the results.

    @param aaFile: A fasta file containing amino acid sequences.
    @type aaFile: C{'file'}
    """
    hmmReadIDs = self.conduct_hmmer_search(aaFile)
    blastReadIDs = self.conduct_BLAST_search(readsFile)
    # Combine the hmm search and the blast search results.
    seqHomolReadIDs = hmmReadIDs + blastReadIDs
    # Filter for uniq FastaIDs.
    uniqFastaIDs = list(set(seqHomolReadIDs))
    # Retrive the reads that passed the seq homology step before.
    seqHomolRecs = get_fasta_seqIDs(readsFile, uniqFastaIDs)
    return seqHomolRecs


  def conduct_hmmer_search(self, aaFile) :
    """Conduct the Hmmersearch.

    """
    fastaIDs = []
    for j in xrange(len(self.hmmProfiles)) :
      # Compile the command line call.
      hmmcmdLine = 'hmmsearch --domE ' + self.hmmEvalues[j] + ' ' + self.hmmProfiles[j] + ' ' + aaFile
      # Run the HMM and collect the read names.
      hmmSearch = subprocess.Popen(shlex.split(hmmcmdLine), bufsize = -1, stdout=subprocess.PIPE).communicate()[0]
      fastaIDs = fastaIDs + self.parse_hmmsearch_output(hmmSearch, self.PEflag)
    # Keep only the unique IDs
    uniqFastaIDs = list(set(fastaIDs))
    return uniqFastaIDs


  def conduct_BLAST_search(self, readsFile) :
    """Comnduct the BLAST search

    """
    return []


  def parse_hmmsearch_output(self, hmmSearch, paired) :
    """Extremly basic parser of the output of the hmmsearch program.

    Return only a list with the unique fasta IDs of the hits, as well as their
    pair end read incase the paired flag is true (by default).
    @param hmmSearch: The output of the hmmsearch program.
    @type hmmSearch: C{'str'}
    @param paired: Flag to specify if we want to return the pair of the read.
    @type paired: C{'bool'}
    @rtype : C{'set'}
    """
    readIDs = []
    # Get the lines that start with >> (hmmsearch output)
    for line in hmmSearch.split("\n") :
      if re.match('>>', line) :
        # Keep the ID except the last two characters (the frame name)
        readID = str(line.split()[1])[:-2]
        readIDs.append(readID)
    # Return two different list depending on paired ends.
    if paired :
      # Remove the paired end ID
      pReadIDs = [read[:-2] for read in readIDs]
      pReadIDs = list(set(pReadIDs))
      readIDs = []
      for readID in pReadIDs :
        readIDs.append(readID + '_1')
        readIDs.append(readID + '_2')
      return readIDs
    else :
      return list(set(readIDs))



class AssemblyFactory(object) :
  """Class to take care of the deNovo assemblies, assess the quality and
  compare with the reference assembly.

  """

  def __init__(self, assParams, libParams) :
    """The constructor.

    """
    self.assembler     = assParams.assemblyProgram
    self.kmers         = assParams.assemblyKmerSizeList
    self.referenceFile = assParams.assemblyReferenceFile
    self.insSize       = libParams.insertSize
    self.libSD         = libParams.standardDeviation


  def __call__(self, readsFastaFile) :
    """Caller It actually prerforms the denovo assembly and selects the one
    with the highest N50.

    """
    N50 = 0
    returnFile = ''
    bestDir    = ''
    if not self.assembler == 'velvet' :
      raise StandardError, 'Only the velvet assembler is supported at the moment.'
    else :
      for kmer in self.kmers :
        velvethCmd = 'velveth velvetAssembly_k' + str(kmer) + ' ' + str(kmer) + '-shortPaired -fasta ' + str(readsFastaFile)
        velvetgCmd = 'velvetg velvetAssembly_k' + str(kmer) + ' -ins_length_sd ' + str(self.libSD) +  ' -ins_length ' + str(self.insSize)
        # Execute the assemblies.
        subprocess.call(shlex.split(velvethCmd), bufsize = -1)
        subprocess.call(shlex.split(velvetgCmd), bufsize = -1)
        # Evaluate the assembly quality.
        assDir = 'velvetAssembly_k' + str(kmer)
        contigsFile = assDir + '/contigs.fa'
        contigN50 = FastxMetrics.calculate_N50(contigsFile)
        if contigN50 > N50 :
          returnFile = contigsFile
          bestDir = assDir
          N50 = contigN50
          print N50
          print assDir
    # Remove the extraenous assembly direcotries.
    assDirs = ['velvetAssembly_k' + str(k) for k in self.kmers]
    assDirs.remove(bestDir)
    for dr in assDirs :
      shutil.rmtree(dr)
    return returnFile


  def assembly_quality(self, contigsFile) :
    """Method to perform quality control of the assembly.

    At the moment only the N50 of the assembly is taking into account.
    """
    pass


  def assembly_comparison(self, contigsFile, referenceFile) :
    """Method to perform a comparison of a given assembly with a refference
    one.

    """
    pass



class TranslationFactory(object) :
  """Class to represent a factory to conduct a 6 frame translation.

  """

  def __init__(self, ngsSeqFile) :
    """Constructor.

    Populate the instance variable ngsParser with the appropriate parser.
    @param ngsSeqFile: A filename or a filehandler containg NGS reads.
    @type ngsSeqFile: C{'str'} or C{'file'}
    """
    self.ngsParser = SeqIO.parse(ngsSeqFile, "fasta")


  def __call__(self, frames) :
    """The caller of the class.

    Implements the translation in the number of frames asked (1-6).
    Returns a list of protein sequences.
    @param frames: The number of ORFs required.
    @type frames: C{int}
    @rtype : C{list} of C{class 'Bio.SeqIO.SeqRecord'}
    """
    # check frames is less than 6 and greater than one.
    if (frames < 1) or (6 < frames) :
      raise StandardError, 'The number of frames should be between 1 to 6.'
    translReads = []
    for read in self.ngsParser :
#      transRead = copy.deepcopy(read)
      readSeq = read.seq
      for i in xrange(frames) :
        if i < 3 :
          transRead = SeqRecord(readSeq)
          transRead.seq = readSeq[i:].translate()
          transRead.id = read.id #+ '_frame%i' % (i + 1)
          transRead.description ='frame%i' % (i + 1)
          translReads.append(transRead)
        if i >= 3 :
          transRead = SeqRecord(readSeq)
          transRead.seq = readSeq[(i - 3):].reverse_complement().translate()
          transRead.id = read.id #+ '_frame%i' % (2 - i)
          transRead.description = 'frame%i' % (2 - i)
          translReads.append(transRead)
    return translReads



# Some usefull functions

def print_SeqRecord_list(seqRecList, fh) :
  """Print the a list of Seqrecord in a file, specified by fh, in FASTA format.

  """
  if not isinstance(fh, file) :
    raise StandardError, 'Method print library requires an open ready to write filehandler to properly write all the libray clones.'
  for seqRec in seqRecList :
    SeqIO.write(seqRec, fh, 'fasta')



def get_fasta_seqIDs(fastaFile, seqIDList) :
  """Return a list of SeqRecord objects containg the sequences of the provided
   list.

  """
  seqRecList = []
  for sq in SeqIO.parse(fastaFile, "fasta") :
    if sq.id in seqIDList :
      seqRecList.append(sq)
  return seqRecList

