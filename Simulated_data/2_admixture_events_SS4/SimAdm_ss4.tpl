//Parameters for the coalescence simulation program : fsimcoal2.exe
5 samples to simulate :
//Population effective sizes (number of genes)
NHG
NHG
NADM
NEF
NNE
//Samples sizes and samples age
0
SSHG
SSADMEF
SSAEG
SSNE
//Growth rates	: negative growth implies population expansion
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
7  historical event
tadm1 2 1 padm1 1 0 0 recordAdmOrigin
tadm1 2 3 1 1 0 0 
tadm1 3 2 1 1 0 0
tadm1 1 0 1 1 0 0
tadm2 2 0 padm2 1 0 0 recordAdmOrigin
tadm2 2 4 1 1 0 0 recordAdmOrigin
tdiv 0 4 1  NANC 0 0 absoluteResize
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
DNA  1e8   1e-8   1.25e-8 OUTEXP

