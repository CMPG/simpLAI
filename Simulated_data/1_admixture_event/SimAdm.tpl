//Parameters for the coalescence simulation program : fsimcoal2.exe
7 samples to simulate :
//Population effective sizes (number of genes)
NS1
NS2
NT
NT
NT
NT
NT
//Samples sizes and samples age
SS1 0 0
SS2 0 0
ST 0 0
ST 0 0
ST 0 0
ST 0 0
ST 0 0
//Growth rates : negative growth implies population expansion
0
0
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
10  historical event
TAdm 3 0 0.05 1 0 0 recordAdmOrigin
TAdm 4 0 0.1 1 0 0 recordAdmOrigin
TAdm 5 0 0.2 1 0 0 recordAdmOrigin
TAdm 6 0 0.3 1 0 0 recordAdmOrigin
TAdm 2 1 1 1 0 0 
TAdm 3 1 1 1 0 0 recordAdmOrigin
TAdm 4 1 1 1 0 0 recordAdmOrigin
TAdm 5 1 1 1 0 0 recordAdmOrigin
TAdm 6 1 1 1 0 0 recordAdmOrigin
TDiv 0 1 1 NAnc 0 0 absoluteResize
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
DNA 1e8 1e-8 1.25e-8 OUTEXP
