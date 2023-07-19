# simpLAI
Local Ancestry Inference (LAI) from small reference panels

## Overview
simpLAI is a window-based approach that leverages information on the matching of haplotypes from admixed and reference individuals to infer the ancestry across the genome.

The program performs three types of inferences. In the first two modes, the ancestry assigned to a target haplotype of predefined size (*-s*) is that of the reference population that displays the lowest average number of mismatches with the target (*mean*) or that includes the chromosome displaying the lowest number of mismatches with the target (*min*). The third mode accounts for intra-population ancestry switches that may have occurred within a haplotype through recombination (*rec*). Each haplotype segment consists of n polymorphic sites (*-n*) and the number of mismatches between a target and each chromosome of the reference population is computed on subsets of t linked polymorphic sites (assumed to be a non-recombining segment). The path consisting of the combination of segments (potentially identified on different chromosomes) among the n polymorphic sites showing the smallest number of mismatches within each reference population is selected to compare reference populations. The reference with the shortest path is selected for the ancestry of the target haplotype. The increment for the sliding window is defined in base pairs (*-i*) for the mean and min mode and in number of polymorphic sites (*-m*) for the rec mode.

#### Usage:
```console
 -h  --help               : prints this usage  
 -g  --genFile test.gen   : name of genotype file (phased)  
     --ss1  10            : number of chromosome in source 1  
     --ss2  10            : number of chromosome in source 2  
     --ss3  10            : number of chromosome in source 3  
     --ss4  10            : number of chromosome in source 4  
     --ssa  20            : number of chromosome in admixed sample  
 -l  --genomeLength  1e8  : total length of genome  
 -s  --winSize       1e6  : size of sliding window (optional: default = 1Mb)  
 -i  --winInc        5e5  : increment for sliding window (optional: default = 500Kb)  
 -n  --numPolSites   1e4  : number of polymorphic sites on which to search for  
                            the shortest recombination path (default = 10,000 sites)  
 -m  --incRec        1e3  : increment for sliding window when searching for  
                            the shortest recombination path (default = 1,000 sites)  
 -t                  5    : no. of linked polymorphic sites between potential  
                            recombination breakpoints (default = 5 sites)  
```
