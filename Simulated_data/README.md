## Instructions to simulate data 

The simulated data used in Oliveira et al., 2024 (Assessing the limits of local ancestry inference from small reference panels) can be generated with *fastsimcoal* v2.8 using the definition file (e.g., *SimAdm.def*) and the template file (e.g., *SimAdm.tpl*).

Here is an example code for the model with one admixture event:

`./ fsc28 -t SimAdm.tpl -f SimAdm.def -n1 -c20 -B20 -I -s0 -q -x -j -k1000000 -G -p --seed 12345`


Each line of the definition file contains a set of parameter values to be used in a particular simulation, whose model is defined by the template file.
The files created (*.gen) are numbered from 1 to the number of parameter combinations (lines).

The parameters used in the example file *SimAdm.def* indicate:


	SS1 – Sample size for source 1
	SS2 – Sample size for source 2
	NS1 – Ne for source 1
	NS2 – Ne for source 2
	ST – Sample size of the target population(s)
	NT – Ne for target (admixed) population
	Nanc – Ne of ancestral population
	TDiv – Time of divergence
	TAdm – Time of admixture

For more information on how to run *fastsimcoal* v2.8 or set up the input files, please consult the software manual.
