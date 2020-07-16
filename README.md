Genotype-phenotype landscapes:

The folder "gml_landscapes" contains the network files in GML format for 746 genotype-phenotype landscapes of TF binding sites. More information about these landscapes can be found in Table S1. Each network file has the following vertex attributes: - id: vertex identification number. - genotype: the nucleotide sequence of the binding site. - score: the enrichment score in protein binding microarrays of the sequence.

Landscape global peaks:

The folder "peaks.zip" contains the text files with the sequences that belong to the global peak (i.e. the peak with the highest bbinding affinity) for each genotype-phenotype landscape.

Evolutionary simulations:

The file "wt-evol.c" contains the basic code used to simulate the evolution of TF binding sites in the polymorphic regime. It requires as input a GML file with landscape details. Additionally, the required initial arguments include: - output file. - number of generations. - number of replicates per initial condition. - population size. - mutation rate. - mutation bias parameter \alpha. - seed for the random number generator.

Composition bias:

The file "biascount.c" contains the code used to compute the composition bias in the entire landscape and along accesible paths to the global peak. It requires as an input a GML file with landscape details. The required initial arguments include: - output file. - noise threshold parameter \delta. - corresponding pk file with the set of sequence in the global peak. The noise threshold parameter \delta corresponding to each landscape can be found in the Supplementary Table 1.
