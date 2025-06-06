# AttenuationNC

These files contain the code required to replicate observations reported in the paper *RNA-ligand complexes and the attenuation of neutral confinement in the evolution of RNA secondary structures*.

Our paper is currently undergoing review. We will include information on the authors, contact information, institutions and journal whenever the paper becomes accepted for publication.

Our code uses tools from the ViennaRNA library (Lorenz et al., 2011. ViennaRNA package 2.0, Algorithms for molecular biology, **6**:26) and from the GNU Scientific library, GSL (Galassi et al. 2009, GNU Scientific Library reference manual. Network theory ltd). 

##Content

### libs/
This directory contains the libraries that we developed for this project. Remember to link them when compiling.

### Structures/
This directory contains information on the RNA secondary structures that we analyzed for the paper.

### sample/
This directory contains the code that we used to sample sequences that fold into a pre-defined RNA secondary structure.

### properties.cc
Code for the analyses of large sets of individual sequences.

### moderna.cc
Code for simulations of the evolution of RNA secondary structures.

### illust.cc
Code for analyses of simplified hypothetical molecules of RNA.

