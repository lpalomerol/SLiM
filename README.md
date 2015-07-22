# SLiM

**S**election on **Li**nked **M**utations: a forward population genetic simulation for studying linkage effects, such as hitchhiking, background selection, and Hill-Robertson interference.

SLiM can incorporate complex scenarios of demography and population substructure, various models for selection and dominance of new mutations, realistic gene and chromosome structure, and user-defined recombination maps. Emphasis was further placed on the ability to model and track individual selective sweeps – both complete and partial. While retaining all capabilities of a forward simulation, SLiM utilizes sophisticated algorithms and optimized data structures that enable simulations on the scale of entire eukaryotic chromosomes in reasonably large populations. All these features are implemented in an easy-to-use C++ command line program.

SLiM is a product of the Messer Lab at Cornell University. It was developed by Philipp Messer, and this fork has been made by Luis Palomero

GitHub home page for SLiM: [https://github.com/MesserLab/SLiM](https://github.com/MesserLab/SLiM)

Messer Lab home page for SLiM: [http://messerlab.org/software/](http://messerlab.org/software/)

## Building and Running SLiM

For Mac OS X users, an Xcode project is provided that can be used to build SLiM. For users on other platforms, or for those who prefer not to use Xcode, SLiM can be built in Terminal with the commands:

```
cd SLiM
g++ -O3 -Iinclude/ slim.cpp src/*.cpp -lgsl -lgslcblas -o slim
```

Note that SLiM uses C++11 extensions, and thus that standard is specified at compilation in order to suppress warnings.

If your GNU Standard Library headers are not in the default search paths for g++, you will need to supply them on the command line.  You can find out the right command-line arguments to use for this by executing:

```
gsl-config --cflags --libs
```

For example, I have installed gsl using MacPorts, so my compilation command looks like:

```
g++ -O3 ./core/*.cpp -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas -std=c++11 -o slim
```

Once SLiM is built, just run it at Terminal's command line. For example, to run the first example provided in SLiM's distribution, execute:

```
cd SLiM
./slim ./examples/input_example_1.txt
```

If you have made a Release build of SLiM with Xcode, it should be at /usr/local/bin/slim; you can provide that path or ensure that it is part of your shell's default search path.

### Changes of this version

This version adds features related with domestication processes. These functionalities are:
 - Manage different fitnesses by population, simulating different environments. 
 - Use sexed populations also with hermaphrodite ones.
 - Choose a fixed number of reproductive males from all the population males.
 - Choose a ratio of females which will be reproduced at every generation. 

#### Different fitness by population

Original SLiM version assumes that the different populations are at similar environments. This version has been thought for taking account different groups of environments (natural, domestic, etc) is necessary simulate different fitness per environment. Functionally has been made necessary update the `#GENOMIC ELEMENT TYPES` parameter including the population id as optional parameter.

```
#GENOMIC ELEMENT TYPES
p1 g1 m1 0.70 m2 0.25 m3 0.05 / exon (70% deleterious, 25% neutral, 5% positive) for population1
p2 g1 m2 0.50 m4 0.50 / exon (50% neutral, 50% deleterious fixed) for population 2
g1 m1 0.15 m2 0.85 / exon (15% neutral,85% deleterious fixed). Default value for population 3
g2 m1 0.50 m2 0.50 / UTR (50% deleterious, 50% neutral) for all p opulations
g3 m2 1.0 / intron (100% neutral) for all populations
```

#### Include sexed populations 

SLiM simulates hermaphrodite populations and the new requirement requires simulate sexed populations. Functionally a new parameter `#HERMAPRHDOTIES` has been added to the list of parameters. This optional parameter has only 2 possible values: 0 and 1, for sexed populations and hermaphrodites respectively. By default, the populations are hermaphrodite.

```
#HERMAPHRODITES 
0
```

When the simulated individuals are sexed, the system generates the number of males and females randomly following a binomial distribution with p=0.5 allowing different number of males and females every generation.


#### New Demography Event: Include reproductive threshold

One of the most basic techics used at artificial selection is limit the number of males that are going to be reproduced at every generation. This number is fixed and this subset of males is selected randomly in bases of fitness following the same distribution than the used when the next generation individuals are generated. This subset of individuals does not repeat any individuals.

At below sample the two populations are defined and p2 only allows to reproduce only 100 of the males.

```
#DEMOGRAPHY AND STRUCTURE
1 P p1 500 / one population of 500 individuals
1 P p2 500
1 T p2 100 / Threshold: 100 of the individuals will cross with all others
```

#### New Demography Event: Include sex ratio of females

Continuing with ancient techniques related with artificial selection, this fork includes the sex ratio, which implies select a variable number of females related with the males. For example, with a population of 1000 individuals (400 males and 600 females) if the ratio is 2 females per male and the threshold is 200 males, the number of reproductive females will be 400. From this selection of 600 individuals (200 males and 400 females) the next generation of 1000 individuals will be generated.

```
#DEMOGRAPHY AND STRUCTURE
1 P p1 1000 / one population of 1000 individuals
1 R p1 200 / Threshold: 200 of the males will cross with all others
1 R p1 0.5 / Sex ratio: Two famales per female
```

Using the same idea than thresholds, a subset of reproductive females will be defined before generating the new population.



## Credits

Forked by Luis Palomero under guidance of Sebastián Ramos from CRAG (http://www.cragenomica.es/)

## License

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version (http://www.gnu.org/licenses/).

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
