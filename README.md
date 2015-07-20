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
