//translator.h


#ifndef TRANSLATOR_H_
#define TRANSLATOR_H_

#include <iostream>
#include <string>

struct genomicElementTypesLine {
        int id;
	int population;
	std::vector<int> mutations;
	std::vector<double> fractions;
};

genomicElementTypesLine translate_genomic_elements_type_line(std::string);

#endif /* TRANSLATOR_H_ */
