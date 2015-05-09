#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>
#include "translator.h"
using namespace std;

genomicElementTypesLine translate_genomic_elements_type_line(string line){
	genomicElementTypesLine line_translated;
	string sub;
    istringstream iss(line);
    //TODO: Add value for population id
    iss >> sub;
    if(sub.at(0) == 'p'){
      sub.erase(0, 1);
      line_translated.population = atoi(sub.c_str());
      iss >> sub;
    } else {
      line_translated.population = 0;
    }
    sub.erase(0, 1);
    line_translated.id = atoi(sub.c_str());
    while (iss >> sub) {
      sub.erase(0, 1);
      line_translated.mutations.push_back(atoi(sub.c_str()));
      iss >> sub;
      line_translated.fractions.push_back(atof(sub.c_str()));
    }

	return line_translated;
}
