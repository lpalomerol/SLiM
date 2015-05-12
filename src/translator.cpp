#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>

#include "translator.h"
using namespace std;

genomic_element_types_line translate_genomic_elements_type_line(string line){
	genomic_element_types_line line_translated;
	string sub;
    istringstream iss(line);
    //TODO: Add value for population id
    iss >> sub;
    if(sub.at(0) == 'p'){
      sub.erase(0, 1);
      line_translated.pop_id = atoi(sub.c_str());
      iss >> sub;
    } else {
      line_translated.pop_id = 0;
    }
    sub.erase(0, 1);
    line_translated.genomic_id = atoi(sub.c_str());
    while (iss >> sub) {
      sub.erase(0, 1);
      line_translated.mutations.push_back(atoi(sub.c_str()));
      iss >> sub;
      line_translated.fractions.push_back(atof(sub.c_str()));
    }

	return line_translated;
}
