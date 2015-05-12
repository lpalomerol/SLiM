// translator.h

#ifndef TRANSLATOR_H_
#define TRANSLATOR_H_

#include <iostream>
#include <string>

struct genomic_element_types_line {
  int genomic_id;
  int pop_id;
  std::vector<int> mutations;
  std::vector<double> fractions;
};

genomic_element_types_line translate_genomic_elements_type_line(std::string);

#endif /* TRANSLATOR_H_ */
