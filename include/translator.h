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

struct output_line {
  int time;
  char event;
  std::vector<std::string> params;
};

genomic_element_types_line translate_genomic_elements_type_line(std::string);
int translate_hermaphrodite_line(std::string);
output_line translate_R_output(std::string);

#endif /* TRANSLATOR_H_ */
