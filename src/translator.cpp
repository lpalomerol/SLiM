
#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>

#include "translator.h"
using namespace std;

genomic_element_types_line translate_genomic_elements_type_line(string line) {
  genomic_element_types_line line_translated;
  string sub;
  istringstream iss(line);
  // TODO: Add value for population id
  iss >> sub;
  if (sub.at(0) == 'p') {
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

int translate_hermaphrodite_line(string line) {
  int hermaphrodite_code;
  string sub;
  istringstream iss(line);
  iss >> sub;
  hermaphrodite_code = atoi(sub.c_str());
  return hermaphrodite_code;
}

output_line translate_R_output(string line) {

  output_line translated_line;
  bool is_ms = (line.find("MS") != string::npos);
  bool has_seed = false;
  string sub;
  istringstream iss(line);
  iss >> sub;
  translated_line.time = (int)atof(sub.c_str());
  iss >> sub;
  translated_line.event = sub.at(0);

  while (iss >> sub) {

    if (sub != "MS") {
      translated_line.params.push_back(sub.c_str());
    }
  }
  // With even params, lacks seed
  if (translated_line.params.size() % 2 == 0) {
    translated_line.params.push_back("0");
  }
  translated_line.params.push_back((is_ms == true) ? "MS" : "DF");

  return translated_line;
}
