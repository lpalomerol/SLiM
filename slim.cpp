///
/// @mainpage slim.cpp
///  SLiM: a forward population genetic simulation with selection and linkage
///
/// version 1.8 (November 4, 2014)
///
/// Copyright (C) 2013  Philipp Messer
///
/// compile by:
///
/// g++ -O3 -Iinclude/ slim.cpp src/*.cpp -lgsl -lgslcblas -o slim
///
/// This program is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version (http://www.gnu.org/licenses/).
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <unistd.h>
#include <cstring>

#include "./include/translator.h"

using namespace std;

const gsl_rng* rng;

/**
 * @class mutation
 */

class mutation {
 public:
  int t;   /**< mutation type identifier */
  int x;   /**< position */
  float s; /**< selection coefficient */
  int i;   /**< subpopulation in which mutation arose */
  int g;   /**< generation in which mutation arose */

  mutation(void) { ; }

  mutation(int T, int X, float S, int I, int G) {
    t = T;
    x = X;
    s = S;
    i = I;
    g = G;
  }
};

bool operator<(const mutation& M1, const mutation& M2) {
  return M1.x < M2.x;
};

bool operator==(const mutation& M1, const mutation& M2) {
  return (M1.x == M2.x && M1.t == M2.t && M1.s == M2.s);
};

/** @class event
 * type of events:
 *
 * t P i n [j]:  add subpopulation i of size n [drawn from j]
 *
 * t N i n:      set size of subpopulation i to n
 *
 * t M i j x:    set fraction x of subpopulation i that originates as migrants
 *
 * from j
 *
 * t S i s;      set selfing fraction of subpopulation i to s
 *
 * t R i n:      output sample of n randomly drawn genomes from subpopulation
 *
 * i
 *
 * t F:          output list of all mutations that have become fixed so far
 *
 * t A [file]:   output state of entire population [into file]
 *
 * t T m:        follow trajectory of mutation m (specified by mutation type)
 *
 * from generation t on
 */
class event {

 public:
  char t;           /**< event type */
  vector<string> s; /**< vector of strings with parameters of event */
  int np;           /**< number of parameters */

  event(char type, vector<string> event_names) {
    t = type;
    s = event_names;
    np = s.size();

    string options = "PNMStrRFAT";
    if (options.find(t) == string::npos) {
      cerr << "ERROR (initialize): invalid event type \"" << t;
      for (int i = 0; i < np; i++) {
        cerr << " " << s[i];
      }
      cerr << "\"" << endl;
      exit(1);
    }
  }
};

class mutation_type {
  /** @class mutation_type
   * a mutation type is specified by the DFE and the dominance coefficient
   *
   * DFE options: f: fixed (s)
   *              e: exponential (mean s)
   *              g: gamma distribution (mean s,shape)
   *
   * examples: synonymous, nonsynonymous, adaptive, etc.
   */
 public:
  float h;          /**< dominance coefficient */
  char d;           /**< DFE (f: fixed, g: gamma, e: exponential) */
  vector<double> p; /**< DFE parameters */

  /**
   * @brief constructor
   * @param float H Dominance coefficient
   * @param char D DFE (f: fixed, g: gamma, e: exponential)
   * @param vector<double> P extra parameters
   */
  mutation_type(float H, char D, vector<double> P) {
    h = H;
    d = D;
    p = P;

    string s = "fge";

    if (s.find(d) == string::npos) {
      cerr << "ERROR (initialize): invalid mutation type parameters" << endl;
      exit(1);
    }
    if (p.size() == 0) {
      cerr << "ERROR (initialize): invalid mutation type parameters" << endl;
      exit(1);
    }
  }
  /**
   * @brief Returns extra parameter value of this mutaiton type
   * @return float
   */
  float draw_s() {
    switch (d) {
      case 'f':
        return p[0];
      case 'g':
        return gsl_ran_gamma(rng, p[1], p[0] / p[1]);
      case 'e':
        return gsl_ran_exponential(rng, p[0]);
      default:
        exit(1);
    }
  }
};

/** @class genomic_element
 * a genomic element has a genomic element type identifier (i), start (s) and
 * end (e) position
 */
class genomic_element {

 public:
  int i; /**< identifier */
  int s; /**< start */
  int e; /**< end*/

  /**
   * @brief builder
   * @param integer I identifier
   * @param integer S start
   * @param integer E End
   */
  genomic_element(int I, int S, int E) {
    i = I;
    s = S;
    e = E;
  }
};

class genomic_element_type {
  /** @class genomic_element_tpe
   * a genomic element type is specified by a vector of the mutation type
   * identifiers off all
   * mutation types than can occur in such elements and a vector of their
   * relative fractions.
   * examples: exon, intron, utr, intergenic, etc.
   */
 private:
  gsl_ran_discrete_t* LT;

 public:
  vector<int> m;    /**< mutation types identifiers in this element */
  vector<double> g; /**< relative fractions of each mutation type */

  genomic_element_type(vector<int> M, vector<double> G) {
    m = M;
    g = G;

    if (m.size() != g.size()) {
      exit(1);
    }
    double A[m.size()];
    for (int i = 0; i < m.size(); i++) {
      A[i] = g[i];
    }
    LT = gsl_ran_discrete_preproc(G.size(), A);
  }

  int draw_mutation_type() { return m[gsl_ran_discrete(rng, LT)]; }
};

class chromosome : public vector<genomic_element> {
  /** @class chromosome
   * the chromosome is a vector of genomic elements (type, start, end)
   */
 private:
  gsl_ran_discrete_t* LT_M;  // mutation
  gsl_ran_discrete_t* LT_R;  // recombination

 public:
  map<int, mutation_type> mutation_types;
  map<pair<int, int>, genomic_element_type> genomic_element_types;
  vector<int> rec_x;
  vector<double> rec_r;

  int L;      /**< length of chromosome */
  double M;   /**< overall mutation rate */
  double R;   /**< overall recombination rate */
  double G_f; /**< gene conversion fraction */
  double G_l; /**< average stretch length */

  chromosome() {
    L = 0;
    M = 0.0;
    R = 0.0;
    G_f = 0.0;
    G_l = 0.0;
  }

  /**
   * @brief Validates the incoming chromosome
   * @param subpopulations List whith all subpopulations related with the
   * experiment.
   *                       Are defined at "demography and structure" section
   */
  void validate(vector<int> subpopulations) {
    if (size() == 0) {
      cerr << "ERROR (initialize): empty chromosome" << endl;
      exit(1);
    }
    if (rec_r.size() == 0) {
      cerr << "ERROR (initialize): recombination rate not specified" << endl;
      exit(1);
    }
    if (!(M >= 0)) {
      cerr << "ERROR (initialize): invalid mutation rate" << endl;
      exit(1);
    }
    for (int i = 1; i <= subpopulations.size(); i++) {
      for (int j = 0; j < size(); j++) {
        // Here validates that at least exist one default genomic element type.
        pair<int, int> row_default = make_pair(0, operator[](j).i);
        pair<int, int> row_custom_pop = make_pair(i, operator[](j).i);
        if (genomic_element_types.count(row_default) == 0 &&
            genomic_element_types.count(row_custom_pop) == 0) {
          cerr << "ERROR (initialize): genomic element type (" << i << ","
               << operator[](j).i
               << ") not defined, it must have at least a default." << endl;
          exit(1);
        }
      }
    }
    // TODO: Add tests
    // Validates that all the mutation types are well defined.
    for (map<pair<int, int>, genomic_element_type>::iterator it =
             genomic_element_types.begin();
         it != genomic_element_types.end(); it++) {
      for (int j = 0; j < it->second.m.size(); j++) {
        if (mutation_types.count(it->second.m[j]) == 0) {
          cerr << "ERROR (initialize): mutation type " << it->second.m[j]
               << " not defined" << endl;
          exit(1);
        }
      }
    }
  }

  void initialize_rng() {

    L = 0;

    double A[size()];
    int l = 0;
    for (int i = 0; i < size(); i++) {
      if (operator[](i).e > L) {
        L = operator[](i).e;
      }
      int l_i = operator[](i).e - operator[](i).s + 1.0;
      A[i] = (double)l_i;
      l += l_i;
    }

    LT_M = gsl_ran_discrete_preproc(size(), A);
    M = M * (double)l;

    double B[rec_r.size()];
    B[0] = rec_r[0] * (double)rec_x[0];
    R += B[0];
    for (int i = 1; i < rec_r.size(); i++) {
      B[i] = rec_r[i] * (double)(rec_x[i] - rec_x[i - 1]);
      R += B[i];
      if (rec_x[i] > L) {
        L = rec_x[i];
      }
    }
    LT_R = gsl_ran_discrete_preproc(rec_r.size(), B);
  }

  int draw_n_mut() { return gsl_ran_poisson(rng, M); }

  mutation draw_new_mut(int pop_id, int g) {
    int ge = gsl_ran_discrete(rng, LT_M);  // genomic element
    pair<int, int> key;
    if (genomic_element_types.count(make_pair(pop_id, operator[](ge).i)) > 0) {
      key = (std::make_pair(pop_id, operator[](ge).i));
    } else {
      key = (std::make_pair(0, operator[](ge).i));
    }
    genomic_element_type ge_type =
        genomic_element_types.find(key)->second;  // genomic element type

    int mut_type_id = ge_type.draw_mutation_type();  // mutation type id
    mutation_type mut_type =
        mutation_types.find(mut_type_id)->second;  // mutation type

    int x = operator[](ge).s +
            gsl_rng_uniform_int(
                rng, operator[](ge).e - operator[](ge).s + 1);  // position
    float s = mut_type.draw_s();  // selection coefficient

    return mutation(mut_type_id, x, s, pop_id, g);
  }

  vector<int> draw_breakpoints() {
    vector<int> r;

    // draw recombination breakpoints

    int nr = gsl_ran_poisson(rng, R);

    for (int i = 0; i < nr; i++) {
      int x = 0;
      int j = gsl_ran_discrete(rng, LT_R);

      if (j == 0) {
        x = gsl_rng_uniform_int(rng, rec_x[j]);
      } else {
        x = rec_x[j - 1] + gsl_rng_uniform_int(rng, rec_x[j] - rec_x[j - 1]);
      }

      r.push_back(x);

      if (gsl_rng_uniform(rng) <
          G_f)  // recombination results in gene conversion
      {
        int x2 = x + gsl_ran_geometric(rng, 1.0 / G_l);
        r.push_back(x2);
      }
    }

    return r;
  }
};

class polymorphism {
  /**
   * @class polymorphism
   */
 public:
  int id;  /**< mutation id */
  int t;   /**< mutation type */
  float s; /**< selection coefficient */
  int i;   /**< subpopulation in which mutation arose */
  int g;   /**< generation in which mutation arose  */
  int n;   /**< prevalence */

  polymorphism(int ID, int T, float S, int I, int G, int N) {
    id = ID;
    t = T;
    s = S;
    i = I;
    g = G;
    n = N;
  }

  void print(int x, chromosome& chr) {
    float h = chr.mutation_types.find(t)->second.h;
    cout << id << " m" << t << " " << x + 1 << " " << s << " " << h << " p" << i
         << " " << g << " " << n << endl;
  }

  void print(ofstream& outfile, int x, chromosome& chr) {
    float h = chr.mutation_types.find(t)->second.h;
    outfile << id << " m" << t << " " << x + 1 << " " << s << " " << h << " p"
            << i << " " << g << " " << n << endl;
  }

  void print_no_id(int x, chromosome& chr) {
    float h = chr.mutation_types.find(t)->second.h;
    cout << "m" << t << " " << x + 1 << " " << s << " " << h << " p" << i << " "
         << g << " " << n << endl;
  }
};

class substitution {
  /**
   * @class substitution
   */
 public:
  int t;   /**< mutation type */
  int x;   /**< position */
  float s; /**< selection coefficient */
  int i;   /**< subpopulation in which mutation arose */
  int g;   /**< generation in which mutation arose */
  int f;   /**< fixation time */

  substitution(mutation M, int F) {
    t = M.t;
    x = M.x;
    s = M.s;
    i = M.i;
    g = M.g;
    f = F;
  }

  void print(chromosome& chr) {
    float h = chr.mutation_types.find(t)->second.h;
    cout << " m" << t << " " << x + 1 << " " << s << " " << h << " p" << i
         << " " << g << " " << f << endl;
  }
};

/**
 * @class introduced_mutation
 */
class introduced_mutation : public mutation {

 public:
  int i;   /**< subpopulation into which mutation is introduced */
  int g;   /**< generation in which mutation is introduced */
  int nAA; /**< number of homozygotes */
  int nAa; /**< number of heterozygotes */

  introduced_mutation(int T, int X, int I, int G, int NAA, int NAa) {
    t = T;
    x = X;
    i = I;
    g = G;
    nAA = NAA;
    nAa = NAa;
  }
};

/**
 * @class partial_sweep
 */

class partial_sweep {
 public:
  int t;
  int x;
  float p;

  partial_sweep(int T, int X, float P) {
    t = T;
    x = X;
    p = P;
  }
};

class genome : public vector<mutation> {};

/**
 * @brief return genome G consisting only of the mutations that are present in
 * both G1 and G2
 * @param G1 First genome to intersect
 * @param G2 Second genome to intersect
 * @return genome Intersected genome
 */
genome fixed(genome& G1, genome& G2) {

  genome G;

  vector<mutation>::iterator g1 = G1.begin();
  vector<mutation>::iterator g2 = G2.begin();

  vector<mutation>::iterator g1_max = G1.end();
  vector<mutation>::iterator g2_max = G2.end();

  while (g1 != g1_max && g2 != g2_max) {
    // advance g1 while g1.x < g2.x

    while (g1 != g1_max && g2 != g2_max && (*g1).x < (*g2).x) {
      g1++;
    }

    // advance g2 while g1.x < g2.x

    while (g1 != g1_max && g2 != g2_max && (*g2).x < (*g1).x) {
      g2++;
    }

    // identify shared mutations at positions x and add to G

    if (g2 != g2_max && g1 != g1_max && (*g2).x == (*g1).x) {
      int x = (*g1).x;

      vector<mutation>::iterator temp;

      while (g1 != g1_max && (*g1).x == x) {
        temp = g2;
        while (temp != g2_max && (*temp).x == x) {
          if ((*temp).t == (*g1).t && (*temp).s == (*g1).s) {
            G.push_back(*g1);
          }
          temp++;
        }
        g1++;
      }
    }
  }

  return G;
}

/**
 * @brief return genome G consisting only of the mutations in G1 that are not in
 * G2
 * @param G1 genome source 1 (minuendo)
 * @param G2 genome to check (sustraendo)
 * @return genome The substracted genome
 */
genome polymorphic(genome& G1, genome& G2) {
  //

  genome G;

  vector<mutation>::iterator g1 = G1.begin();
  vector<mutation>::iterator g2 = G2.begin();

  vector<mutation>::iterator g1_max = G1.end();
  vector<mutation>::iterator g2_max = G2.end();

  while (g1 != g1_max && g2 != g2_max) {
    // advance g1 while g1.x < g2.x

    while (g1 != g1_max && g2 != g2_max && (*g1).x < (*g2).x) {
      G.push_back(*g1);
      g1++;
    }

    // advance g2 while g1.x < g2.x

    while (g2 != g2_max && g1 != g1_max && (*g2).x < (*g1).x) {
      g2++;
    }

    // identify polymorphic mutations at positions x and add to G

    if (g2 != g2_max && g1 != g1_max && (*g2).x == (*g1).x) {
      int x = (*g1).x;

      // go through g1 and check for those mutations that are not present in g2

      vector<mutation>::iterator temp = g2;

      while (g1 != g1_max && (*g1).x == x) {
        bool poly = 1;

        while (temp != g2_max && (*temp).x == x) {
          if ((*g1).t == (*temp).t && (*g1).s == (*temp).s) {
            poly = 0;
          }
          temp++;
        }
        if (poly == 1) {
          G.push_back(*g1);
        }
        g1++;
      }

      while (g2 != g2_max && (*g2).x == x) {
        g2++;
      }
    }
  }

  while (g1 != g1_max) {
    G.push_back(*g1);
    g1++;
  }

  return G;
}

class subpopulation {
 protected:
  gsl_ran_discrete_t* LT;
  gsl_ran_discrete_t* LT_TH;
  gsl_ran_discrete_t* LT_males;

 public:
  int N;          /**< population size */
  double S;       /**< selfing fraction */
  double S_ratio; /**< sex_ratio fraction */
  int T;          /**< threshold ratio */

  vector<genome> G_parent;
  vector<genome> G_child;

  map<int, double> m;  // m[i]: fraction made up of migrants from subpopulation
                       // i per generation

  virtual void set_selfing(double s) = 0;

  int draw_individual() { return gsl_ran_discrete(rng, LT); }

  int draw_individual_threshold() { return gsl_ran_discrete(rng, LT_TH); }
  /**
   * @brief calculate fitnesses in parent population and create new lookup table
   * @param chr Chromosome population
   * @return void
   */
  virtual void update_fitness(chromosome& chr) = 0;

  /**
   * @brief updates the fitness of selected males (or acting as male) for
   * threshold purposes
   * @param double* List of current fitness for all individuals
   * @param int number of individuals
   */
  void update_threshold_fitness(double* All, int parent_size) {
    if (T > 0) {
      double Males[parent_size];
      int n_males = T;
      memset(Males, 0, sizeof(Males));
      int j = 0;
      int male;
      while (j < n_males) {
        male = draw_individual_by_sex(true);
        if (Males[male] == 0) {
          Males[male] = All[male];
          j++;
        }
      }
      LT_males = gsl_ran_discrete_preproc(parent_size, Males);
    } else {
      LT_males = LT;
    }
  }

  /**
   * @brief calculate the fitness of the individual constituted by genomes i and
   * j in the parent population
   * @param i First genome id to compute
   * @param j Second genome id to compute
   * @param chr chromosome to check
   * @return double The computed fitness
   */
  double W(int i, int j, chromosome& chr) {

    double w = 1.0;

    vector<mutation>::iterator pi = G_parent[i].begin();
    vector<mutation>::iterator pj = G_parent[j].begin();

    vector<mutation>::iterator pi_max = G_parent[i].end();
    vector<mutation>::iterator pj_max = G_parent[j].end();

    while (w > 0 && (pi != pi_max || pj != pj_max)) {
      // advance i while pi.x < pj.x

      while (pi != pi_max && (pj == pj_max || (*pi).x < (*pj).x) && w > 0) {
        if ((*pi).s != 0) {
          w = w * (1.0 + chr.mutation_types.find((*pi).t)->second.h * (*pi).s);
        }
        pi++;
      }

      // advance j while pj.x < pi.x

      while (pj != pj_max && (pi == pi_max || (*pj).x < (*pi).x) && w > 0) {
        if ((*pj).s != 0) {
          w = w * (1.0 + chr.mutation_types.find((*pj).t)->second.h * (*pj).s);
        }
        pj++;
      }

      // check for homozygotes and heterozygotes at x

      if (pi != pi_max && pj != pj_max && (*pj).x == (*pi).x) {
        int x = (*pi).x;

        vector<mutation>::iterator pi_start = pi;

        // advance through pi

        while (pi != pi_max && (*pi).x == x && w > 0) {
          if ((*pi).s != 0.0) {
            vector<mutation>::iterator temp_j = pj;
            bool homo = 0;

            while (homo == 0 && temp_j != pj_max && (*temp_j).x == x && w > 0) {
              if ((*pi).t == (*temp_j).t && (*pi).s == (*temp_j).s) {
                w = w * (1.0 + (*pi).s);
                homo = 1;
              }
              temp_j++;
            }

            if (homo == 0) {
              w = w *
                  (1.0 + chr.mutation_types.find((*pi).t)->second.h * (*pi).s);
            }
          }
          pi++;
        }

        // advance through pj

        while (pj != pj_max && (*pj).x == x && w > 0) {
          if ((*pj).s != 0.0) {
            vector<mutation>::iterator temp_i = pi_start;
            bool homo = 0;

            while (homo == 0 && temp_i != pi_max && (*temp_i).x == x && w > 0) {
              if ((*pj).t == (*temp_i).t && (*pj).s == (*temp_i).s) {
                homo = 1;
              }
              temp_i++;
            }
            if (homo == 0) {
              w = w *
                  (1.0 + chr.mutation_types.find((*pj).t)->second.h * (*pj).s);
            }
          }
          pj++;
        }
      }
    }

    if (w < 0) {
      w = 0.0;
    }

    return w;
  }
  virtual void swap() = 0;
  subpopulation() {};
  ~subpopulation() {};

  virtual int draw_individual_by_sex(bool male) = 0;
  virtual void child_is_male(int i, bool male) {}
  virtual void parent_is_male(int i, bool male) {};
  virtual bool parent_is_male(int i) = 0;
};

class subpopulation_sexed : public subpopulation {
  /**
   * @class subpopulation_sexed
   *  a subpopulation is described by the vector G of 2N genomes
   *  individual i is constituted by the two genomes 2*i and 2*i+1
   */

 private:
  gsl_ran_discrete_t* LT_females;

 public:
  vector<bool> males_parent;
  vector<bool> males_child;

  static int get_reproductive_females(int N, int males, int threshold,
                                      double ratio) {
    if (ratio < 0 || ratio > 1.0) {
      cerr << "ERROR (get_reproductive_females): Invalid sex ratio value "
           << ratio << endl;
      exit(1);
    }

    int all_females = N - males;
    if (ratio == 0) {
      return all_females;
    }

    int threshold_males = (threshold == 0) ? males : min(threshold, males);
    double ratio_females = 1 / ratio;
    return min(int(threshold_males * ratio_females), all_females);
  }

  subpopulation_sexed(int n) {
    N = n;
    S = 0.0;
    T = 0;

    S_ratio = 0;

    G_parent.resize(2 * N);
    G_child.resize(2 * N);

    males_parent.resize(N);
    males_child.resize(N);

    double All[N];
    double Males[N];
    double Females[N];
    for (int i = 0; i < N; i++) {
      All[i] = 1.0;
      males_parent[i] = (gsl_ran_binomial(rng, 0.5, 1) == 1);
      if (males_parent[i] == true) {
        Males[i] = 1.0;
        Females[i] = 0;
      } else {
        Males[i] = 0;
        Females[i] = 1.0;
      }
    }

    LT = gsl_ran_discrete_preproc(N, All);
    LT_TH = gsl_ran_discrete_preproc(N, All);
    LT_males = gsl_ran_discrete_preproc(N, Males);
    LT_females = gsl_ran_discrete_preproc(N, Females);
  }

  void set_selfing(double s) {}

  void child_is_male(int i, bool male) { males_child[i] = male; }

  bool parent_is_male(int i) { return males_parent[i] == true; }
  void parent_is_male(int i, bool male) { males_parent[i] = male; }

  ~subpopulation_sexed() {
    G_parent.clear();
    G_child.clear();
    males_child.clear();
    males_parent.clear();
  }

  int draw_individual_by_sex(bool male) {
    return gsl_ran_discrete(rng, (male == true) ? LT_males : LT_females);
  }

  int draw_individual_threshold() { return gsl_ran_discrete(rng, LT_TH); }
  /**
   * @brief calculate fitnesses in parent population and create new lookup table
   * @param chr Chromosome population
   * @return void
   */
  void update_fitness(chromosome& chr) {

    gsl_ran_discrete_free(LT);
    int parent_size = (int)(G_parent.size() / 2);
    int n_males = 0;
    double All[parent_size];
    double Males[parent_size];
    double Females[parent_size];
    for (int i = 0; i < parent_size; i++) {
      All[i] = W(2 * i, 2 * i + 1, chr);
      if (males_parent[i] == true) {
        n_males++;
        Males[i] = All[i];
        Females[i] = 0;
      } else {
        Males[i] = 0;
        Females[i] = All[i];
      }
    }
    LT = gsl_ran_discrete_preproc(parent_size, All);
    LT_males = gsl_ran_discrete_preproc(parent_size, Males);
    LT_females = gsl_ran_discrete_preproc(parent_size, Females);
    if (T > 0) {
      update_threshold_fitness(All, parent_size);
    }
    if (S_ratio > 0) {
      update_sex_ratio_fitness(All, parent_size, n_males);
    }
  }

  void update_sex_ratio_fitness(double* All, int parent_size, int n_males) {

    double Females[parent_size];
    int n_females =
        subpopulation_sexed::get_reproductive_females(N, n_males, T, S_ratio);
    memset(Females, 0, sizeof(Females));

    int j = 0;
    int female;
    while (j < n_females) {
      female = draw_individual_by_sex(false);
      if (Females[female] == 0) {
        Females[female] = All[female];
        j++;
      }
    }
    LT_females = gsl_ran_discrete_preproc(parent_size, Females);
  }

  void swap() {
    G_child.swap(G_parent);
    G_child.resize(2 * N);
    males_parent.swap(males_child);
    males_child.resize(N);
  }
};

class subpopulation_hermaphrodite : public subpopulation {
  /**
   * @class subpopulation_hermaphrodite
   *  a subpopulation is described by the vector G of 2N genomes
   *  individual i is constituted by the two genomes 2*i and 2*i+1
   */

 public:
  subpopulation_hermaphrodite(int n) {
    N = n;
    S = 0.0;
    T = 0;
    G_parent.resize(2 * N);
    G_child.resize(2 * N);

    double A[N];
    for (int i = 0; i < N; i++) {
      A[i] = 1.0;
    }
    LT = gsl_ran_discrete_preproc(N, A);
    LT_males = LT;
  }

  void set_selfing(double s) { S = s; }

  bool parent_is_male(int i) { return true; }

  void child_is_male(int i, bool male) {}

  void parent_is_male(int i, bool male) {}

  ~subpopulation_hermaphrodite() {
    G_parent.clear();
    G_child.clear();
  }

  /**
   * @brief Retrieves an individual taking account the selected sex.
   *        In case of hermaphrodites is not necessary consider this,
   *        so is an alias of draw_individual
   */
  int draw_individual_by_sex(bool male) {
    return gsl_ran_discrete(rng, (male == true && T > 0) ? LT_males : LT);
  }

  /**
   * @brief calculate fitnesses in parent population and create new lookup table
   * @param chr Chromosome population
   * @return void
   */
  void update_fitness(chromosome& chr) {

    gsl_ran_discrete_free(LT);
    int parent_size = (int)(G_parent.size() / 2);
    double All[parent_size];
    for (int i = 0; i < parent_size; i++) {
      All[i] = W(2 * i, 2 * i + 1, chr);
    }
    LT = gsl_ran_discrete_preproc(parent_size, All);

    update_threshold_fitness(All, parent_size);
  }

  void swap() {
    G_child.swap(G_parent);
    G_child.resize(2 * N);
  }
};

class population : public map<int, subpopulation*> {
  /**
   * @class population
   * the population is a map of subpopulations
   */

 public:
  vector<substitution> Substitutions; /**< The population substitutions */
  map<int, subpopulation*>::iterator it;
  vector<string> parameters;

  int hermaphrodites;
  int sex_ratio;
  population() { hermaphrodites = 1; }
  /**
   * @brief add new empty subpopulation i of size N
   * @param i Population identifier
   * @param N Population number
   * @return void
   */
  void add_subpopulation(int i, unsigned int N) {

    if (count(i) != 0) {
      cerr << "ERROR (add subpopulation): subpopulation p" << i
           << " already exists" << endl;
      exit(1);
    }
    if (N < 1) {
      cerr << "ERROR (add subpopulation): subpopulation p" << i << " empty"
           << endl;
      exit(1);
    }

    if (hermaphrodites == 1) {
      insert(pair<int, subpopulation*>(i, new subpopulation_hermaphrodite(N)));
    } else {
      insert(pair<int, subpopulation*>(i, new subpopulation_sexed(N)));
    }
  }

  /**
   * @brief   add new subpopulation i of size N individuals drawn from source
   * population j
   * @param i New population id
   * @param j Populaition source
   * @param N Individuals drawn
   * @return void
   */

  void add_subpopulation(int i, int j, unsigned int N) {

    if (count(i) != 0) {
      cerr << "ERROR (add subpopulation): subpopulation p" << i
           << " already exists" << endl;
      exit(1);
    }
    if (count(j) == 0) {
      cerr << "ERROR (add subpopulation): source subpopulation p" << j
           << " does not exists" << endl;
      exit(1);
    }
    if (N < 1) {
      cerr << "ERROR (add subpopulation): subpopulation p" << i << " empty"
           << endl;
      exit(1);
    }

    if (hermaphrodites == 1) {
      insert(pair<int, subpopulation*>(i, new subpopulation_hermaphrodite(N)));
    } else {
      insert(pair<int, subpopulation*>(i, new subpopulation_sexed(N)));
    }

    for (int p = 0; p < find(i)->second->N; p++) {
      // draw individual from subpopulation j and assign to be a parent in i

      int m = find(j)->second->draw_individual();

      find(i)->second->G_parent[2 * p] = find(j)->second->G_parent[2 * m];
      find(i)->second->G_parent[2 * p + 1] =
          find(j)->second->G_parent[2 * m + 1];
      bool male = find(j)->second->parent_is_male(m);
      find(i)->second->parent_is_male(p, male);
    }
  }
  /**
   * @Brief set size of subpopulation i to N
   * @param i Population identifier
   * @param N New size of this population, if is 0, removes it.
   * @return void
   */
  void set_size(int i, unsigned int N) {

    if (count(i) == 0) {
      cerr << "ERROR (change size): no subpopulation p" << i << endl;
      exit(1);
    }

    if (N == 0)  // remove subpopulation i
    {
      erase(i);
      for (it = begin(); it != end(); it++) {
        it->second->m.erase(i);
      }
    } else {
      find(i)->second->N = N;
      find(i)->second->G_child.resize(2 * N);
    }
  }

  /**
   * @brief set fraction s of i that reproduces by selfing
   * @param i Population id
   * @param s Selfing fraction, must be a double between 0 and 1
   */
  void set_selfing(int i, double s) {

    if (count(i) == 0) {
      cerr << "ERROR (set selfing): no subpopulation p" << i << endl;
      exit(1);
    }
    if (s < 0.0 || s > 1.0) {
      cerr << "ERROR (set selfing): selfing fraction has to be within [0,1]"
           << endl;
      exit(1);
    }

    if (hermaphrodites == 0) {
      cerr << "ERROR (set selfing): Can't set selfing fraction with "
              "no hermaphrodites population" << endl;
      exit(1);
    }

    find(i)->second->S = s;
  }

  /**
   * @brief  set fraction m of i that originates as migrants from j per
   * generation
   * @param destination Destination Population Id
   * @param source Source population Id
   * @param ratio Fraction of migrants per generation, between 0 and 1
   * @return void
   */
  void set_migration(int destination, int source, double ratio) {

    if (count(destination) == 0) {
      cerr << "ERROR (set migration): no subpopulation p" << destination
           << endl;
      exit(1);
    }
    if (count(source) == 0) {
      cerr << "ERROR (set migration): no subpopulation p" << source << endl;
      exit(1);
    }
    if (ratio < 0.0 || ratio > 1.0) {
      cerr << "ERROR (set migration): migration fraction has to be within [0,1]"
           << endl;
      exit(1);
    }

    if (find(destination)->second->m.count(source) != 0) {
      find(destination)->second->m.erase(source);
    }

    find(destination)->second->m.insert(pair<int, double>(source, ratio));
  }

  /**
   * Defines the reproduction threshold for a population.
   * The thresold is like chosing a group of individuals for being the parents
   * of the next generation
   * @param population The subpopulation id
   * @param threshold  The new population threshold
   */
  void set_threshold(int subpopulation, int threshold) {
    if (count(subpopulation) == 0) {
      cerr << "ERROR (set threshold): no subpopulation source " << subpopulation
           << endl;
      exit(1);
    }
    if (threshold < 0) {
      cerr << "ERROR (set threshold): migration threshold has to be equal or "
              "greater than 1" << endl;
      exit(1);
    }
    // If threshold is greater than population means we want cross all
    // individuals.
    find(subpopulation)->second->T =
        (find(subpopulation)->second->N > threshold ? threshold : 0);
  }

  /**
   * Defines the sex threshold of this non-hermaphrodite population
   * The ratio is the definition of the number of females per reproductive male
   * Only available for non hremaprhrodites
   * @param population The subpopulation id
   * @param ratio The sex ratio
   */
  void set_sex_ratio(int subpopulation, double ratio) {
    if (hermaphrodites == 1) {
      cerr << "ERROR (set sex_ratio): Can't set sex ratio with "
              "hermaphrodites population" << endl;
      exit(1);
    }
    find(subpopulation)->second->S_ratio = ratio;
  }

  /**
   * Executes the main events of the system (new population, migrations, etc)
   * @param E an event
   * @param g the current generation
   * @param chr The current chromosome host of this event
   * @param FM tracked mutation-types
   * @return void
   */
  void execute_event(event& E, int g, chromosome& chr, vector<int>& FM) {
    char type = E.t; /* Type of event (P, N, S...) */

    if (type == 'P')  // add subpopulation
    {
      if (E.np == 2)  // empty subpopulation
      {
        string sub = E.s[0];
        sub.erase(0, 1);

        int i = atoi(sub.c_str());
        int n = (int)atof(E.s[1].c_str());
        add_subpopulation(i, n);
      }

      if (E.np == 3)  // drawn from source population
      {
        string sub1 = E.s[0];
        sub1.erase(0, 1);
        string sub2 = E.s[2];
        sub2.erase(0, 1);

        int i = atoi(sub1.c_str());
        int j = atoi(sub2.c_str());
        int n = (int)atof(E.s[1].c_str());
        add_subpopulation(i, j, n);
      }
    }

    if (type == 'N')  // set subpopulation size
    {
      string sub = E.s[0];
      sub.erase(0, 1);

      int i = atoi(sub.c_str());
      int n = (int)atof(E.s[1].c_str());

      set_size(i, n);
    }

    if (type == 'S')  // set selfing rate
    {
      string sub = E.s[0];
      sub.erase(0, 1);

      int i = atoi(sub.c_str());
      double s = atof(E.s[1].c_str());

      set_selfing(i, s);
    }

    if (type == 'M')  // change migration rate
    {
      string sub1 = E.s[0];
      sub1.erase(0, 1);
      string sub2 = E.s[1];
      sub2.erase(0, 1);

      int destination = atoi(sub1.c_str());
      int source = atoi(sub2.c_str());
      double rate = atof(E.s[2].c_str());

      set_migration(destination, source, rate);
    }

    if (type == 't')  // set threshold rate, means compute the threshold again
    {
      string sub1 = E.s[0];
      sub1.erase(0, 1);
      int subpopulation = atoi(sub1.c_str());
      int threshold = atoi(E.s[1].c_str());
      set_threshold(subpopulation, threshold);
    }

    if (type == 'r')  // change sex ratio
    {
      string sub1 = E.s[0];
      sub1.erase(0, 1);
      int subpopulation = atoi(sub1.c_str());
      double ratio = atof(E.s[1].c_str());
      set_sex_ratio(subpopulation, ratio);
    }

    if (type == 'A')  // output state of entire population
    {
      if (E.s.size() == 0) {
        cout << "#OUT: " << g << " A" << endl;
        print_all(chr);
      }
      if (E.s.size() == 1) {
        ofstream outfile;
        outfile.open(E.s[0].c_str());

        for (int i = 0; i < parameters.size(); i++) {
          outfile << parameters[i] << endl;
        }

        if (outfile.is_open()) {
          outfile << "#OUT: " << g << " A " << E.s[0].c_str() << endl;
          print_all(outfile, chr);
          outfile.close();
        } else {
          cerr << "ERROR (output): could not open " << E.s[0].c_str() << endl;
          exit(1);
        }
      }
    }

    if (type == 'R')  // output random subpopulation sample
    {
      string sub = E.s[0];
      sub.erase(0, 1);

      int i = atoi(sub.c_str());
      int n = atoi(E.s[1].c_str());
      cout << "#OUT: " << g << " R p" << i << " " << n << endl;

      if (E.s.size() == 3 && E.s[2] == "MS") {
        print_sample_ms(i, n, chr);
      } else {
        print_sample(i, n, chr);
      }
    }

    if (type == 'F')  // output list of fixed mutations
    {
      cout << "#OUT: " << g << " F " << endl;
      cout << "Mutations:" << endl;
      for (int i = 0; i < Substitutions.size(); i++) {
        cout << i + 1;
        Substitutions[i].print(chr);
      }
    }

    if (type == 'T')  // track mutation-types
    {
      string sub = E.s[0];
      sub.erase(0, 1);
      FM.push_back(atoi(sub.c_str()));
    }
  }
  /**
   * @brief introduce user-defined mutation
   * @param M introduced mutation object (a new mutation)
   * @param chr Chromosome host of the mutation
   * @return void
   */
  void introduce_mutation(introduced_mutation M, chromosome& chr) {

    if (count(M.i) == 0) {
      cerr << "ERROR (predetermined mutation): subpopulation " << M.i
           << " does not exists" << endl;
      exit(1);
    }
    if (chr.mutation_types.count(M.t) == 0) {
      cerr << "ERROR (predetermined mutation): mutation type m" << M.t
           << " has not been defined" << endl;
      exit(1);
    }
    if (find(M.i)->second->G_child.size() / 2 < M.nAA + M.nAa) {
      cerr << "ERROR (predetermined mutation): not enough individuals in "
              "subpopulation " << M.i << endl;
      exit(1);
    }

    mutation m;

    m.x = M.x;
    m.t = M.t;
    m.s = chr.mutation_types.find(M.t)->second.draw_s();
    m.i = M.i;
    m.g = M.g;

    // introduce homozygotes

    for (int j = 0; j < M.nAA; j++) {
      genome* g1 = &find(M.i)->second->G_child[2 * j];
      genome* g2 = &find(M.i)->second->G_child[2 * j + 1];
      (*g1).push_back(m);
      (*g2).push_back(m);
      sort((*g1).begin(), (*g1).end());
      sort((*g2).begin(), (*g2).end());
      (*g1).erase(unique((*g1).begin(), (*g1).end()), (*g1).end());
      (*g2).erase(unique((*g2).begin(), (*g2).end()), (*g2).end());
    }

    // introduce heterozygotes

    for (int j = M.nAA; j < M.nAA + M.nAa; j++) {
      genome* g1 = &find(M.i)->second->G_child[2 * j];
      (*g1).push_back(m);
      sort((*g1).begin(), (*g1).end());
      (*g1).erase(unique((*g1).begin(), (*g1).end()), (*g1).end());
    }
  }

  /**
   * @brief output trajectories of followed mutations and set s=0 for partial
   * sweeps
   * @param g Generation id
   * @param TM Mutations to track
   * @param PS partial_sweep to track
   * @param chr Chromosome host of the mutation
   */
  void track_mutations(int g, vector<int>& TM, vector<partial_sweep>& PS,
                       chromosome& chr) {

    // find all polymorphism of the types that are to be tracked

    for (it = begin(); it != end(); it++)  // go through all subpopulations
    {
      multimap<int, polymorphism> P;
      multimap<int, polymorphism>::iterator P_it;

      for (int i = 0; i < 2 * it->second->N; i++)  // go through all children
      {
        for (int k = 0; k < it->second->G_child[i].size();
             k++)  // go through all mutations
        {
          for (int j = 0; j < TM.size(); j++) {
            if (it->second->G_child[i][k].t == TM[j]) {
              add_mut(P, it->second->G_child[i][k]);
            }
          }
        }
      }

      // out put the frequencies of these mutations in each subpopulation

      for (P_it = P.begin(); P_it != P.end(); P_it++) {
        cout << "#OUT: " << g << " T p" << it->first << " ";
        P_it->second.print_no_id(P_it->first, chr);
      }
    }

    // check partial sweeps

    multimap<int, polymorphism> P;
    multimap<int, polymorphism>::iterator P_it;

    if (PS.size() > 0) {
      P.clear();

      int N = 0;
      for (it = begin(); it != end(); it++) {
        N += it->second->N;
      }

      // find all polymorphism that are supposed to undergo partial sweeps

      for (it = begin(); it != end(); it++)  // go through all subpopulations
      {
        for (int i = 0; i < 2 * it->second->N; i++)  // go through all children
        {
          for (int k = 0; k < it->second->G_child[i].size();
               k++)  // go through all mutations
          {
            for (int j = 0; j < PS.size(); j++) {
              if (it->second->G_child[i][k].x == PS[j].x &&
                  it->second->G_child[i][k].t == PS[j].t) {
                add_mut(P, it->second->G_child[i][k]);
              }
            }
          }
        }
      }

      // check whether a partial sweep has reached its target frequency

      for (P_it = P.begin(); P_it != P.end(); P_it++) {
        for (int j = 0; j < PS.size(); j++) {
          if (P_it->first == PS[j].x && P_it->second.t == PS[j].t) {
            if (((float)P_it->second.n) / (2 * N) >= PS[j].p) {
              // sweep has reached target frequency, set all s to zero

              for (it = begin(); it != end();
                   it++)  // go through all subpopulations
              {
                for (int i = 0; i < 2 * it->second->N;
                     i++)  // go through all children
                {
                  for (int k = 0; k < it->second->G_child[i].size();
                       k++)  // go through all mutations
                  {
                    if (it->second->G_child[i][k].x == PS[j].x &&
                        it->second->G_child[i][k].t == PS[j].t) {
                      it->second->G_child[i][k].s = 0.0;
                    }
                  }
                }
              }
              PS.erase(PS.begin() + j);
              j--;
            }
          }
        }
      }
    }
  }
  /**
   * @brief Evolves subpopulation a new generation.
   *        First of all are included the
   *        migrant populations. Later, the rest of the population group is
   *included.
   *        For example, when there are 2 subpopulations P and Q, with 10
   *individuals
   *        each one, if the population P has a migration rate of 0.2 from
   *population
   *        Q, two of the new individuals of population P will be chosen
   *randomly (taking
   *        care of fitness), later, the other 8 are individuals will be
   *reproduced from
   *        population P.
   *
   * @param i Subpopulation id
   * @param chr Chromosome related
   * @param g Generation number
   * @return void
   */
  void evolve_subpopulation(int i, chromosome& chr, int g) {
    int g1, g2, p1, p2, n_mut_1, n_mut_2;

    // create map of shuffled children ids

    int child_map[find(i)->second->N];
    for (int j = 0; j < find(i)->second->N; j++) {
      child_map[j] = j;
    }
    gsl_ran_shuffle(rng, child_map, find(i)->second->N, sizeof(int));

    int child_count = 0;  // counter over all N children (will get mapped to
    // child_map[child_count])

    // draw number of migrant individuals

    map<int, double>::iterator it;

    double m_rates[find(i)->second->m.size() + 1];
    unsigned int n_migrants[find(i)->second->m.size() + 1];

    double m_sum = 0.0;
    int pop_count = 0;

    for (map<int, double>::iterator it = find(i)->second->m.begin();
         it != find(i)->second->m.end(); it++) {
      m_rates[pop_count] = it->second;
      m_sum += it->second;
      pop_count++;
    }

    if (m_sum <= 1.0) {
      m_rates[pop_count] = 1.0 - m_sum;
    } else {
      cerr
          << "ERROR (evolve subpopulation): too many migrants in subpopulation "
          << i << endl;
      exit(1);
    }

    gsl_ran_multinomial(rng, find(i)->second->m.size() + 1, find(i)->second->N,
                        m_rates, n_migrants);

    // loop over all migration source populations

    pop_count = 0;

    for (map<int, double>::iterator it = find(i)->second->m.begin();
         it != find(i)->second->m.end(); it++) {
      int migrant_count = 0;

      while (migrant_count < n_migrants[pop_count]) {
        g1 = 2 * child_map[child_count];      // child genome 1
        g2 = 2 * child_map[child_count] + 1;  // child genome 2

        // draw parents in source population
        p1 = gsl_rng_uniform_int(rng,
                                 find(it->first)->second->G_parent.size() / 2);

        if (gsl_rng_uniform(rng) <
            find(it->first)->second->S) {  // Check if is at selfing rate
          p2 = p1;
        } else {
          p2 = gsl_rng_uniform_int(
              rng, find(it->first)->second->G_parent.size() / 2);
        }

        // recombination, gene-conversion, mutation
        crossover_mutation(i, g1, it->first, 2 * p1, 2 * p1 + 1, chr, g);
        crossover_mutation(i, g2, it->first, 2 * p2, 2 * p2 + 1, chr, g);
        find(it->first)->second->child_is_male(
            child_count, gsl_ran_binomial(rng, 0.5, 1) == 1);

        migrant_count++;
        child_count++;
      }
      pop_count++;
    }

    // remainder
    while (child_count < find(i)->second->N) {

      g1 = 2 * child_map[child_count];      // child genome 1
      g2 = 2 * child_map[child_count] + 1;  // child genome 2

      p1 =
          find(i)->second->draw_individual_by_sex(true);  // parent 1, male role

      if (gsl_rng_uniform(rng) < find(i)->second->S) {
        p2 = p1;
      }  // parent 2
      else {
        p2 = find(i)->second->draw_individual_by_sex(
            false);  // parent 1, female role
      }

      crossover_mutation(i, g1, i, 2 * p1, 2 * p1 + 1, chr, g);
      crossover_mutation(i, g2, i, 2 * p2, 2 * p2 + 1, chr, g);
      find(i)->second->child_is_male(child_count,
                                     gsl_ran_binomial(rng, 0.5, 1) == 1);
      child_count++;
    }
  }

  /**
   *
   * @param i subpopulation id destination
   * @param c genome c
   * @param j subpopulation id source
   * @param P1 genome from population 1
   * @param P2 genome from population 2
   * @param chr Chromosome
   * @param g generation number
   * @return void
   * child genome c in subpopulation i is assigned outcome of cross-overs at
   * breakpoints r
   * between parent genomes p1 and p2 from subpopulation j and new mutations
   * added
   *
   *  example R = (r1,r2)
   *
   * mutations (      x < r1) assigned from p1
   * mutations (r1 <= x < r2) assigned from p2
   * mutations (r2 <= x     ) assigned from p1
   *
   * p1 and p2 are swapped in half of the cases to assure random assortement.
   *   This swapping is only done without threshold.
   */
  void crossover_mutation(int i, int c, int j, int P1, int P2, chromosome& chr,
                          int g) {
    if (gsl_rng_uniform_int(rng, 2) == 0) {
      int swap = P1;
      P1 = P2;
      P2 = swap;
    }  // swap p1 and p2 only without threshold.

    find(i)->second->G_child[c].clear();

    // create vector with the mutations to be added

    vector<mutation> M;
    int n_mut = chr.draw_n_mut();
    for (int k = 0; k < n_mut; k++) {
      M.push_back(chr.draw_new_mut(j, g));
    }
    sort(M.begin(), M.end());

    // create vector with recombination breakpoints

    vector<int> R = chr.draw_breakpoints();
    R.push_back(chr.L + 1);
    sort(R.begin(), R.end());
    R.erase(unique(R.begin(), R.end()), R.end());

    vector<mutation>::iterator p1;
    vector<mutation>::iterator p1_max;
    vector<mutation>::iterator p2;
    vector<mutation>::iterator p2_max;
    // TODO: Refactorize me!

    p1 = find(j)->second->G_parent[P1].begin();
    p1_max = find(j)->second->G_parent[P1].end();
    p2 = find(j)->second->G_parent[P2].begin();
    p2_max = find(j)->second->G_parent[P2].end();

    vector<mutation>::iterator m = M.begin();
    vector<mutation>::iterator m_max = M.end();

    vector<mutation>::iterator p = p1;
    vector<mutation>::iterator p_max = p1_max;

    int r = 0;
    int r_max = R.size();
    int n = 0;
    bool present;

    while (r != r_max) {

      while ((p != p_max && (*p).x < R[r]) || (m != m_max && (*m).x < R[r])) {
        while (p != p_max && (*p).x < R[r] &&
               (m == m_max || (*p).x <= (*m).x)) {
          present = 0;
          if (n != 0 && find(i)->second->G_child[c].back().x == (*p).x) {
            int k = n - 1;

            while (present == 0 && k >= 0) {
              if (find(i)->second->G_child[c][k] == (*p)) {
                present = 1;
              }
              k--;
            }
          }

          if (present == 0) {
            find(i)->second->G_child[c].push_back(*p);
            n++;
          }
          p++;
        }
        while (m != m_max && (*m).x < R[r] &&
               (p == p_max || (*m).x <= (*p).x)) {
          present = 0;
          if (n != 0 && find(i)->second->G_child[c].back().x == (*m).x) {
            int k = n - 1;
            while (present == 0 && k >= 0) {
              if (find(i)->second->G_child[c][k] == (*m)) {
                present = 1;
              }
              k--;
            }
          }
          if (present == 0) {
            find(i)->second->G_child[c].push_back(*m);
            n++;
          }
          m++;
        }
      }
      // swap parents

      p1 = p2;
      p1_max = p2_max;
      p2 = p;
      p2_max = p_max;
      p = p1;
      p_max = p1_max;

      while (p != p_max && (*p).x < R[r]) {
        p++;
      }

      r++;
    }
  }

  /**
   * @brief find and remove fixed mutations from the children in all
   * subpopulations
   * @param g Generation id
   * @param chr chromosome
   */
  void swap_generations(int g, chromosome& chr) {
    remove_fixed(g);

    // make children the new parents and update fitnesses

    for (it = begin(); it != end(); it++) {
      it->second->swap();
      it->second->update_fitness(chr);
    }
  }
  /**
   * @brief find mutations that are fixed in all child subpopulations and return
   * vector with their ids
   *        (possibly the doc may be wrong)
   * @param int g generation number
   * return void
   */
  void remove_fixed(int g) {

    genome G = begin()->second->G_child[0];

    for (it = begin(); it != end(); it++)  // subpopulations
    {
      for (int i = 0; i < 2 * it->second->N; i++)  // child genomes
      {
        G = fixed(it->second->G_child[i], G);
      }
    }

    if (G.size() > 0) {
      for (it = begin(); it != end(); it++)  // subpopulations
      {
        for (int i = 0; i < 2 * it->second->N; i++)  // child genomes
        {
          it->second->G_child[i] = polymorphic(it->second->G_child[i], G);
        }
      }
      for (int i = 0; i < G.size(); i++) {
        Substitutions.push_back(substitution(G[i], g));
      }
    }
  }

  /**
   * @brief print all mutations and all genomes
   * @param chr Chromosome to print
   */
  void print_all(chromosome& chr) {

    cout << "Populations:" << endl;
    for (it = begin(); it != end(); it++) {
      cout << "p" << it->first << " " << it->second->N << endl;
    }

    multimap<int, polymorphism> P;
    multimap<int, polymorphism>::iterator P_it;

    // add all polymorphisms

    for (it = begin(); it != end(); it++)  // go through all subpopulations
    {
      for (int i = 0; i < 2 * it->second->N; i++)  // go through all children
      {
        for (int k = 0; k < it->second->G_child[i].size();
             k++)  // go through all mutations
        {
          add_mut(P, it->second->G_child[i][k]);
        }
      }
    }

    cout << "Mutations:" << endl;

    for (P_it = P.begin(); P_it != P.end(); P_it++) {
      P_it->second.print(P_it->first, chr);
    }

    cout << "Genomes:" << endl;

    // print all genomes

    for (it = begin(); it != end(); it++)  // go through all subpopulations
    {
      for (int i = 0; i < 2 * it->second->N; i++)  // go through all children
      {
        cout << "p" << it->first << ":" << i + 1;

        for (int k = 0; k < it->second->G_child[i].size();
             k++)  // go through all mutations
        {
          int id = find_mut(P, it->second->G_child[i][k]);
          cout << " " << id;
        }
        cout << endl;
      }
    }
  }

  /**
   * @brief print all mutations and all genomes to a file
   * @param outfile Output file stream
   * @param chromosome
   * @return void
   */
  void print_all(ofstream& outfile, chromosome& chr) {

    outfile << "Populations:" << endl;
    for (it = begin(); it != end(); it++) {
      outfile << "p" << it->first << " " << it->second->N << endl;
    }

    multimap<int, polymorphism> P;
    multimap<int, polymorphism>::iterator P_it;

    // add all polymorphisms

    for (it = begin(); it != end(); it++)  // go through all subpopulations
    {
      for (int i = 0; i < 2 * it->second->N; i++)  // go through all children
      {
        for (int k = 0; k < it->second->G_child[i].size();
             k++)  // go through all mutations
        {
          add_mut(P, it->second->G_child[i][k]);
        }
      }
    }

    outfile << "Mutations:" << endl;

    for (P_it = P.begin(); P_it != P.end(); P_it++) {
      P_it->second.print(outfile, P_it->first, chr);
    }

    outfile << "Genomes:" << endl;

    // print all genomes

    for (it = begin(); it != end(); it++)  // go through all subpopulations
    {
      for (int i = 0; i < 2 * it->second->N; i++)  // go through all children
      {
        outfile << "p" << it->first << ":" << i + 1;

        for (int k = 0; k < it->second->G_child[i].size();
             k++)  // go through all mutations
        {
          int id = find_mut(P, it->second->G_child[i][k]);
          outfile << " " << id;
        }
        outfile << endl;
      }
    }
  }

  /**
   * @brief print sample of n genomes from subpopulation  i
   * @param i subpopulation to print
   * @param n number of genomes to print
   * @param chr Chromosome source
   * @param void
   */
  void print_sample(int i, int n, chromosome& chr) {

    if (count(i) == 0) {
      cerr << "ERROR (output): subpopulation p" << i << " does not exists"
           << endl;
      exit(1);
    }

    vector<int> sample;

    multimap<int, polymorphism> P;
    multimap<int, polymorphism>::iterator P_it;

    for (int s = 0; s < n; s++) {
      int j = gsl_rng_uniform_int(rng, find(i)->second->G_child.size());
      sample.push_back(j);

      for (int k = 0; k < find(i)->second->G_child[j].size();
           k++)  // go through all mutations
      {
        add_mut(P, find(i)->second->G_child[j][k]);
      }
    }

    cout << "Mutations:" << endl;

    for (P_it = P.begin(); P_it != P.end(); P_it++) {
      P_it->second.print(P_it->first, chr);
    }

    cout << "Genomes:" << endl;

    // print all genomes

    for (int j = 0; j < sample.size(); j++)  // go through all children
    {
      cout << "p" << find(i)->first << ":" << sample[j] + 1;

      for (int k = 0; k < find(i)->second->G_child[sample[j]].size();
           k++)  // go through all mutations
      {
        int id = find_mut(P, find(i)->second->G_child[sample[j]][k]);
        cout << " " << id;
      }
      cout << endl;
    }
  }
  /**
   * @brief print sample of n genomes from subpopulation  i
   * @param i Subpopulation id
   * @param n Number of genomes to print
   * @param chr Chromomsome id
   * @return void
   */
  void print_sample_ms(int i, int n, chromosome& chr) {

    if (count(i) == 0) {
      cerr << "ERROR (output): subpopulation p" << i << " does not exists"
           << endl;
      exit(1);
    }

    vector<int> sample;

    multimap<int, polymorphism> P;
    multimap<int, polymorphism>::iterator P_it;

    for (int s = 0; s < n; s++) {
      int j = gsl_rng_uniform_int(rng, find(i)->second->G_child.size());
      sample.push_back(j);

      for (int k = 0; k < find(i)->second->G_child[j].size();
           k++)  // go through all mutations
      {
        add_mut(P, find(i)->second->G_child[j][k]);
      }
    }

    // print header

    cout << endl << "//" << endl << "segsites: " << P.size() << endl;

    // print all positions

    if (P.size() > 0) {
      cout << "positions:";
      for (P_it = P.begin(); P_it != P.end(); P_it++) {
        cout << " " << fixed << setprecision(7)
             << (double)(P_it->first + 1) / (chr.L + 1);
      }
      cout << endl;
    }

    // print genotypes

    for (int j = 0; j < sample.size(); j++)  // go through all children
    {
      string genotype(P.size(), '0');

      for (int k = 0; k < find(i)->second->G_child[sample[j]].size();
           k++)  // go through all mutations
      {
        int pos = 0;
        mutation m = find(i)->second->G_child[sample[j]][k];

        for (P_it = P.begin(); P_it != P.end(); P_it++) {
          if (P_it->first == m.x && P_it->second.t == m.t &&
              P_it->second.s == m.s) {
            genotype.replace(pos, 1, "1");
            break;
          }
          pos++;
        }
      }
      cout << genotype << endl;
    }
  }
  /**
   * @brief find m in P and return its id
   * @param P Population list
   * @param m mutation object
   * @return integer Mutation id
   */
  int find_mut(multimap<int, polymorphism>& P, mutation m) {

    int id = 0;

    // iterate through all mutations with same position

    multimap<int, polymorphism>::iterator it;
    pair<multimap<int, polymorphism>::iterator,
         multimap<int, polymorphism>::iterator> range = P.equal_range(m.x);
    it = range.first;

    while (it != range.second) {
      if (it->second.t == m.t && it->second.s == m.s) {
        id = it->second.id;
        it = range.second;
      } else {
        it++;
      }
    }

    return id;
  }
  /**
   * @brief if mutation is present in P increase prevalence, otherwise add it
   * @param P Population
   * @param m Mutation object
   * @return void
   */
  void add_mut(multimap<int, polymorphism>& P, mutation m) {

    int id = 0;

    // iterate through all mutations with same position

    multimap<int, polymorphism>::iterator it;
    pair<multimap<int, polymorphism>::iterator,
         multimap<int, polymorphism>::iterator> range = P.equal_range(m.x);
    it = range.first;

    while (it != range.second) {
      if (it->second.t == m.t && it->second.s == m.s) {
        id = it->second.id;
        it->second.n++;
        it = range.second;
      } else {
        it++;
      }
    }

    // if not already present, add mutation to P

    if (id == 0) {
      id = P.size() + 1;
      P.insert(pair<int, polymorphism>(
          m.x, polymorphism(id, m.t, m.s, m.i, m.g, 1)));
    }
  }
};

void get_line(ifstream& infile, string& line) {
  getline(infile, line);
  if (line.find("/") != string::npos) {
    line.erase(line.find("/"));
  }                                            // remove all after "/"
  line.erase(0, line.find_first_not_of(' '));  // remove leading whitespaces
  line.erase(line.find_last_not_of(' ') + 1);  // remove trailing whitespaces
};

void input_error(int type, string line) {
  cerr << endl;

  if (type == -2)  // no population defined
  {
    cerr << "ERROR (parameter file): no population to simulate:" << endl
         << endl;
  } else if (type == -1)  // unknown parameter
  {
    cerr << "ERROR (parameter file): unknown parameter: " << line << endl
         << endl;
  } else if (type == 0)  // invalid parameter file
  {
    cerr << "ERROR (parameter file): could not open: " << line << endl << endl;
  } else if (type == 1)  // mutation rate
  {
    cerr << "ERROR (parameter file): invalid mutation rate: " << line << endl
         << endl;
    cerr << "Required syntax:" << endl << endl;
    cerr << "#MUTATION RATE" << endl;
    cerr << "<u>" << endl << endl;
    cerr << "Example:" << endl << endl;
    cerr << "#MUTATION RATE" << endl;
    cerr << "1.5e-8" << endl << endl;
  } else if (type == 2)  // mutation type
  {
    cerr << "ERROR (parameter file): invalid mutation type: " << line << endl
         << endl;
    cerr << "Required syntax:" << endl << endl;
    cerr << "#MUTATION TYPES" << endl;
    cerr << "<mutation-type-id> <h> <DFE-type> [DFE parameters]" << endl;
    cerr << "..." << endl << endl;
    cerr << "Example:" << endl << endl;
    cerr << "#MUTATION TYPES" << endl;
    cerr << "m1 0.2 g -0.05 0.2" << endl;
    cerr << "m2 0.0 f 0.0" << endl;
    cerr << "m3 0.5 e 0.01" << endl << endl;
  } else if (type == 3)  // genomic element type
  {
    cerr << "ERROR (parameter file): invalid genomic element type: " << line
         << endl << endl;
    cerr << "Required syntax:" << endl << endl;
    cerr << "#GENOMIC ELEMENT TYPES" << endl;
    cerr << "<element-type-id> <mut-type> <x> [<mut-type> <x>...]" << endl;
    cerr << "..." << endl << endl;
    cerr << "Example:" << endl << endl;
    cerr << "#GENOMIC ELEMENT TYPES" << endl;
    cerr << "g1 m3 0.8 m2 0.01 m1 0.19" << endl << endl;
  } else if (type == 4)  // chromosome organization
  {
    cerr << "ERROR (parameter file): invalid chromosome organization: " << line
         << endl << endl;
    cerr << "Required syntax:" << endl << endl;
    cerr << "#CHROMOSOME ORGANIZATION" << endl;
    cerr << "<element-type> <start> <end>" << endl;
    cerr << "..." << endl << endl;
    cerr << "Example:" << endl << endl;
    cerr << "#CHROMOSOME ORGANIZATION" << endl;
    cerr << "g1 1000 1999" << endl << endl;
  } else if (type == 5)  // recombination rate
  {
    cerr << "ERROR (parameter file): invalid recombination rate: " << line
         << endl << endl;
    cerr << "Required syntax:" << endl << endl;
    cerr << "#RECOMBINATION RATE" << endl;
    cerr << "<interval-end> <r>" << endl;
    cerr << "..." << endl << endl;
    cerr << "Example:" << endl << endl;
    cerr << "#RECOMBINATION RATE" << endl;
    cerr << "10000 1e-8" << endl;
    cerr << "20000 4.5e-8" << endl << endl;
  } else if (type == 6)  // generations
  {
    cerr << "ERROR (parameter file): invalid generations: " << line << endl
         << endl;
    cerr << "Required syntax:" << endl << endl;
    cerr << "#GENERATIONS" << endl;
    cerr << "<t> [<start>]" << endl << endl;
    cerr << "Example:" << endl << endl;
    cerr << "#GENERATIONS" << endl;
    cerr << "10000" << endl << endl;
  } else if (type == 7)  // demography and structure
  {
    cerr << "ERROR (parameter file): invalid demography and structure: " << line
         << endl << endl;
    cerr << "Required syntax:" << endl << endl;
    cerr << "#DEMOGRAPHY AND STRUCTURE" << endl;
    cerr << "<time> <event-type> [event parameters]" << endl;
    cerr << "..." << endl << endl;
    cerr << "Example:" << endl << endl;
    cerr << "DEMOGRAPHY AND STRUCTURE" << endl;
    cerr << "1 P p1 1000" << endl;
    cerr << "1 T p1 500" << endl;
    cerr << "1 S p1 0.05 /only when hermaphrodites" << endl;
    cerr << "1 R p1 0.5 /only when no hermaphrodites" << endl;
    cerr << "1000 P p2 100 p1" << endl;
    cerr << "1000 T p1 500" << endl;
    cerr << "1000 S p2 0.05 /only when hermaphrodites" << endl;
    cerr << "1000 R p1 0.5 /only when no hermaphrodites" << endl;
    cerr << "2000 N p1 1e4" << endl;
    cerr << "2000 M p2 p1 0.01" << endl << endl;
  } else if (type == 8)  // output
  {
    cerr << "ERROR (parameter file): invalid output: " << line << endl << endl;
    cerr << "Required syntax:" << endl << endl;
    cerr << "#OUTPUT" << endl;
    cerr << "<time> <output-type> [output parameters]" << endl;
    cerr << "..." << endl << endl;
    cerr << "Example:" << endl << endl;
    cerr << "OUTPUT" << endl;
    cerr << "2000 A outfile" << endl;
    cerr << "1000 R p1 10" << endl;
    cerr << "1000 R p1 10 MS" << endl;
    cerr << "2000 F" << endl;
    cerr << "1 T m3" << endl << endl;
  } else if (type == 9)  // initialization
  {
    cerr << "ERROR (parameter file): invalid initialization: " << line << endl
         << endl;
    cerr << "Required syntax:" << endl << endl;
    cerr << "#INITIALIZATION" << endl;
    cerr << "<filename>" << endl << endl;
    cerr << "Example:" << endl << endl;
    cerr << "#INITIALIZATION" << endl;
    cerr << "outfile" << endl << endl;
  } else if (type == 10)  // seed
  {
    cerr << "ERROR (parameter file): invalid seed: " << line << endl << endl;
    cerr << "Required syntax:" << endl << endl;
    cerr << "#SEED" << endl;
    cerr << "<seed>" << endl << endl;
    cerr << "Example:" << endl << endl;
    cerr << "#SEED" << endl;
    cerr << "141235" << endl << endl;
  } else if (type == 11)  // predetermined mutation
  {
    cerr << "ERROR (parameter file): invalid predetermined mutations: " << line
         << endl << endl;
    cerr << "Required syntax:" << endl << endl;
    cerr << "#PREDETERMINED MUTATIONS" << endl;
    cerr << "<time> <mut-type> <x> <pop> <nAA> <nAa>" << endl << endl;
    cerr << "Example:" << endl << endl;
    cerr << "#PREDETERMINED MUTATIONS" << endl;
    cerr << "5000 m7 45000 p1 0 1" << endl << endl;
  } else if (type == 12)  // gene conversion
  {
    cerr << "ERROR (parameter file): invalid gene conversion: " << line << endl
         << endl;
    cerr << "Required syntax:" << endl << endl;
    cerr << "#GENE CONVERSION" << endl;
    cerr << "<fraction> <average-length>" << endl << endl;
    cerr << "Example:" << endl << endl;
    cerr << "#GENE CONVERSION" << endl;
    cerr << "0.5 20" << endl << endl;
  }

  exit(1);
};

void check_input_file(char* file) {
  int mutation_types = 0;
  int mutation_rate = 0;
  int genomic_element_types = 0;
  int chromosome_organization = 0;
  int recombination_rate = 0;
  int generations = 0;
  int population = 0;

  ifstream infile(file);
  if (!infile.is_open()) {
    input_error(0, string(file));
  }

  string line;
  string sub;

  get_line(infile, line);

  while (!infile.eof()) {
    if (line.find('#') != string::npos) {
      if (line.find("MUTATION RATE") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            if (line.find_first_not_of("1234567890.e-") != string::npos) {
              input_error(1, line);
            } else {
              mutation_rate++;
            }
          }
          get_line(infile, line);
        }
      } else if (line.find("MUTATION TYPES") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            int good = 1;

            istringstream iss(line);
            iss >> sub;

            if (sub.compare(0, 1, "m") != 0) {
              good = 0;
            }
            sub.erase(0, 1);

            if (sub.find_first_not_of("1234567890") != string::npos) {
              good = 0;
            }  // id
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;
            if (sub.find_first_not_of("1234567890.-") != string::npos) {
              good = 0;
            }  // h
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;
            if (sub.find_first_not_of("fge") != string::npos) {
              good = 0;
            }  // DFE-type

            if (sub.compare("f") == 0 ||
                sub.compare("e") == 0)  // one parameter
            {
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;
              if (sub.find_first_not_of("1234567890.-") != string::npos) {
                good = 0;
              }
              if (!iss.eof()) {
                good = 0;
              }
            }
            if (sub.compare("g") == 0)  // two parameters
            {
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;
              if (sub.find_first_not_of("1234567890.-") != string::npos) {
                good = 0;
              }
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;
              if (sub.find_first_not_of("1234567890.-") != string::npos) {
                good = 0;
              }
              if (!iss.eof()) {
                good = 0;
              }
            }

            if (good == 0) {
              input_error(2, line);
            } else {
              mutation_types++;
            }
          }
          get_line(infile, line);
        }
      } else if (line.find("GENOMIC ELEMENT TYPES") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            int good = 1;

            istringstream iss(line);
            iss >> sub;

            if (sub.compare(0, 1, "p") == 0) {
              sub.erase(0, 1);

              if (sub.find_first_not_of("1234567890") != string::npos) {
                good = 0;
              }
              iss >> sub;
            }  // optional population id

            if (sub.compare(0, 1, "g") != 0) {
              good = 0;
            }
            sub.erase(0, 1);

            if (sub.find_first_not_of("1234567890") != string::npos) {
              good = 0;
            }  // id

            if (iss.eof()) {
              good = 0;
            }

            while (!iss.eof()) {
              iss >> sub;
              if (sub.compare(0, 1, "m") != 0) {
                good = 0;
              }
              sub.erase(0, 1);  // mutation type id
              if (sub.find_first_not_of("1234567890") != string::npos) {
                good = 0;
              }
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;
              if (sub.find_first_not_of("1234567890e.-") != string::npos) {
                good = 0;
              }  // fraction
            }

            if (good == 0) {
              input_error(3, line);
            } else {
              genomic_element_types++;
            }
          }
          get_line(infile, line);
        }
      } else if (line.find("CHROMOSOME ORGANIZATION") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            int good = 1;

            istringstream iss(line);
            iss >> sub;

            if (sub.compare(0, 1, "g") != 0) {
              good = 0;
            }
            sub.erase(0, 1);

            if (sub.find_first_not_of("1234567890") != string::npos) {
              good = 0;
            }  // id
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;
            if (sub.find_first_not_of("1234567890e") != string::npos) {
              good = 0;
            }  // start
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;
            if (sub.find_first_not_of("1234567890e") != string::npos) {
              good = 0;
            }  // end
            if (!iss.eof()) {
              good = 0;
            }

            if (good == 0) {
              input_error(4, line);
            } else {
              chromosome_organization++;
            }
          }
          get_line(infile, line);
        }
      } else if (line.find("RECOMBINATION RATE") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            int good = 1;

            istringstream iss(line);
            iss >> sub;

            if (sub.find_first_not_of("1234567890e") != string::npos) {
              good = 0;
            }  // end
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;
            if (sub.find_first_not_of("1234567890e.-") != string::npos) {
              good = 0;
            }  // rate
            if (!iss.eof()) {
              good = 0;
            }

            if (good == 0) {
              input_error(5, line);
            } else {
              recombination_rate++;
            }
          }
          get_line(infile, line);
        }
      } else if (line.find("GENE CONVERSION") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            int good = 1;

            istringstream iss(line);
            iss >> sub;

            if (sub.find_first_not_of("1234567890e.-") != string::npos) {
              good = 0;
            }  // fraction
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;
            if (sub.find_first_not_of("1234567890e.-") != string::npos) {
              good = 0;
            }  // average length
            if (!iss.eof()) {
              good = 0;
            }

            if (good == 0) {
              input_error(12, line);
            }
          }
          get_line(infile, line);
        }
      } else if (line.find("GENERATIONS") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            int good = 1;

            istringstream iss(line);
            iss >> sub;

            if (sub.find_first_not_of("1234567890e") != string::npos) {
              good = 0;
            }  // t
            if (!iss.eof()) {
              iss >> sub;
              if (sub.find_first_not_of("1234567890e") != string::npos) {
                good = 0;
              }  // start
            }
            if (!iss.eof()) {
              good = 0;
            }

            if (good == 0) {
              input_error(6, line);
            } else {
              generations++;
            }
          }
          get_line(infile, line);
        }
      } else if (line.find("DEMOGRAPHY AND STRUCTURE") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            int good = 1;

            istringstream iss(line);
            iss >> sub;

            if (sub.find_first_not_of("1234567890e") != string::npos) {
              good = 0;
            }  // t
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;
            if (sub.find_first_not_of("PSMNTR") != string::npos) {
              good = 0;
            }  // event type

            if (sub.compare("P") == 0 ||
                sub.compare("T") == 0)  // two or three positive integers
            {
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // p1
              if (sub.compare(0, 1, "p") != 0) {
                good = 0;
              }
              sub.erase(0, 1);
              if (sub.find_first_not_of("1234567890") != string::npos) {
                good = 0;
              }

              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // N
              if (sub.find_first_not_of("1234567890e") != string::npos) {
                good = 0;
              }

              if (!iss.eof())  // p2
              {
                iss >> sub;
                if (sub.compare(0, 1, "p") != 0) {
                  good = 0;
                }
                sub.erase(0, 1);
                if (sub.find_first_not_of("1234567890") != string::npos) {
                  good = 0;
                }
                if (!iss.eof()) {
                  good = 0;
                }
              }

              population++;
            }

            if (sub.compare("N") == 0)  // two positive integers
            {
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // p
              if (sub.compare(0, 1, "p") != 0) {
                good = 0;
              }
              sub.erase(0, 1);
              if (sub.find_first_not_of("1234567890") != string::npos) {
                good = 0;
              }
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // N
              if (sub.find_first_not_of("1234567890e") != string::npos) {
                good = 0;
              }
              if (!iss.eof()) {
                good = 0;
              }
            }

            if (sub.compare("S") == 0)  // one positive integer and a double
            {
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // p
              if (sub.compare(0, 1, "p") != 0) {
                good = 0;
              }
              sub.erase(0, 1);
              if (sub.find_first_not_of("1234567890") != string::npos) {
                good = 0;
              }
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // sigma
              if (sub.find_first_not_of("1234567890.-e") != string::npos) {
                good = 0;
              }
              if (!iss.eof()) {
                good = 0;
              }
            }

            if (sub.compare("M") == 0)  // two positive integers and a double
            {
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // p
              if (sub.compare(0, 1, "p") != 0) {
                good = 0;
              }
              sub.erase(0, 1);
              if (sub.find_first_not_of("1234567890") != string::npos) {
                good = 0;
              }
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // p
              if (sub.compare(0, 1, "p") != 0) {
                good = 0;
              }
              sub.erase(0, 1);
              if (sub.find_first_not_of("1234567890") != string::npos) {
                good = 0;
              }
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // M
              if (sub.find_first_not_of("1234567890.-e") != string::npos) {
                good = 0;
              }
              if (!iss.eof()) {
                good = 0;
              }
            }

            if (sub.compare("R") == 0)  // A double
            {
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // p
              if (sub.compare(0, 1, "p") != 0) {
                good = 0;
              }
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // ratio
              if (sub.find_first_not_of("1234567890.") != string::npos) {
                good = 0;
              }
              if (!iss.eof()) {
                good = 0;
              }
            }

            if (good == 0) {
              input_error(7, line);
            }
          }
          get_line(infile, line);
        }
      } else if (line.find("OUTPUT") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            int good = 1;

            istringstream iss(line);
            iss >> sub;

            if (sub.find_first_not_of("1234567890e") != string::npos) {
              good = 0;
            }  // t
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;
            if (sub.find_first_not_of("ARFT") != string::npos) {
              good = 0;
            }  // event type

            if (sub.compare("A") == 0)  // no parameter of filename
            {
              if (!iss.eof()) {
                iss >> sub;
                if (!iss.eof()) {
                  good = 0;
                }
              }
            }

            if (sub.compare("R") == 0)  // two positive integers
            {
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // p1
              if (sub.compare(0, 1, "p") != 0) {
                good = 0;
              }
              sub.erase(0, 1);
              if (sub.find_first_not_of("1234567890") != string::npos) {
                good = 0;
              }
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // p2
              if (sub.find_first_not_of("1234567890") != string::npos) {
                good = 0;
              }
              if (!iss.eof()) {
                iss >> sub;  // MS
                if (sub != "MS") {
                  good = 0;
                }
              }
              if (!iss.eof()) {
                good = 0;
              }
            }

            if (sub.compare("F") == 0)  // no parameter
            {
              if (!iss.eof()) {
                good = 0;
              }
            }

            if (good == 0) {
              input_error(8, line);
            }
          }
          get_line(infile, line);
        }
      } else if (line.find("INITIALIZATION") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            int good = 1;

            istringstream iss(line);
            iss >> sub;
            if (!iss.eof()) {
              good = 0;
            }

            if (good == 0) {
              input_error(9, line);
            }

            population++;
          }
          get_line(infile, line);
        }
      } else if (line.find("SEED") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            int good = 1;

            istringstream iss(line);
            iss >> sub;
            if (sub.find_first_not_of("1234567890-") != string::npos) {
              good = 0;
            }
            if (!iss.eof()) {
              good = 0;
            }

            if (good == 0) {
              input_error(10, line);
            }
          }
          get_line(infile, line);
        }
      } else if (line.find("HERMAPHRODITES") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            int good = 1;

            istringstream iss(line);
            iss >> sub;
            if (sub.find_first_not_of("10") != string::npos) {
              good = 0;
            }
            if (!iss.eof()) {
              good = 0;
            }

            if (good == 0) {
              input_error(10, line);
            }
          }
          get_line(infile, line);
        }

      } else if (line.find("PREDETERMINED MUTATIONS") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            int good = 1;

            istringstream iss(line);
            iss >> sub;  // time
            if (sub.find_first_not_of("1234567890e") != string::npos) {
              good = 0;
            }
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;  // id
            if (sub.compare(0, 1, "m") != 0) {
              good = 0;
            }
            sub.erase(0, 1);
            if (sub.find_first_not_of("1234567890") != string::npos) {
              good = 0;
            }
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;  // x
            if (sub.find_first_not_of("1234567890e") != string::npos) {
              good = 0;
            }
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;  // sub
            if (sub.compare(0, 1, "p") != 0) {
              good = 0;
            }
            sub.erase(0, 1);
            if (sub.find_first_not_of("1234567890") != string::npos) {
              good = 0;
            }
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;  // nAA
            if (sub.find_first_not_of("1234567890") != string::npos) {
              good = 0;
            }
            if (iss.eof()) {
              good = 0;
            }
            iss >> sub;  // nAa
            if (sub.find_first_not_of("1234567890") != string::npos) {
              good = 0;
            }

            if (!iss.eof()) {
              iss >> sub;
              if (sub.find_first_not_of("P") != string::npos) {
                good = 0;
              }
              if (iss.eof()) {
                good = 0;
              }
              iss >> sub;  // freq
              if (sub.find_first_not_of("1234567890.-e") != string::npos) {
                good = 0;
              }
            }

            if (!iss.eof()) {
              good = 0;
            }

            if (good == 0) {
              input_error(11, line);
            }
          }
          get_line(infile, line);
        }
      } else {
        input_error(-1, line);
      }
    } else {
      get_line(infile, line);
    }
  }

  if (mutation_rate != 1) {
    input_error(1, string());
  }
  if (mutation_types < 1) {
    input_error(2, string());
  }
  if (genomic_element_types < 1) {
    input_error(3, string());
  }
  if (chromosome_organization < 1) {
    input_error(4, string());
  }
  if (recombination_rate < 1) {
    input_error(5, string());
  }
  if (generations < 1) {
    input_error(6, string());
  }
  if (population < 1) {
    input_error(-2, string());
  }
};

void initialize_from_file(population& P, const char* file, chromosome& chr) {
  // initialize population from file

  map<int, mutation> M;

  string line;
  string sub;

  ifstream infile(file);

  if (!infile.is_open()) {
    cerr << "ERROR (initialize): could not open initialization file" << endl;
    exit(1);
  }

  get_line(infile, line);

  while (line.find("Populations") == string::npos && !infile.eof()) {
    get_line(infile, line);
  }

  get_line(infile, line);

  while (line.find("Mutations") == string::npos && !infile.eof()) {
    istringstream iss(line);
    iss >> sub;
    sub.erase(0, 1);
    int i = atoi(sub.c_str());
    iss >> sub;
    int n = atoi(sub.c_str());
    P.add_subpopulation(i, n);
    get_line(infile, line);
  }

  get_line(infile, line);

  while (line.find("Genomes") == string::npos && !infile.eof()) {
    istringstream iss(line);
    iss >> sub;
    int id = atoi(sub.c_str());
    iss >> sub;
    sub.erase(0, 1);
    int t = atoi(sub.c_str());
    iss >> sub;
    int x = atoi(sub.c_str()) - 1;
    iss >> sub;
    float s = atof(sub.c_str());
    iss >> sub;
    iss >> sub;
    sub.erase(0, 1);
    int i = atoi(sub.c_str());
    iss >> sub;
    int g = atoi(sub.c_str());

    M.insert(pair<int, mutation>(id, mutation(t, x, s, i, g)));
    get_line(infile, line);
  }

  get_line(infile, line);

  while (!infile.eof()) {
    istringstream iss(line);
    iss >> sub;
    sub.erase(0, 1);
    int pos = sub.find_first_of(":");
    int p = atoi(sub.substr(0, pos + 1).c_str());
    sub.erase(0, pos + 1);
    int i = atoi(sub.c_str());

    while (iss >> sub) {
      int id = atoi(sub.c_str());
      P.find(p)->second->G_parent[i - 1].push_back(M.find(id)->second);
    }
    get_line(infile, line);
  }

  for (P.it = P.begin(); P.it != P.end(); P.it++) {
    P.it->second->update_fitness(chr);
  }
};

void initialize(population& P, char* file, chromosome& chr, int& t_start,
                int& t_duration, multimap<int, event>& E,
                multimap<int, event>& O, multimap<int, introduced_mutation>& IM,
                vector<partial_sweep>& PS, vector<string>& parameters) {
  string line;
  string sub;
  vector<int> subpopulations;
  ifstream infile(file);

  int hermaphrodites = 1;

  long pid = getpid();
  time_t* tp, t;
  tp = &t;
  time(tp);
  t += pid;
  int seed = t;

  get_line(infile, line);

  while (!infile.eof()) {
    if (line.find('#') != string::npos) {
      if (line.find("MUTATION RATE") != string::npos) {
        get_line(infile, line);
        parameters.push_back("#MUTATION RATE");
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            parameters.push_back(line);
            istringstream iss(line);
            iss >> sub;
            chr.M = atof(sub.c_str());
          }
          get_line(infile, line);
        }
      } else if (line.find("MUTATION TYPES") != string::npos) {
        get_line(infile, line);
        parameters.push_back("#MUTATION TYPES");
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            parameters.push_back(line);

            // FORMAT: i h t p1 p2 ... (identifier, dominance coefficient DFE
            // type, DFE parameters)

            int i;
            float h;
            char t;
            vector<double> p;
            istringstream iss(line);
            iss >> sub;
            sub.erase(0, 1);
            i = atoi(sub.c_str());

            if (chr.mutation_types.count(i) > 0) {
              cerr << "ERROR (initialize): mutation type " << i
                   << " already defined" << endl;
              exit(1);
            }

            iss >> sub;
            h = atof(sub.c_str());
            iss >> sub;
            t = sub.at(0);
            while (iss >> sub) {
              p.push_back(atof(sub.c_str()));
            }
            chr.mutation_types.insert(
                pair<int, mutation_type>(i, mutation_type(h, t, p)));
          }
          get_line(infile, line);
        }
      } else if (line.find("GENOMIC ELEMENT TYPES") != string::npos) {
        get_line(infile, line);
        parameters.push_back("#GENOMIC ELEMENT TYPES");
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            parameters.push_back(line);
            genomic_element_types_line genomic_line;
            genomic_line = translate_genomic_elements_type_line(line);
            pair<int, int> key =
                make_pair(genomic_line.pop_id, genomic_line.genomic_id);
            if (chr.genomic_element_types.count(key) > 0) {
              cerr << "ERROR (initialize): genomic element type "
                   << genomic_line.genomic_id << " already defined" << endl;
              exit(1);
            }

            chr.genomic_element_types.insert(
                pair<pair<int, int>, genomic_element_type>(
                    key, genomic_element_type(genomic_line.mutations,
                                              genomic_line.fractions)));
          }
          get_line(infile, line);
        }
      } else if (line.find("CHROMOSOME ORGANIZATION") != string::npos) {
        get_line(infile, line);
        parameters.push_back("#CHROMOSOME ORGANIZATION");
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            parameters.push_back(line);

            // FORMAT: i s e (genomic element type identifier, start, end)

            int i, s, e;
            istringstream iss(line);
            iss >> sub;
            sub.erase(0, 1);
            i = atoi(sub.c_str());
            iss >> sub;
            s = (int)atof(sub.c_str()) - 1;
            iss >> sub;
            e = (int)atof(sub.c_str()) - 1;
            chr.push_back(genomic_element(i, s, e));
          }
          get_line(infile, line);
        }
      } else if (line.find("RECOMBINATION RATE") != string::npos) {
        get_line(infile, line);
        parameters.push_back("#RECOMBINATION RATE");
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            parameters.push_back(line);

            // FORMAT: x r (interval end, rec rate in events per bp)

            int x;
            double r;
            istringstream iss(line);
            iss >> sub;
            x = (int)atof(sub.c_str()) - 1;
            iss >> sub;
            r = atof(sub.c_str());
            chr.rec_x.push_back(x);
            chr.rec_r.push_back(r);
          }
          get_line(infile, line);
        }
      } else if (line.find("GENE CONVERSION") != string::npos) {
        get_line(infile, line);
        parameters.push_back("#GENE CONVERSION");
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            parameters.push_back(line);

            // FORMAT: G_f G_l (gene conversion fraction, average stretch length
            // in bp)

            istringstream iss(line);
            iss >> sub;
            chr.G_f = atof(sub.c_str());
            iss >> sub;
            chr.G_l = atof(sub.c_str());
          }
          get_line(infile, line);
        }
      } else if (line.find("GENERATIONS") != string::npos) {
        get_line(infile, line);
        parameters.push_back("#GENERATIONS");
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            parameters.push_back(line);
            istringstream iss(line);
            iss >> sub;
            t_duration = (int)atof(sub.c_str());
            if (iss >> sub) {
              t_start = (int)atof(sub.c_str());
            } else {
              t_start = 1;
            }
          }
          get_line(infile, line);
        }
      } else if (line.find("DEMOGRAPHY AND STRUCTURE") != string::npos) {
        get_line(infile, line);
        parameters.push_back("#DEMOGRAPHY AND STRUCTURE");
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            parameters.push_back(line);

            // FORMAT: t event_type event_parameters

            int time;
            char event_type;
            vector<string> s;
            istringstream iss(line);
            iss >> sub;
            time = (int)atof(sub.c_str());
            iss >> sub;
            event_type = sub.at(0);

            if (event_type == 'T') {
              event_type = 't';
            } else if (event_type == 'R') {
              event_type = 'r';
            }

            while (iss >> sub) {
              if (event_type == 'P' &&
                  sub.at(0) == 'p') {  // Add subpopulation to list, wee need to
                                       // keep them for validate genomic
                                       // populations
                int subpopulation;
                istringstream(sub.substr(1, sub.length())) >> subpopulation;
                subpopulations.push_back(subpopulation);
              }
              s.push_back(sub.c_str());
            }
            E.insert(pair<int, event>(time, event(event_type, s)));
          }
          get_line(infile, line);
        }
      } else if (line.find("HERMAPHRODITES") != string::npos) {
        parameters.push_back("#HERMAPHRODITES");
        get_line(infile, line);
        parameters.push_back(line);
        P.hermaphrodites = translate_hermaphrodite_line(line);
      } else if (line.find("OUTPUT") != string::npos) {
        get_line(infile, line);
        parameters.push_back("#OUTPUT");
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            parameters.push_back(line);

            // FORMAT: t event_type event_paramaters

            int t;
            char c;
            vector<string> s;
            istringstream iss(line);
            iss >> sub;
            t = (int)atof(sub.c_str());
            iss >> sub;
            c = sub.at(0);

            while (iss >> sub) {
              s.push_back(sub.c_str());
            }
            O.insert(pair<int, event>(t, event(c, s)));
          }
          get_line(infile, line);
        }
      } else if (line.find("PREDETERMINED MUTATIONS") != string::npos) {
        get_line(infile, line);
        parameters.push_back("#PREDETERMINED MUTATIONS");
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            parameters.push_back(line);

            // FORMAT: t type x i nAA nAa [P <freq>]

            istringstream iss(line);

            iss >> sub;
            int time = (int)atof(sub.c_str());
            iss >> sub;
            sub.erase(0, 1);
            int t = atoi(sub.c_str());
            iss >> sub;
            int x = (int)atof(sub.c_str()) - 1;
            iss >> sub;
            sub.erase(0, 1);
            int i = atoi(sub.c_str());
            iss >> sub;
            int nAA = (int)atof(sub.c_str());
            iss >> sub;
            int nAa = (int)atof(sub.c_str());

            introduced_mutation M(t, x, i, time, nAA, nAa);

            IM.insert(pair<int, introduced_mutation>(time, M));

            while (iss >> sub) {
              if (sub.find('P') != string::npos) {
                iss >> sub;
                float p = atof(sub.c_str());
                PS.push_back(partial_sweep(t, x, p));
              }
            }
          }
          get_line(infile, line);
        }
      } else if (line.find("SEED") != string::npos) {
        get_line(infile, line);
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            istringstream iss(line);
            iss >> sub;
            seed = atoi(sub.c_str());
          }
          get_line(infile, line);
        }
      } else if (line.find("INITIALIZATION") != string::npos) {
        get_line(infile, line);
        parameters.push_back("#INITIALIZATION");
        while (line.find('#') == string::npos && !infile.eof()) {
          if (line.length() > 0) {
            parameters.push_back(line);

            istringstream iss(line);
            iss >> sub;
            initialize_from_file(P, sub.c_str(), chr);
          }
          get_line(infile, line);
        }
      }
    } else {
      get_line(infile, line);
    }
  }

  // initialize chromosome
  chr.validate(subpopulations);
  chr.initialize_rng();

  // initialize rng

  rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng, (long)seed);

  parameters.push_back("#SEED");
  stringstream ss;
  ss << seed;
  parameters.push_back(ss.str());

  // parameter output

  for (int i = 0; i < P.parameters.size(); i++) {
    cout << parameters[i] << endl;
  }
  // Cleaning variables
  subpopulations.clear();
}

int main(int argc, char* argv[]) {
  // initialize simulation parameters

  if (argc <= 1) {
    cerr << "usage: slim <parameter file>" << endl;
    exit(1);
  }

  char* input_file = argv[1];
  check_input_file(input_file);

  int t_start; /**< Time start of the simulation, retrieved from config file  */
  int t_duration; /**< Time end of the simulation, retrieved from config file */
  chromosome chr;

  population P;
  map<int, subpopulation*>::iterator itP;

  P.parameters.push_back("#INPUT PARAMETER FILE");
  P.parameters.push_back(input_file);

  // demographic and structure events

  multimap<int, event> E;
  multimap<int, event>::iterator itE;

  // output events (time, output)

  multimap<int, event> O;
  multimap<int, event>::iterator itO;

  // user-defined mutations that will be introduced (time, mutation)

  multimap<int, introduced_mutation> IM;
  multimap<int, introduced_mutation>::iterator itIM;

  // tracked mutation-types

  vector<int> TM;

  // mutations undergoing partial sweeps

  vector<partial_sweep> PS;

  initialize(P, input_file, chr, t_start, t_duration, E, O, IM, PS,
             P.parameters);

  cout << t_start << " " << t_duration << endl;

  // evolve over t generations

  for (int g = t_start; g < (t_start + t_duration); g++) {
    // execute demographic and substructure events in this generation

    pair<multimap<int, event>::iterator, multimap<int, event>::iterator>
        rangeE = E.equal_range(g);
    for (itE = rangeE.first; itE != rangeE.second; itE++) {
      P.execute_event(itE->second, g, chr, TM);
    }

    // evolve all subpopulations

    for (itP = P.begin(); itP != P.end(); itP++) {
      P.evolve_subpopulation(itP->first, chr, g);
    }

    // introduce user-defined mutations

    pair<multimap<int, introduced_mutation>::iterator,
         multimap<int, introduced_mutation>::iterator> rangeIM =
        IM.equal_range(g);
    for (itIM = rangeIM.first; itIM != rangeIM.second; itIM++) {
      P.introduce_mutation(itIM->second, chr);
    }

    // execute output events

    pair<multimap<int, event>::iterator, multimap<int, event>::iterator>
        rangeO = O.equal_range(g);

    for (itO = rangeO.first; itO != rangeO.second; itO++) {
      P.execute_event(itO->second, g, chr, TM);
    }

    // track particular mutation-types and set s=0 for partial sweeps when
    // completed

    if (TM.size() > 0 || PS.size() > 0) {
      P.track_mutations(g, TM, PS, chr);
    }

    // swap generations

    P.swap_generations(g, chr);
  }
  return 0;
}
