#include "gtest/gtest.h"
#include "../slim.cpp"

TEST(test_chromosome_validation, should_crash_when_empty) {
  chromosome chr;
  vector<int> subpopulations;

  ASSERT_DEATH(chr.validate(subpopulations), "empty chromosome");
}

TEST(test_chromosome_validation, should_crash_when_rec_r_empty) {
  chromosome chr;
  vector<int> subpopulations;

  chr.push_back(genomic_element(1, 1, 1));
  ASSERT_DEATH(chr.validate(subpopulations),
               "recombination rate not specified");
}

TEST(test_chromosome_validation, sould_crash_when_invalid_mutation_rate) {
  chromosome chr;
  vector<int> subpopulations;
  chr.push_back(genomic_element(1, 1, 1));
  chr.rec_r.push_back(1);
  chr.M = -1.4;
  ASSERT_DEATH(chr.validate(subpopulations), "mutation rate");
}

TEST(test_chromosome_validation,
     should_pass_when_all_genomic_elements_are_default) {
  chromosome chr;
  vector<int> subpopulations;
  chr.push_back(genomic_element(1, 1, 1));
  chr.rec_r.push_back(1);
  chr.M = 1.4;

  genomic_element_types_line line1;
  line1.genomic_id = 1;
  line1.pop_id = 0;
  line1.mutations.push_back(1);
  line1.fractions.push_back(1);

  genomic_element_types_line line2;
  line2.genomic_id = 2;
  line2.pop_id = 0;
  line2.mutations.push_back(1);
  line2.fractions.push_back(1);

  subpopulations.push_back(1);
  subpopulations.push_back(2);

  chr.genomic_element_types.insert(pair<pair<int, int>, genomic_element_type>(
      make_pair(line1.pop_id, line1.genomic_id),
      genomic_element_type(line1.mutations, line1.fractions)));
  chr.genomic_element_types.insert(pair<pair<int, int>, genomic_element_type>(
      make_pair(line2.pop_id, line2.genomic_id),
      genomic_element_type(line2.mutations, line2.fractions)));
  ASSERT_DEATH(chr.validate(subpopulations), "mutation type");
}

TEST(test_chromosome_validation,
     should_pass_when_all_genomic_elements_are_custom) {
  chromosome chr;
  vector<int> subpopulations;
  chr.push_back(genomic_element(1, 1, 1));
  chr.rec_r.push_back(1);
  chr.M = 1.4;

  genomic_element_types_line line1;
  line1.genomic_id = 1;
  line1.pop_id = 1;
  line1.mutations.push_back(1);
  line1.fractions.push_back(1);

  genomic_element_types_line line2;
  line2.genomic_id = 2;
  line2.pop_id = 1;
  line2.mutations.push_back(1);
  line2.fractions.push_back(1);

  genomic_element_types_line line3;
  line3.genomic_id = 1;
  line3.pop_id = 2;
  line3.mutations.push_back(1);
  line3.fractions.push_back(1);

  genomic_element_types_line line4;
  line4.genomic_id = 2;
  line4.pop_id = 2;
  line4.mutations.push_back(1);
  line4.fractions.push_back(1);

  subpopulations.push_back(1);
  subpopulations.push_back(2);

  chr.genomic_element_types.insert(pair<pair<int, int>, genomic_element_type>(
      make_pair(line1.pop_id, line1.genomic_id),
      genomic_element_type(line1.mutations, line1.fractions)));
  chr.genomic_element_types.insert(pair<pair<int, int>, genomic_element_type>(
      make_pair(line2.pop_id, line2.genomic_id),
      genomic_element_type(line2.mutations, line2.fractions)));

  chr.genomic_element_types.insert(pair<pair<int, int>, genomic_element_type>(
      make_pair(line3.pop_id, line3.genomic_id),
      genomic_element_type(line3.mutations, line3.fractions)));
  chr.genomic_element_types.insert(pair<pair<int, int>, genomic_element_type>(
      make_pair(line4.pop_id, line4.genomic_id),
      genomic_element_type(line4.mutations, line4.fractions)));

  ASSERT_DEATH(chr.validate(subpopulations), "mutation type");
}

TEST(test_chromosome_validation, should_pass_when_one_is_custom_other_default) {
  chromosome chr;
  vector<int> subpopulations;
  chr.push_back(genomic_element(1, 1, 1));
  chr.rec_r.push_back(1);
  chr.M = 1.4;

  genomic_element_types_line line1;
  line1.genomic_id = 1;
  line1.pop_id = 1;
  line1.mutations.push_back(1);
  line1.fractions.push_back(1);

  genomic_element_types_line line2;
  line2.genomic_id = 2;
  line2.pop_id = 1;
  line2.mutations.push_back(1);
  line2.fractions.push_back(1);

  genomic_element_types_line line3;
  line3.genomic_id = 1;
  line3.pop_id = 0;
  line3.mutations.push_back(1);
  line3.fractions.push_back(1);

  subpopulations.push_back(1);
  subpopulations.push_back(2);

  chr.genomic_element_types.insert(pair<pair<int, int>, genomic_element_type>(
      make_pair(line1.pop_id, line1.genomic_id),
      genomic_element_type(line1.mutations, line1.fractions)));
  chr.genomic_element_types.insert(pair<pair<int, int>, genomic_element_type>(
      make_pair(line2.pop_id, line2.genomic_id),
      genomic_element_type(line2.mutations, line2.fractions)));

  chr.genomic_element_types.insert(pair<pair<int, int>, genomic_element_type>(
      make_pair(line3.pop_id, line3.genomic_id),
      genomic_element_type(line3.mutations, line3.fractions)));

  ASSERT_DEATH(chr.validate(subpopulations), "mutation type");
}

TEST(test_chromosome_validation,
     should_fail_pass_when_one_is_custom_other_not_exists) {
  chromosome chr;
  vector<int> subpopulations;
  chr.push_back(genomic_element(1, 1, 1));
  chr.rec_r.push_back(1);
  chr.M = 1.4;

  genomic_element_types_line line1;
  line1.genomic_id = 1;
  line1.pop_id = 1;
  line1.mutations.push_back(1);
  line1.fractions.push_back(1);

  genomic_element_types_line line2;
  line2.genomic_id = 2;
  line2.pop_id = 1;
  line2.mutations.push_back(1);
  line2.fractions.push_back(1);

  subpopulations.push_back(1);
  subpopulations.push_back(2);

  chr.genomic_element_types.insert(pair<pair<int, int>, genomic_element_type>(
      make_pair(line1.pop_id, line1.genomic_id),
      genomic_element_type(line1.mutations, line1.fractions)));
  chr.genomic_element_types.insert(pair<pair<int, int>, genomic_element_type>(
      make_pair(line2.pop_id, line2.genomic_id),
      genomic_element_type(line2.mutations, line2.fractions)));

  ASSERT_DEATH(chr.validate(subpopulations), "2,1) not defined");
}

TEST(test_population_set_threshold, should_fail_with_invalid_population) {
  population pop;
  ASSERT_DEATH(pop.set_threshold(0, 10), "no subpopulation source 0");
}

TEST(test_population_set_threshold, should_fail_with_negative_threshold) {
  population pop;
  pop.add_subpopulation(1, 100);
  ASSERT_DEATH(pop.set_threshold(1, -1),
               "migration threshold has to be equal or greater than 1");
}

TEST(test_population_set_threshold, should_disable_threshold_when_0) {
  population pop;
  pop.add_subpopulation(1, 100);
  pop.set_threshold(1, 0);
  EXPECT_EQ(0, pop.find(1)->second->T);
}

TEST(test_population_set_threshold,
     should_update_threshold_when_lower_than_pop_number) {
  population pop;
  pop.add_subpopulation(1, 100);
  pop.set_threshold(1, 50);
  EXPECT_EQ(50, pop.find(1)->second->T);
}

TEST(test_population_set_threshold,
     should_disable_threshold_when_greater_than_pop) {
  population pop;
  pop.add_subpopulation(1, 100);
  pop.set_threshold(1, 200);
  EXPECT_EQ(0, pop.find(1)->second->T);
}

TEST(test_select_threshold, should_define_threshold_genome_when_defined) {
  chromosome chr;
  population pop;
  pop.add_subpopulation(1, 1000);
  pop.find(1)->second->update_fitness(chr);
  pop.find(1)->second->T = 500;

  rng = gsl_rng_alloc(gsl_rng_taus2);

  pop.find(1)->second->select_threshold();
  EXPECT_EQ(500, pop.find(1)->second->G_parent_threshold.size());
}

// TEST(test_select_threshold, should_fail_when_threshold_is_0);

// Selfing definition tests
TEST(no_hermaphrodites, should_define_selfing_with_hermaprhrodites) {
  population pop;
  pop.add_subpopulation(1, 1000);
  pop.set_selfing(1, 0.5);
  EXPECT_EQ(pop.find(1)->second->S, 0.5);
}

TEST(no_hermaphrodites, should_crash_defining_selfing_with_no_hermaprhrodites) {
  population pop;
  pop.hermaphrodites = 0;
  pop.add_subpopulation(1, 1000);
  ASSERT_DEATH(pop.set_selfing(1, 0.5), "selfing fraction");
}

////Adding children tests
TEST(no_hermaphrodites, should_add_child_defining_sex) {
  population pop;
  pop.hermaphrodites = 0;
  pop.add_subpopulation(1, 100);
  // Add a male
  pop.find(1)->second->parent_is_male(1, true);
  EXPECT_EQ(pop.find(1)->second->parent_is_male(1), true);
  // Override a female
  pop.find(1)->second->parent_is_male(1, false);
  EXPECT_EQ(pop.find(1)->second->parent_is_male(1), false);

  // Override a male again
  pop.find(1)->second->parent_is_male(1, true);
  EXPECT_EQ(pop.find(1)->second->parent_is_male(1), true);
}

TEST(no_hermaphrodites, should_add_individual_from_p2_to_p1_with_sex_fixed) {
  population pop;
  pop.hermaphrodites = 0;
  pop.add_subpopulation(1, 5);
  // Forcing all males for avoiding randomization
  for (int i = 0; i < 5; i++) {
    pop.find(1)->second->parent_is_male(i, true);
  }
  pop.add_subpopulation(2, 1, 5);
  // After defining a second subpopulation, the sex of the individuals are the
  // same than pop1
  for (int i = 0; i < 5; i++) {
    EXPECT_EQ(pop.find(2)->second->parent_is_male(i), true);
  }
}

TEST(no_hermaphrodites, check_draw_individual_by_sex_with_no_hermaprh) {
  population pop;
  pop.hermaphrodites = 0;
  pop.add_subpopulation(1, 2);
  pop.find(1)->second->parent_is_male(0, true);
  pop.find(1)->second->parent_is_male(1, false);
  bool male_id = pop.find(1)->second->draw_individual_by_sex(true);
  bool female_id = pop.find(1)->second->draw_individual_by_sex(false);
  EXPECT_EQ(pop.find(1)->second->parent_is_male(male_id), true);
  EXPECT_EQ(pop.find(1)->second->parent_is_male(female_id), false);
}

TEST(no_hermaphrodites, check_draw_individual_by_sex_with_hermaprh) {
  population pop;
  pop.hermaphrodites = 1;
  pop.add_subpopulation(1, 2);
  pop.find(1)->second->parent_is_male(0, true);
  pop.find(1)->second->parent_is_male(1, false);
  bool male_id = pop.find(1)->second->draw_individual_by_sex(true);
  bool female_id = pop.find(1)->second->draw_individual_by_sex(false);
  EXPECT_EQ(pop.find(1)->second->parent_is_male(male_id), true);
  // With hermaphrodites, both can be "male"
  EXPECT_EQ(pop.find(1)->second->parent_is_male(female_id), true);
}
