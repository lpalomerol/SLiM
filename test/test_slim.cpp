#include "gtest/gtest.h"
#include "../slim.cpp"

class slim_test : public ::testing::Test {

 protected:
  static void SetUpTestCase() {
    rng = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng, (long)12345);
  }
};

typedef slim_test test_chromosome_validation;
typedef slim_test test_population_set_threshold;
typedef slim_test no_hermaphrodites;
typedef slim_test set_sex_ratio;
typedef slim_test get_reproductive_females;
typedef slim_test update_threshold_fitness;
typedef slim_test equiv_genome_test;

TEST_F(test_chromosome_validation, should_crash_when_empty) {
  chromosome chr;
  vector<int> subpopulations;

  ASSERT_DEATH(chr.validate(subpopulations), "empty chromosome");
}

TEST_F(test_chromosome_validation, should_crash_when_rec_r_empty) {
  chromosome chr;
  vector<int> subpopulations;

  chr.push_back(genomic_element(1, 1, 1));
  ASSERT_DEATH(chr.validate(subpopulations),
               "recombination rate not specified");
}

TEST_F(test_chromosome_validation, sould_crash_when_invalid_mutation_rate) {
  chromosome chr;
  vector<int> subpopulations;
  chr.push_back(genomic_element(1, 1, 1));
  chr.rec_r.push_back(1);
  chr.M = -1.4;
  ASSERT_DEATH(chr.validate(subpopulations), "mutation rate");
}

TEST_F(test_chromosome_validation,
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

TEST_F(test_chromosome_validation,
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

TEST_F(test_chromosome_validation,
       should_pass_when_one_is_custom_other_default) {
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

TEST_F(test_chromosome_validation,
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

TEST_F(test_population_set_threshold, should_fail_with_invalid_population) {
  population pop;
  ASSERT_DEATH(pop.set_threshold(0, 10), "no subpopulation source 0");
}

TEST_F(test_population_set_threshold, should_fail_with_negative_threshold) {
  population pop;
  pop.add_subpopulation(1, 100);
  ASSERT_DEATH(pop.set_threshold(1, -1),
               "migration threshold has to be equal or greater than 1");
}

TEST_F(test_population_set_threshold, should_disable_threshold_when_0) {
  population pop;
  pop.add_subpopulation(1, 100);
  pop.set_threshold(1, 0);
  EXPECT_EQ(0, pop.find(1)->second->T);
}

TEST_F(test_population_set_threshold,
       should_update_threshold_when_lower_than_pop_number) {
  population pop;
  pop.add_subpopulation(1, 100);
  pop.set_threshold(1, 50);
  EXPECT_EQ(50, pop.find(1)->second->T);
}

TEST_F(test_population_set_threshold,
       should_disable_threshold_when_greater_than_pop) {
  population pop;
  pop.add_subpopulation(1, 100);
  pop.set_threshold(1, 200);
  EXPECT_EQ(0, pop.find(1)->second->T);
}

// Selfing definition tests
TEST_F(no_hermaphrodites, should_define_selfing_with_hermaprhrodites) {
  population pop;
  pop.add_subpopulation(1, 1000);
  pop.set_selfing(1, 0.5);
  EXPECT_EQ(pop.find(1)->second->S, 0.5);
}

TEST_F(no_hermaphrodites,
       should_crash_defining_selfing_with_no_hermaprhrodites) {

  population pop;
  pop.hermaphrodites = 0;
  pop.add_subpopulation(1, 1000);
  ASSERT_DEATH(pop.set_selfing(1, 0.5), "selfing fraction");
}

////Adding children tests
TEST_F(no_hermaphrodites, should_add_child_defining_sex) {
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

TEST_F(no_hermaphrodites, should_add_individual_from_p2_to_p1_with_sex_fixed) {
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

TEST_F(no_hermaphrodites, check_draw_individual_by_sex_with_no_hermaprh) {
  population pop;
  chromosome chrom;
  pop.hermaphrodites = 0;
  pop.add_subpopulation(1, 2);
  pop.find(1)->second->parent_is_male(0, true);
  pop.find(1)->second->parent_is_male(1, false);
  pop.find(1)->second->update_fitness(chrom);
  int male_id = pop.find(1)->second->draw_individual_by_sex(true);
  int female_id = pop.find(1)->second->draw_individual_by_sex(false);
  EXPECT_EQ(pop.find(1)->second->parent_is_male(male_id), true);
  EXPECT_EQ(pop.find(1)->second->parent_is_male(female_id), false);
}

TEST_F(no_hermaphrodites, check_draw_individual_by_sex_with_hermaprh) {
  population pop;
  pop.hermaphrodites = 1;
  pop.add_subpopulation(1, 2);
  pop.find(1)->second->parent_is_male(1, true);
  pop.find(1)->second->parent_is_male(0, false);
  int male_id = pop.find(1)->second->draw_individual_by_sex(true);
  int female_id = pop.find(1)->second->draw_individual_by_sex(false);
  EXPECT_EQ(pop.find(1)->second->parent_is_male(male_id), true);
  // With hermaphrodites, both can be "male"
  EXPECT_EQ(pop.find(1)->second->parent_is_male(female_id), true);
}

TEST_F(set_sex_ratio, should_crash_when_sex_ratio_hermaprhodites) {
  population pop;
  pop.hermaphrodites = 1;
  pop.add_subpopulation(1, 100);
  ASSERT_DEATH(pop.set_sex_ratio(1, 0.5), "sex ratio");
}

TEST_F(set_sex_ratio, should_crash_when_sex_ratio_non_existent_population) {
  population pop;
  pop.hermaphrodites = 0;
  pop.add_subpopulation(1, 100);
  pop.set_sex_ratio(1, 0.15);
  EXPECT_EQ(pop.find(1)->second->S_ratio, 0.15);
}

TEST_F(get_reproductive_females, reproductive_females_tests) {

  /*
   * 100 individuals, 50 males, no threshold, ratio 1.0 -> 50 reproductive
   * females
   */
  EXPECT_EQ(subpopulation_sexed::get_reproductive_females(100, 50, 0, 1.0), 50);

  /*
   * 100 individuals, 55 males, no threshold, ratio 1.0 -> 45 reproductive
   * females
   */
  EXPECT_EQ(subpopulation_sexed::get_reproductive_females(100, 55, 0, 1.0), 45);

  /*
   * 100 individuals, 45 males, no threshold, ratio 1.0 -> 45 reproductive
   * females
   */
  EXPECT_EQ(subpopulation_sexed::get_reproductive_females(100, 45, 0, 1.0), 45);

  /*
   * 100 indiv, 50 males, 10 males threshold, ratio 1.0 -> 10 reproductive
   * females
   */
  EXPECT_EQ(subpopulation_sexed::get_reproductive_females(100, 50, 10, 1.0),
            10);

  /*
   * 100 indiv, 50 males, 10 males threshold, ratio 0.5 -> 20  reproductive
   * females
   */
  EXPECT_EQ(subpopulation_sexed::get_reproductive_females(100, 50, 10, 0.5),
            20);

  /*
   * 2 individuals, 1 male, 10 males threshold, ratio 1.0 -> 1 reproductive
   * female
   */
  EXPECT_EQ(subpopulation_sexed::get_reproductive_females(2, 1, 10, 1.0), 1);

  /*
   * 100 individuals, 1 male, 10 males threshold, ratio 0.5 -> 2 reproductive
   * females
   */
  EXPECT_EQ(subpopulation_sexed::get_reproductive_females(100, 1, 10, 0.5), 2);

  /*
   * 100 individuals, 50 males, 10 males threshold, ratio 0.25 -> 40
   * reproductive females
   */
  EXPECT_EQ(subpopulation_sexed::get_reproductive_females(100, 50, 10, 0.25),
            40);

  /*
   * 100 individuals, 50 males, 10 males threshold, ratio 0.125 -> 50
   * reproductive females ...
   *    (by threshold should be 80 but there are not enough females at the
   * group)
   */
  EXPECT_EQ(subpopulation_sexed::get_reproductive_females(100, 50, 10, 0.125),
            50);

  /*
   * 110 individuals, 50 males, 10 males threshold, ratio 0.125 -> 60
   * reproductive females ...
   *    (by threshold should be 80 but there are not enough females at the
   * group)
   */
  EXPECT_EQ(subpopulation_sexed::get_reproductive_females(110, 50, 10, 0.125),
            60);

  /*
   * Without sex ratio (ratio = 0) should be selected all females
   */
  EXPECT_EQ(subpopulation_sexed::get_reproductive_females(100, 20, 10, 0), 80);

  /*
   * Sex ratio can't be greater than 1
   */
  ASSERT_DEATH(subpopulation_sexed::get_reproductive_females(100, 50, 10, 1.5),
               "Invalid sex ratio");
}

TEST_F(update_threshold_fitness,
       hermaphrodite_should_make_selectable_threshold_subset) {
  chromosome chr;
  population pop;
  pop.hermaphrodites = 1;
  pop.add_subpopulation(1, 10);

  pop.find(1)->second->T = 5;

  pop.find(1)->second->update_fitness(chr);

  // Now we will retrieve a lot of males, only 5 can be selected due threshold
  int selected[10];
  memset(selected, 0, sizeof(selected));

  for (int j = 0; j < 100; j++) {
    selected[pop.find(1)->second->draw_individual_by_sex(true)] = 1;
  }
  int found = 0;
  for (int j = 0; j < 10; j++) {
    if (selected[j] == 1) {
      found++;
    }
  }
  EXPECT_EQ(found, 5);
}

TEST_F(update_threshold_fitness,
       sexed_should_make_selectable_threshold_subset) {
  chromosome chr;
  population pop;
  pop.hermaphrodites = 0;
  pop.add_subpopulation(1, 10);

  // Force 5 males and 5 females
  for (int i = 0; i < 10; i++) {
    pop.find(1)->second->parent_is_male(i, (i < 5) ? true : false);
  }

  pop.find(1)->second->T = 3;

  pop.find(1)->second->update_fitness(chr);

  // Now we will retrieve a lot of males, only 5 can be selected due threshold
  int males[10];
  int females[10];
  memset(males, 0, sizeof(males));
  memset(females, 0, sizeof(females));
  for (int j = 0; j < 100; j++) {
    males[pop.find(1)->second->draw_individual_by_sex(true)] = 1;
    females[pop.find(1)->second->draw_individual_by_sex(false)] = 1;
  }
  int found_males = 0;
  int found_females = 0;
  for (int j = 0; j < 10; j++) {
    if (males[j] == 1) {
      found_males++;
    }
    if (females[j] == 1) {
      found_females++;
    }
  }
  EXPECT_EQ(found_males, 3);
  EXPECT_EQ(found_females, 5);
}

TEST_F(update_threshold_fitness,
       sexed_should_make_selectable_threshold_subset_with_ratio) {
  chromosome chr;
  population pop;
  pop.hermaphrodites = 0;
  pop.add_subpopulation(1, 10);

  // Force 5 males and 5 females
  for (int i = 0; i < 10; i++) {
    pop.find(1)->second->parent_is_male(i, (i < 5) ? true : false);
  }

  pop.find(1)->second->T = 1;
  pop.find(1)->second->S_ratio = 0.5;

  pop.find(1)->second->update_fitness(chr);

  // Now we will retrieve a lot of males, only 5 can be selected due threshold
  int males[10];
  int females[10];
  memset(males, 0, sizeof(males));
  memset(females, 0, sizeof(females));
  for (int j = 0; j < 100; j++) {
    males[pop.find(1)->second->draw_individual_by_sex(true)] = 1;
    females[pop.find(1)->second->draw_individual_by_sex(false)] = 1;
  }
  int found_males = 0;
  int found_females = 0;
  for (int j = 0; j < 10; j++) {
    if (males[j] == 1) {
      found_males++;
    }
    if (females[j] == 1) {
      found_females++;
    }
  }
  // Threshold is 1 male and sex_ratio is 2 females per male.
  EXPECT_EQ(found_males, 1);
  EXPECT_EQ(found_females, 2);
}

TEST_F(equiv_genome_test, should_return_the_same_id_when_found){
	vector<int> m;
	vector<double> g;
	m.push_back(1);
	g.push_back(1.0);
	genomic_element_type origin = genomic_element_type(m,g);
	EXPECT_EQ(origin.find_equivalent(1), 1);
}

TEST_F(equiv_genome_test, should_return_the_eixsten_id_when_not_found_and_only_one_ratio_1){
	vector<int> m;
	vector<double> g;
	m.push_back(2);
	g.push_back(1.0);
	genomic_element_type origin = genomic_element_type(m,g);
	EXPECT_EQ(origin.find_equivalent(1), 2);
}

TEST_F(equiv_genome_test, should_throw_error_when_no_equiv_and_more_than_one_candidate){
	vector<int> m;
	vector<double> g;
	m.push_back(2);
	g.push_back(.5);
	m.push_back(3);
	g.push_back(.5);
	genomic_element_type origin = genomic_element_type(m,g);
	ASSERT_DEATH(origin.find_equivalent(1), "Invalid genomic configuration");
}
