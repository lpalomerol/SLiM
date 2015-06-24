#include "gtest/gtest.h"
#include "translator.h"

TEST(translate_genomic_elements_type_line_test, valid_line_no_pop) {
  genomic_element_types_line the_line;
  the_line = translate_genomic_elements_type_line("g2 m3 0.5");
  EXPECT_EQ(2, the_line.genomic_id);
  EXPECT_EQ(0, the_line.pop_id);
  EXPECT_EQ(3, the_line.mutations.front());
  EXPECT_EQ(0.5, the_line.fractions.front());
}

TEST(translate_genomic_elements_type_line_test, valid_line_with_pop) {
  genomic_element_types_line the_line;
  the_line = translate_genomic_elements_type_line("p1 g2 m3 0.5");
  EXPECT_EQ(2, the_line.genomic_id);
  EXPECT_EQ(1, the_line.pop_id);
  EXPECT_EQ(3, the_line.mutations.front());
  EXPECT_EQ(0.5, the_line.fractions.front());
}

TEST(translate_hermaphrodite_line_test, valid_line_0) {
  int hermaph;
  hermaph = translate_hermaphrodite_line("0");
  EXPECT_EQ(0, hermaph);
}

TEST(translate_hermaphrodite_line_test, valid_line_1_comment) {
  int hermaph;
  hermaph = translate_hermaphrodite_line("1 / I am a comment");
  EXPECT_EQ(1, hermaph);
}


TEST(translate_R_output_line, check_valid_sample_p1_10) {
  output_line translated;
  translated = translate_R_output("500 R p1 10");
  EXPECT_EQ(translated.time, 500);
  EXPECT_EQ(translated.event, 'R');
  EXPECT_EQ(translated.params.size(), 4);
  EXPECT_EQ(translated.params[0], "p1");
  EXPECT_EQ(translated.params[1], "10");
  EXPECT_EQ(translated.params[2], "0");
  EXPECT_EQ(translated.params[3], "DF");
}

TEST(translate_R_output_line, check_valid_sample_p1_10_MS) {
  output_line translated;
  translated = translate_R_output("500 R p1 10 MS");
  EXPECT_EQ(translated.time, 500);
  EXPECT_EQ(translated.event, 'R');
  EXPECT_EQ(translated.params.size(), 4);
  EXPECT_EQ(translated.params[0], "p1");
  EXPECT_EQ(translated.params[1], "10");
  EXPECT_EQ(translated.params[2], "0");
  EXPECT_EQ(translated.params[3], "MS");
}

TEST(translate_R_output_line, check_valid_sample_p1_10_123) {
  output_line translated;
  translated = translate_R_output("500 R p1 10 123");
  EXPECT_EQ(translated.time, 500);
  EXPECT_EQ(translated.event, 'R');
  EXPECT_EQ(translated.params.size(), 4);
  EXPECT_EQ(translated.params[0], "p1");
  EXPECT_EQ(translated.params[1], "10");
  EXPECT_EQ(translated.params[2], "123");
  EXPECT_EQ(translated.params[3], "DF");
}

TEST(translate_R_output_line, check_valid_sample_p1_10_MS_123) {
  output_line translated;
  translated = translate_R_output("500 R p1 10 MS 123");
  EXPECT_EQ(translated.time, 500);
  EXPECT_EQ(translated.event, 'R');
  EXPECT_EQ(translated.params.size(), 4);
  EXPECT_EQ(translated.params[0], "p1");
  EXPECT_EQ(translated.params[1], "10");
  EXPECT_EQ(translated.params[2], "123");
  EXPECT_EQ(translated.params[3], "MS");
}

TEST(translate_R_output_line, check_valid_sample_p1_10_p2_20) {
  output_line translated;
  translated = translate_R_output("500 R p1 10 p2 20");
  EXPECT_EQ(translated.time, 500);
  EXPECT_EQ(translated.event, 'R');
  EXPECT_EQ(translated.params.size(), 6);
  EXPECT_EQ(translated.params[0], "p1");
  EXPECT_EQ(translated.params[1], "10");
  EXPECT_EQ(translated.params[2], "p2");
  EXPECT_EQ(translated.params[3], "20");
  EXPECT_EQ(translated.params[4], "0");
  EXPECT_EQ(translated.params[5], "DF");

}

TEST(translate_R_output_line, check_valid_sample_p1_10_p2_20_MS) {
  output_line translated;
  translated = translate_R_output("500 R p1 10 p2 20 MS");
  EXPECT_EQ(translated.time, 500);
  EXPECT_EQ(translated.event, 'R');
  EXPECT_EQ(translated.params.size(), 6);
  EXPECT_EQ(translated.params[0], "p1");
  EXPECT_EQ(translated.params[1], "10");
  EXPECT_EQ(translated.params[2], "p2");
  EXPECT_EQ(translated.params[3], "20");
  EXPECT_EQ(translated.params[4], "0");
  EXPECT_EQ(translated.params[5], "MS");
}

TEST(translate_R_output_line, check_valid_sample_p1_10_p2_20_MS_1234) {
  output_line translated;
  translated = translate_R_output("500 R p1 10 p2 20 MS 1234");
  EXPECT_EQ(translated.time, 500);
  EXPECT_EQ(translated.event, 'R');
  EXPECT_EQ(translated.params.size(), 6);
  EXPECT_EQ(translated.params[0], "p1");
  EXPECT_EQ(translated.params[1], "10");
  EXPECT_EQ(translated.params[2], "p2");
  EXPECT_EQ(translated.params[3], "20");
  EXPECT_EQ(translated.params[4], "1234");
  EXPECT_EQ(translated.params[5], "MS");
}
TEST(translate_R_output_line, check_valid_sample_p1_10_p2_20_1234) {
  output_line translated;
  translated = translate_R_output("500 R p1 10 p2 20 1234");
  EXPECT_EQ(translated.time, 500);
  EXPECT_EQ(translated.event, 'R');
  EXPECT_EQ(translated.params.size(), 6);
  EXPECT_EQ(translated.params[0], "p1");
  EXPECT_EQ(translated.params[1], "10");
  EXPECT_EQ(translated.params[2], "p2");
  EXPECT_EQ(translated.params[3], "20");
  EXPECT_EQ(translated.params[4], "1234");
  EXPECT_EQ(translated.params[5], "DF");
}
