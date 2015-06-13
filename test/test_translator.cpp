# include "gtest/gtest.h"
# include "translator.h"

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


TEST(translate_hermaphrodite_line_test, valid_line_0){
	int hermaph;
	hermaph = translate_hermaphrodite_line("0");
	EXPECT_EQ(0, hermaph);
}


TEST(translate_hermaphrodite_line_test, valid_line_1_comment){
	int hermaph;
	hermaph = translate_hermaphrodite_line("1 / I am a comment");
	EXPECT_EQ(1, hermaph);
}

