# include "gtest/gtest.h"
# include "translator.h"

TEST(translate_genomic_elements_type_line_test, DISABLED_validLineNoPop) {
	genomicElementTypesLine the_line;
	the_line = translate_genomic_elements_type_line("g2 m3 0.5");
	EXPECT_EQ(2, the_line.id);
	EXPECT_EQ(0, the_line.population);
	EXPECT_EQ(3, the_line.mutations.front());
	EXPECT_EQ(0.5, the_line.fractions.front());
}

TEST(translate_genomic_elements_type_line_test, validLineWithPop) {
	genomicElementTypesLine the_line;
	the_line = translate_genomic_elements_type_line("p1 g2 m3 0.5");
	EXPECT_EQ(2, the_line.id);
	EXPECT_EQ(1, the_line.population);
	EXPECT_EQ(3, the_line.mutations.front());
	EXPECT_EQ(0.5, the_line.fractions.front());
}
