#include <limits.h>
#include "libh5md_test.hpp"
#include "gtest.h"

TEST(h5mdWriteTest, test_create_file) {
	// This test is named "createFile", and belongs to the "h5mdWriteTest"
	// test case.
	EXPECT_EQ(0, test_create_file()); //EXPECT_EQ checks equivalence, EXPECT_GT checks for one value being greater
}

TEST(h5mdWriteTest, test_write_file_with_datasets) {
	EXPECT_EQ(0, test_write_file_with_datasets());
}

//TEST(h5mdWriteTest, writeFixedLengthString) {
//	EXPECT_EQ(0, test_write_fixed_length_string_dataset());
//}

TEST(h5mdWriteTest, read_string_dataset) {
	EXPECT_EQ(0, read_string_dataset()); 
}

TEST(h5mdWriteTest, test_appending_to_dataset_position_style) {
	EXPECT_EQ(0, test_appending_to_dataset_position_style());
}

TEST(h5mdWriteTest, read_real_dataset_float) {
	EXPECT_EQ(0, read_real_dataset_float());
}


TEST(h5mdWriteTest, read_real_dataset_string) {
	EXPECT_EQ(0, read_real_dataset_string());
}


TEST(h5mdWriteTest, test_delete_dataset) {
	EXPECT_EQ(0, test_delete_dataset("/dataset/go/jjjj"));
}

TEST(h5mdWriteTest, test_appending_to_onedimensional_dataset_float) {
	EXPECT_EQ(0, test_appending_to_onedimensional_dataset_float());
}

TEST(h5mdWriteTest, test_appending_to_onedimensional_string_dataset) {
	EXPECT_EQ(0, test_appending_to_onedimensional_string_dataset());
}

TEST(h5mdWriteTest, test_appending_to_dataset_position_style_grandcanonical) {
	EXPECT_EQ(0, test_appending_to_dataset_position_style_grandcanonical());
}

TEST(h5mdWriteTest, test_get_fill_value_float) {
	EXPECT_EQ(-1, test_get_fill_value_float());
}

TEST(h5mdReadTest, test_h5md_get_timestep){
	EXPECT_EQ(0, test_h5md_get_timestep());
}



