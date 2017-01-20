#ifndef GTEST_LIBH5MD_H_
#define GTEST_LIBH5MD_H_

int test_h5md_get_timestep();
int test_create_file();
int test_write_file_with_datasets();
int test_write_fixed_length_string_dataset();
int read_string_dataset();
int read_real_dataset_float();
int read_real_dataset_string();
int test_delete_dataset(char* absolute_name_of_dataset);
int test_appending_to_dataset_position_style();
int test_appending_to_onedimensional_dataset_float();
int test_appending_to_onedimensional_string_dataset();
int test_appending_to_dataset_position_style_grandcanonical();
int test_get_fill_value_float();

#endif  // GTEST_LIBH5MD_H_
