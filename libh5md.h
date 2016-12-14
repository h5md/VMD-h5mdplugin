#ifndef __H5MD_H
#define __H5MD_H

struct h5md_file;

// all functions return 0 upon success and != 0 upon failure

//hide hdf5 error messages
void h5md_hide_hdf5_error_messages();

//show hdf5 error messages
void h5md_show_hdf5_error_messages();

/*read operations*/

// opens the file, iff it exists and creates the internal structure and goes to the first timestep
int h5md_open(struct h5md_file** _file, const char *filename, int can_open);

// close the file and frees the internal structure
int h5md_close(struct h5md_file* file);

// return the current error message 
const char* h5md_error(struct h5md_file* file);

// go to the i'th timestep
int h5md_seek_timestep(struct h5md_file* file, int i);

// reads the next timestep of all groups, allocates coords and sets natoms to the number of atoms
int h5md_get_timestep(struct h5md_file* file, int* natoms, float **coords);

//reads all box informations of all groups, returns only the box information of the first group (since VMD doesn't support more than one box per file), the returned array has length 6 (3 lengths: A, B, C, and 3 angles alpha, beta, gamma)
int h5md_get_box_information(struct h5md_file* file, float* out_box_information);

//reads species information of all groups, allocates species_infromation_out
int h5md_get_all_species_infromation(struct h5md_file *file, int** species_infromation_out);

//reads mass information of all groups, allocates mass_infromation_out
int h5md_get_all_mass_infromation(struct h5md_file *file, float** mass_infromation_out);

//reads charge information of all groups, allocates charge_infromation_out
int h5md_get_all_charge_infromation(struct h5md_file *file, float** charge_infromation_out);

// reads the next timestep of all groups and writes the data to coords iff natoms is the number of atoms in the timestep
int h5md_read_timestep(struct h5md_file* file, int natoms, float* coords);

//reads the number of atoms
int h5md_get_natoms(struct h5md_file* file, int* natoms);

//set the number of atoms for a dataset
int h5md_set_natoms(struct h5md_file* file, int natoms);

//reads the current time
int h5md_get_current_time(struct h5md_file* file, int* current_time);

//reads the total number of timesteps
int h5md_get_ntime(struct h5md_file* file, int* ntime);

//read time-independent dataset automatically
int h5md_read_timeindependent_dataset_automatically(struct h5md_file* file, char* dataset_name, void** _data_out, H5T_class_t* type_class_out);

//read time-independent integer dataset automatically
int h5md_read_timeindependent_dataset_int(struct h5md_file* file, char* dataset_name, int** _data_out);

//free time-independent dataset automatically
int h5md_free_timeindependent_dataset_automatically(H5T_class_t type_class, void* old_data_out, int length_array_of_strings);

//get length of one-dimensional dataset
int h5md_get_length_of_one_dimensional_dataset(struct h5md_file *file,char *dataset_name, int *length_of_dataset);

/*write operations*/

//creates a h5 file and sets the author (if not overwritten later, the username of the currently logged in user is used) and saves the reference to this bare file in a h5md_file struct, fails if a file with the provided filename already exists
int h5md_create_file(struct h5md_file **_file, const char* filename);

//sets the author's name and email. if name==NULL then the username of the currently logged in user is used, if email_address==NULL, then it is remarked, that no email was provided
int h5md_set_author(struct h5md_file* file, char* name, char* email_address);

//sets the creator program name and version. if name==NULL then the string "libh5md" is set. If version==NULL, then it is remarked, that no version was provided
int h5md_set_creator(struct h5md_file* file, char* name, char* version);

//deletes a file
int h5md_delete_file(char* filename);

//writes a dataset of given datatype
int h5md_write_dataset(struct h5md_file *file, char* absolute_name_of_dataset, hid_t datatype, void* data_in, int rank_in, hsize_t* dims_in);

//appends to an existing dataset or creates it
int h5md_append_dataset(struct h5md_file *file, char* absolute_name_of_dataset, hid_t datatype, void* positions_in, int rank_in, hsize_t* dims_in);

//gets the fill value of a dataset
int get_fill_value(struct h5md_file* file, char* absolute_name_of_dataset, void* filler);

//deletes a dataset
int h5md_delete_dataset(struct h5md_file* file, char* absolute_name_of_dataset);

/* correction for VMD counting timesteps from 1 onwards */
void h5md_set_correction_for_VMD_counting_timesteps(struct h5md_file* file);
int h5md_get_correction_for_VMD_counting_timesteps(struct h5md_file* file);

#endif
