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
int h5md_open(struct h5md_file** _file, const char *filename);

// close the file and frees the internal structure
int h5md_close(struct h5md_file* file);

// return the current error message 
const char* h5md_error(struct h5md_file* file);

// go to the i'th timestep
int h5md_seek_timestep(struct h5md_file* file, int i);

// reads the next timestep, allocates coords and sets natoms to the number of atoms
int h5md_get_timestep(struct h5md_file* file, int* natoms, float **coords);

// reads the next timestep and writes the data to coords iff natoms is the number of atoms in the timestep
int h5md_read_timestep(struct h5md_file* file, int natoms, float* coords);

//reads the number of atoms
int h5md_get_natoms(struct h5md_file* file, int* natoms);

//reads the current time
int h5md_get_current_time(struct h5md_file* file, int* current_time);

//reads the total number of timesteps
int h5md_get_ntime(struct h5md_file* file, int* ntime);

//read time-independent dataset automatically
int h5md_read_timeindependent_dataset_automatically(struct h5md_file* file, char* dataset_name, void** _data_out, H5T_class_t* type_class_out);

int h5md_read_timeindependent_dataset_int(struct h5md_file* file, char* dataset_name, int** _data_out);

//free time-independent dataset automatically
int h5md_free_timeindependent_dataset_automatically(H5T_class_t type_class, void* old_data_out);

//get length of one-dimensional dataset
int h5md_get_length_of_one_dimensional_dataset(struct h5md_file *file,char *dataset_name, int *length_of_dataset);



/*write operations*/
int h5md_create_file(struct h5md_file **_file, char* filename);
int h5md_delete_file(char* filename);
int h5md_write_dataset(struct h5md_file *file, char* absolute_name_of_dataset, hid_t datatype, void* data_in, int rank_in, hsize_t* dims_in);
int h5md_delete_dataset(struct h5md_file* file, char* absolute_name_of_dataset);
int h5md_set_author(struct h5md_file* file, char* name);



#endif
