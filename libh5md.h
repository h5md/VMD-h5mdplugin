#ifndef __H5MD_H
#define __H5MD_H

struct h5md_file;

// all functions return 0 upon success and != 0 upon failure

// opens the file, creates the internal structure and goes to the first timestep
int h5md_open(struct h5md_file** file, const char *filename);

// close the file and frees the internal structure
int h5md_close(struct h5md_file* file);

// return the current error message 
const char* h5md_error(struct h5md_file* file);

// go to the i'th timestep
int h5md_seek_timestep(struct h5md_file* file, int i);

// reads the next timestep, allocates coords and sets natoms to the number of atoms
int h5md_get_timestep(struct h5md_file*, int* natoms, double **coords);

// reads the next timestep and writes the data to coords iff natoms is the number of atoms in the timestep
int h5md_read_timestep(struct h5md_file*, int natoms, double* coords);

typedef struct {
	int nspacedims;
	int natoms;
	int ntime;
	//file handling
	char *file_name;
	hid_t file_id;
	hid_t dataset_id;
	//bonds
	int *nbonds;
	int *bond_from;
	int *bond_to;
	//positions
	double ***data_xyz;
	int ndims;
}h5mddata_lib;

void open_h5md_read_lib(const char *filename, h5mddata_lib* data);
void read_position_lib(h5mddata_lib *data);
int read_h5md_timestep_lib(h5mddata_lib*data, int natoms, double ***data_xyz, int time_i);
void read_position_of_file_lib(const char *filename, h5mddata_lib* data);

typedef struct{
	hid_t dataset_id_species;
	int *data_species;
}h5mdspecies;

void read_species(h5mddata_lib *data, h5mdspecies* species, int* data_speciess);

#endif
