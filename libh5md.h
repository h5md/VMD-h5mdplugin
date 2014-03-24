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
