#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "libh5md.h"

//global varible
int position_already_read_lib =-1;

//gathers important information about positions dataset
void open_h5md_read_lib(const char *filename,h5mddata_lib *data){
	hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	hid_t dataset_id = H5Dopen(file_id, "/particles/atoms/position/value",H5P_DEFAULT);
	data->file_id=file_id;
	data->dataset_id=dataset_id;
	//get dims
	hid_t dataspace = H5Dget_space(dataset_id);
	unsigned long long int dims_out[3];
	H5Sget_simple_extent_dims(dataspace, dims_out, NULL );
	//save repeatedly needed files in h5mddata_lib struct 
	data->file_name = filename;
	data->ntime = dims_out[0];
	data->natoms = dims_out[1];
	data->nspacedims = dims_out[2];
        data->data_xyz = NULL;
}

//reads all positions into heap and stores triple pointer in data->data_xyz
void read_position_lib(h5mddata_lib *data){
  //read position data to data_xyz_read
  data->data_xyz = malloc(sizeof(double)*data->ntime*data->natoms*data->nspacedims);
  hid_t file_id = H5Fopen(data->file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
  H5Dread(data->dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data->data_xyz);
}

//writes positions of all atoms at time time_i into the array data_xyz
int read_h5md_timestep_lib(h5mddata_lib *data, double *data_xyz, int time_i) {
	if(position_already_read_lib<0){
		read_position_lib(data);
		position_already_read_lib=1;
	}
	/* read the coordinates */
	unsigned int ntime = data->ntime;
	if (time_i >= ntime - 1) {
		//reset library
		position_already_read_lib=-1;
		return -1;
	}
	for (int i = 0; i < data->natoms; i++) {
          memcpy(data_xyz, data->data_xyz + time_i*data->natoms*3, data->natoms*3);
	}
	return 1;
}

void read_position_of_file_lib(const char *filename, h5mddata_lib* data){
	open_h5md_read_lib(filename, data);
	read_position_lib(data);
}

void read_species(h5mddata_lib* data, h5mdspecies* species, int *data_species){
	//load species
	species->data_species=data_species;
	species->dataset_id_species = H5Dopen(data->file_id, "/particles/atoms/species/value",H5P_DEFAULT);
	H5Dread(species->dataset_id_species, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, species->data_species);
	H5Dclose (species->dataset_id_species);
}
