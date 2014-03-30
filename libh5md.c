#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "libh5md.h"

#define TRUE	1
#define FALSE	0

/*
Return value 0 means everything is fine
*/

//file handling
struct h5md_file{
	hid_t file_id;
	hid_t pos_dataset_id;
	hid_t datatype;
	H5T_class_t t_class;
	H5T_order_t order;
	size_t size_datatype;
	hid_t dataspace;
	int rank_dataset;
	int nspacedims;
	int natoms;
	int ntime;
	char *last_error_message;
	int current_time;	//for h5md_seek_timestep()
};

// opens the file, creates the internal structure and goes to the first timestep
// you have to use double pointers in order to be able to change a pointer in a foreing function
int h5md_open(struct h5md_file** _file, const char *filename){
	struct h5md_file *file = malloc(sizeof(struct h5md_file));

	file->file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT); // H5F_ACC_RDONLY <-> read only
	file->pos_dataset_id = H5Dopen2(file->file_id, "/particles/atoms/position/value",H5P_DEFAULT);
	
	/*
	* Get datatype and dataspace handles and then query
	* dataset class, order, size, rank and dimensions.
	*/
	file->datatype  = H5Dget_type(file->pos_dataset_id);     /* datatype handle */
	file->t_class     = H5Tget_class(file->datatype);
	file->order     = H5Tget_order(file->datatype);
	file->size_datatype  = H5Tget_size(file->datatype);
	file->dataspace = H5Dget_space(file->pos_dataset_id);    /* dataspace handle */
	file->rank_dataset      = H5Sget_simple_extent_ndims(file->dataspace);
	hsize_t dims_out[file->rank_dataset];
	H5Sget_simple_extent_dims(file->dataspace, dims_out, NULL);
	
	file->ntime = dims_out[0];
	file->natoms = dims_out[1];
	file->nspacedims = dims_out[2];
	
	file->current_time=0;	//set current time to 0
	
	*_file = file;
	if(file->pos_dataset_id <0){
		return -1;
	}else{
		return 0;
	}	
}


// close the file and frees the internal structure
int h5md_close(struct h5md_file* file){
	if(file!=NULL){
		H5Fclose(file->file_id);
		free(file);	//free h5md_file struct
		return 0;	
	}
	else{
		return -1;	
	}

}

// return the current error message 
const char* h5md_error(struct h5md_file* file){
	return file->last_error_message;
}

// get number of timesteps
int h5md_get_ntime(struct h5md_file* file,int* ntime){
	(*ntime)=file->ntime;
	return 0;
}

// get number of atoms iff this number is constant during time
int h5md_get_natoms(struct h5md_file* file, int* natoms){
	*natoms=file->natoms;
	return 0;
}

//get current time
int h5md_get_current_time(struct h5md_file* file, int* current_time){
	(*current_time)=file->current_time;
	return 0;
}


// go to the i'th timestep
int h5md_seek_timestep(struct h5md_file* file, int i){
	int ntime;
	h5md_get_ntime(file,&ntime);
	if(i<ntime){
		file->current_time=i;
		return 0;
	}else{
		return -1;
	}
}



// reads the next timestep, allocates coords and sets natoms to the number of atoms
int h5md_get_timestep(struct h5md_file* file, int* natoms, double **coords){ //TODO rewrite should not store the data in h5md_file
	/* 
	* Define dataset dataspace in file.
	*/
	hid_t dataspace_id=H5Dget_space(file->pos_dataset_id);
	/* 
	* Define hyperslab in the dataset. 
	*/
	hsize_t dataset_slab_offset[file->rank_dataset];
	hsize_t dataset_slab_count[file->rank_dataset];
	dataset_slab_offset[0] = file->current_time;
	dataset_slab_offset[1] = 0;
	dataset_slab_offset[2]=0;
	dataset_slab_count[0] = 1;
	dataset_slab_count[1] = file->natoms;
	dataset_slab_count[2] = file->nspacedims;
	H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, dataset_slab_offset, NULL, dataset_slab_count, NULL);
	
	/*
	* Define memory dataspace.
	*/
	int rank=1;
	int dimsm[rank];
	dimsm[0]=file->natoms* file->nspacedims;

	hid_t memspace_id = H5Screate_simple(rank,dimsm,NULL);

	/* 
	* Define memory hyperslab. 
	*/
	//hsize_t offset_out[1]={0};
	//hsize_t count_out[1]={file->natoms*file->nspacedims};
	//herr_t status = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);

	/*
	* Read data from hyperslab in the file into the hyperslab in 
	* memory and display.
	*/
	float* data_out= malloc(sizeof(float)*file->natoms * file->nspacedims); //data in memory have the order data_out[3*atom_nr+coord_i]
	hid_t wanted_memory_datatype = H5T_NATIVE_FLOAT;
	H5Dread (file->pos_dataset_id, wanted_memory_datatype, memspace_id, dataspace_id, H5P_DEFAULT, data_out);
	*coords=data_out;

	
	int current_time;
	h5md_get_current_time(file,&current_time);
	int status_seek= h5md_seek_timestep(file, current_time+1); //modify timestep
	
	
	natoms=&(file->natoms); //set natoms to the actual number of atoms
	
	
	H5Sclose(memspace_id); //close resources
	
	if(status_seek==0){
		return 0;
	}else{
		return -1;
	}
}


// reads the next timestep and writes the data to coords iff natoms is the number of atoms in the timestep
int h5md_read_timestep(struct h5md_file* file, int natoms, double* coords){
	h5md_get_timestep(file,&natoms,&coords);

}


//read timeindependent dataset automatically
int h5md_read_timeindependent_dataset_automatically(struct h5md_file* file, char* dataset_name, void** _data_out, H5T_class_t* type_class_out){
	hid_t dataset_id = H5Dopen2(file->file_id, dataset_name,H5P_DEFAULT);	
	/*
	* Get datatype and dataspace handles and then query
	* dataset class, order, size, rank and dimensions.
	*/
	hid_t datatype  = H5Dget_type(dataset_id);     /* datatype handle */
	*type_class_out     = H5Tget_class(datatype);
	//H5T_order_t order     = H5Tget_order(datatype);
	size_t size_datatype  = H5Tget_size(datatype);
	hid_t dataspace_id = H5Dget_space(dataset_id);    /* dataspace handle */
	int rank_dataset      = H5Sget_simple_extent_ndims(dataspace_id);
	unsigned long long int dims_dataset[rank_dataset];
	H5Sget_simple_extent_dims(dataspace_id, dims_dataset, NULL);

	if(dataset_id<0){
		printf("Dataset could not be opened.\n");
		return -1;
	}else{
		switch (*type_class_out) {
		case H5T_INTEGER:{
			//determine needed size
			int needed_size=dims_dataset[0];
			int len_dims_dataset=sizeof(dims_dataset)/sizeof(dims_dataset[0]);
			for(int i=1; i<len_dims_dataset; i++){
				needed_size*=dims_dataset[i];
			}
			int* data_out=(int*) malloc(sizeof(size_datatype)*needed_size);
			H5Dread(dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out);
			H5Dclose (dataset_id);
			*(_data_out)=data_out;
			return 0;
			}
		break;
		case H5T_FLOAT:{
			//determine needed size
			int needed_size=dims_dataset[0];
			int len_dims_dataset=sizeof(dims_dataset)/sizeof(dims_dataset[0]);
			for(int i=1; i<len_dims_dataset; i++){
				needed_size*=dims_dataset[i];
			}
			float* data_out=(float*) malloc(sizeof(size_datatype)*needed_size);
			H5Dread(dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out);
			H5Dclose (dataset_id);
			*(_data_out)=data_out;
			return 0;
			}
		break;
		case H5T_STRING:{
				//load content of name
				size_datatype++;// Make room for null terminator with which strings are terminated in C
				/*
				* Get dataspace and allocate memory for read buffer.  This is a
				* two dimensional dataset so the dynamic allocation must be done
				* in steps for fixed-length strings, for variable-length strings HDF5 does some stuff
				*/
				H5Sget_simple_extent_dims(dataspace_id, dims_dataset, NULL);
				// Allocate array of pointers to rows.
				char **data_out = (char **) malloc (dims_dataset[0] * sizeof (char *));		

				if(H5Tis_variable_str(datatype)==0){
					//string length is variable
					// Allocate space for data.
					data_out[0] = (char *) malloc (dims_dataset[0] * size_datatype);
					// Set the rest of the pointers to rows to the correct addresses.
					for (int i=1; i<dims_dataset[0]; i++)
						data_out[i] = data_out[0] + i * size_datatype;
					//Create the memory datatype.
					H5T_class_t type_class_out = H5Tget_class(datatype);
					hid_t memtype = H5Tcopy (H5T_C_S1);
					H5Tset_size (memtype, size_datatype);
					H5Dread (dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out[0]);
					H5Tclose (memtype); //close memtype
				}else{
					//string length is fixed
					//Create the memory datatype.
					H5T_class_t type_class_out = H5Tget_class(datatype);
					hid_t memtype = H5Tcopy (H5T_C_S1);
					H5Tset_size (memtype, H5T_VARIABLE);
					H5Dread (dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out);
					H5Tclose (memtype); //close memtype
				}

				*(_data_out)=data_out;
				return 0;
			}
		break;

		default:{
				printf("Dataset contains datatype that is not H5T_INTEGER, H5T_FLOAT or H5T_STRING. Not implemented case.\n");
			}
		break;
		}
	}


	//close resources
	H5Dclose (dataset_id);
	H5Sclose (dataspace_id);
	H5Tclose (datatype);


}

int h5md_free_timeindependent_dataset_automatically(H5T_class_t type_class, void* old_data_out){
	switch (type_class) {
	case H5T_INTEGER:
		old_data_out=(int*) old_data_out;
		free(old_data_out);
	break;
	case H5T_FLOAT:
		free(old_data_out);
	break;
	case H5T_STRING:
		old_data_out=(char *) old_data_out; //TODO not correct free for an array of strings which should have void **old_data_out -> missing free(*old_data_out);
		free(old_data_out);
	break;
	}

	return 0;

}

int h5md_get_length_of_one_dimensional_dataset(struct h5md_file *file,char *dataset_name, int *length_of_dataset){
	hid_t dataset_id = H5Dopen2(file->file_id, dataset_name,H5P_DEFAULT);
	hid_t dataspace_id = H5Dget_space(dataset_id);    /* dataspace handle */
	int rank_dataset      = H5Sget_simple_extent_ndims(dataspace_id);
	unsigned long long int dims_dataset[rank_dataset];
	H5Sget_simple_extent_dims(dataspace_id, dims_dataset, NULL);
	if(rank_dataset==1){
		
		*length_of_dataset=dims_dataset[0];
		return 0;
	}else{
		printf("Dataset is not one dimensional.\n");	
		return -1;	
	}
}


void h5md_hide_hdf5_error_messages(){
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
}

void h5md_show_hdf5_error_messages(){
	H5Eset_auto(H5E_DEFAULT, H5Eprint, NULL);
}






/* write operations */
//creates a h5md_file if it does not exist yet

int h5md_create_file(struct h5md_file **_file, char* filename){
	struct h5md_file *file = (struct h5md_file*) malloc(sizeof(struct h5md_file));
	/* Create a new file using default properties. */
	hid_t file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT); //H5F_ACC_EXCL <-> fail if file already exists, alternative H5F_ACC_TRUNC <-> overwrite file
	file->file_id=file_id;
	*(_file)=file;
	if(file_id<0){
		printf("ERROR: A file with this filename already exists.\n");
		return -1;	
	}else{
		return 0;
	}
}

int h5md_delete_file(char* filename){
	int status = remove(filename);
	if (status == 0)
		return 0;
	else
		return -1;
}

int h5md_write_dataset(struct h5md_file *file, char* absolute_name_of_dataset, hid_t datatype, void* data_in, int rank_in, hsize_t* dims_in){
	hid_t dataspace_id=H5Screate_simple(rank_in, dims_in, NULL);
	hid_t link_crt_plist = H5Pcreate(H5P_LINK_CREATE);
	H5Pset_create_intermediate_group(link_crt_plist, TRUE);	// Set flag for intermediate group creation
	hid_t dataset_id = H5Dcreate(file->file_id, absolute_name_of_dataset, datatype, dataspace_id, link_crt_plist, H5P_DEFAULT, H5P_DEFAULT ); //create dataset in place
	herr_t status_write = H5Dwrite(dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_in);

	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	if(status_write>0)
		return 0;
	else
		return -1;
	return 0;
}


int h5md_delete_dataset(struct h5md_file* file, char* absolute_name_of_dataset){ //absolute_name_of_dataset="absolute_path/name_of_dataset"
	hid_t dataset_id=H5Dopen(file->file_id, absolute_name_of_dataset,H5P_DEFAULT); //get dataset_id

	herr_t status_delete=H5Ldelete(dataset_id, absolute_name_of_dataset ,H5P_DEFAULT);
	if(status_delete>0)
		return 0;
	else
		return -1;
}

int h5md_write_attribute(){

}


int h5md_delete_attribute(){

}

int h5md_copy_dataset(){

}

