#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "libh5md.h"

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

	//TODO time dependent data is stored in attributes -> H5A http://www.hdfgroup.org/HDF5/doc/RM/RM_H5A.html#Annot-Open
	//int is_dataset;	//for time independent data
	//int is_attribute;	//for time dependent data
};


// opens the file, creates the internal structure and goes to the first timestep
// you have to use double pointers in order to be able to change a pointer in a foreing function
int h5md_open(struct h5md_file** _file, const char *filename){
	struct h5md_file *file = malloc(sizeof(struct h5md_file));

	file->file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
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
		free(file);
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
int h5md_get_timestep(struct h5md_file* file, int* natoms, double **coords){
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
	//data in memory have the order data_out[3*atom_nr+coord_i]
	double data_out [file->natoms * file->nspacedims];
	
	H5Dread (file->pos_dataset_id, file->datatype, memspace_id, dataspace_id, H5P_DEFAULT, data_out);

	float* data_out_float= malloc(sizeof(float)*file->natoms * file->nspacedims);

	//nasty conversion from double to float
	for(int i=0;i<file->natoms*file->nspacedims;i++){
		data_out_float[i]=(float) data_out[i] ;
	}
	*coords=data_out_float;

	//modify timestep
	int current_time;
	h5md_get_current_time(file,&current_time);
	int status_seek= h5md_seek_timestep(file, current_time+1);
	
	//set natoms to the actual number of atoms
	natoms=&(file->natoms);
	
	//close resources
	H5Sclose(memspace_id);
	
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
int h5md_read_timeindependent_dataset_automatically(struct h5md_file* file, char* dataset_name, void** _data_out, H5T_class_t type_class_out){
	hid_t dataset_id = H5Dopen2(file->file_id, dataset_name,H5P_DEFAULT);	
	/*
	* Get datatype and dataspace handles and then query
	* dataset class, order, size, rank and dimensions.
	*/
	hid_t datatype  = H5Dget_type(dataset_id);     /* datatype handle */
	type_class_out     = H5Tget_class(datatype);
	//H5T_order_t order     = H5Tget_order(datatype);
	size_t size_datatype  = H5Tget_size(datatype);
	hid_t dataspace_id = H5Dget_space(dataset_id);    /* dataspace handle */
	int rank_dataset      = H5Sget_simple_extent_ndims(dataspace_id);
	unsigned long long int dims_dataset[rank_dataset];
	H5Sget_simple_extent_dims(dataspace_id, dims_dataset, NULL);

	switch (type_class_out) {
	case H5T_INTEGER:
		if(dataset_id<0){
			return -1;
		}else{
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
	case H5T_FLOAT:
		if(dataset_id<0){
			return -1;
		}else{
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
	case H5T_STRING:
		if(dataset_id<0){
			return -1;
		}else{
			//load content of name
			//get datatype and its size
			hid_t filetype= H5Dget_type (dataset_id);
			size_t sdim = H5Tget_size (filetype);
			sdim++;// Make room for null terminator
			/*
			* Get dataspace and allocate memory for read buffer.  This is a
			* two dimensional dataset so the dynamic allocation must be done
			* in steps.
			*/
			H5Sget_simple_extent_dims(dataspace_id, dims_dataset, NULL);
			// Allocate array of pointers to rows.
			char **data_out = (char **) malloc (dims_dataset[0] * sizeof (char *));		
			// Allocate space for integer data.
			data_out[0] = (char *) malloc (dims_dataset[0] * sdim * sizeof (char));
			// Set the rest of the pointers to rows to the correct addresses.
			for (int i=1; i<dims_dataset[0]; i++)
				data_out[i] = data_out[0] + i * sdim;
			// Create the memory datatype.
			hid_t memtype = H5Tcopy (H5T_C_S1);
			H5Tset_size (memtype, sdim);
			// Read the data.
			H5Dread (dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out[0]);
			//close
			H5Dclose (dataset_id);
			H5Sclose (dataspace_id);
			H5Tclose (filetype);
			H5Tclose (memtype);
			
			*(_data_out)=data_out;
			return 0;
		}
	break;

	default:
		if(dataset_id<0){
			return -1;
			printf("Dataset contains datatype that is not H5T_INTEGER, H5T_FLOAT or H5T_STRING. Not implemented case.\n");
		}
	break;
	}
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


