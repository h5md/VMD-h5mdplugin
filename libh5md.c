#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "libh5md.h"

#include <unistd.h> //for set author to username

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
int h5md_open(struct h5md_file** _file, const char *filename, int can_write){
	struct h5md_file *file = malloc(sizeof(struct h5md_file));

	if(can_write==0)
		file->file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT); //read&write access
	else
		file->file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	h5md_hide_hdf5_error_messages();
	file->pos_dataset_id = H5Dopen2(file->file_id, "/particles/atoms/position/value",H5P_DEFAULT);
	
	if(file->pos_dataset_id>0){
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
	}
	
	*_file = file;
	if(file->file_id <0){
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

//set number of atoms iff this number is constant during time
int h5md_set_natoms(struct h5md_file* file, int natoms){
	file->natoms=natoms;
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
	if(i<=ntime){
		file->current_time=i;
		return 0;
	}else{
		return -1;
	}
}



// reads the next timestep, allocates coords and sets natoms to the number of atoms
int h5md_get_timestep(struct h5md_file* file, int* natoms, float **coords){ //TODO rewrite should not store the data in h5md_file
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
	hsize_t dimsm[rank];
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

	natoms=&(file->natoms); //set natoms to the actual number of atoms
	
	int current_time;
	h5md_get_current_time(file,&current_time);
	int status_seek= h5md_seek_timestep(file, current_time+1); //modify timestep
	H5Sclose(memspace_id); //close resources
	
	if(status_seek==0){
		return 0;
	}else{
		return -1;
	}
}


// reads the next timestep and writes the data to coords iff natoms is the number of atoms in the timestep
int h5md_read_timestep(struct h5md_file* file, int natoms, float* coords){
	int status_read =h5md_get_timestep(file,&natoms,&coords); //TODO h5md_get timestep allocates new memory, h5md_read_timestep should not allocate new memory but write into existing memory -> this function may not use h5md_get_timestep if done correctly
	if(status_read==0)
		return 0;
	else
		return -1;
}


//read timeindependent dataset automatically
int h5md_read_timeindependent_dataset_automatically(struct h5md_file* file, char* dataset_name, void** _data_out, H5T_class_t* type_class_out){
	int status;
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
		status=-1;
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
			int status_read=H5Dread(dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out);
			status=status_read;
			*(_data_out)=data_out;
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
			int status_read=H5Dread(dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out);
			status=status_read;
			*(_data_out)=data_out;
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
					hid_t memtype = H5Tcopy (H5T_C_S1);
					H5Tset_size (memtype, size_datatype);
					int status_read=H5Dread (dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out[0]);
					status=status_read;
					H5Tclose (memtype); //close memtype
				}else{
					//string length is fixed
					//Create the memory datatype.
					hid_t memtype = H5Tcopy (H5T_C_S1);
					H5Tset_size (memtype, H5T_VARIABLE);
					int status_read=H5Dread (dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out);
					status=status_read;
					H5Tclose (memtype); //close memtype
				}

				*(_data_out)=data_out;
			}
		break;

		default:{
				printf("Dataset contains datatype that is not H5T_INTEGER, H5T_FLOAT or H5T_STRING. Not implemented case.\n");
				status=-1;
			}
		break;
		}
	}

	//close resources
	H5Dclose (dataset_id);
	H5Sclose (dataspace_id);
	H5Tclose (datatype);
	return status;
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

int h5md_read_timeindependent_dataset_int(struct h5md_file* file, char* dataset_name, int ** _data_out){
	H5T_class_t type_class_out;
	int status_read=h5md_read_timeindependent_dataset_automatically(file, dataset_name, (void**)_data_out, &type_class_out);

	if( (status_read!=0) || (type_class_out =! H5Tget_class(H5T_NATIVE_INT)))
		return -1;
	else
		return 0;
}

int h5md_get_length_of_one_dimensional_dataset(struct h5md_file *file,char *dataset_name, int *length_of_dataset){
	hid_t dataset_id = H5Dopen2(file->file_id, dataset_name,H5P_DEFAULT);
	hid_t dataspace_id = H5Dget_space(dataset_id);    /* dataspace handle */
	hsize_t rank_dataset      = H5Sget_simple_extent_ndims(dataspace_id);
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

int h5md_get_file_id(struct h5md_file *file, hid_t *file_id){
	*file_id=file->file_id;
	if(file_id !=0)
		return 0;
	else
		return -1;
}




/* write operations */
//creates a h5md_file iff it does not exist yet
int h5md_create_file(struct h5md_file **_file, const char* filename){
	struct h5md_file *file = (struct h5md_file*) malloc(sizeof(struct h5md_file));
	/* Create a new file using default properties. */
	hid_t file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT); //H5F_ACC_EXCL <-> fail if file already exists, alternative H5F_ACC_TRUNC <-> overwrite file
	file->file_id=file_id;
	h5md_set_author(file, NULL, NULL); //sets author name by default to the user account name, can be overwritten by another call of h5md_set_author(file,"username");
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
	int status;
	hid_t link_crt_plist = H5Pcreate(H5P_LINK_CREATE);
	H5Pset_create_intermediate_group(link_crt_plist, TRUE);	// Set flag for intermediate group creation

	//Modify dataset creation properties, i.e. enable chunking. 
	//Chunking has to be enabled at the creation time of the dataset!
	//"The storage properties cannot be changed after the dataset is created"
	hid_t cparms = H5Pcreate (H5P_DATASET_CREATE);
	H5Pset_chunk( cparms, rank_in, dims_in); // enable chunking
	
	
	switch(H5Tget_class(datatype)){
		case H5T_INTEGER:
		{	
			int filler =-1;
			H5Pset_fill_value(cparms, datatype, &filler);//set fill value
		break;
		}	
		case H5T_FLOAT:
		{
			float filler =-1.0;
			H5Pset_fill_value(cparms, datatype, &filler);//set fill value
		break;
		}	
		case H5T_STRING:
		{	
			char* filler ="filler";
			H5Pset_fill_value(cparms, datatype, &filler);//set fill value
		break;
		}
		default:
			printf("datatype not implemented\n");
		break;
	}
	
	
	
	hsize_t* maxdims=malloc(sizeof(hsize_t)*rank_in);
	for(int i=0;i<rank_in;i++){
		maxdims[i]=H5S_UNLIMITED;	
	}
	hid_t dataspace_id=H5Screate_simple(rank_in, dims_in, maxdims);

	hid_t dataset_id = H5Dcreate(file->file_id, absolute_name_of_dataset, datatype, dataspace_id, link_crt_plist, cparms, H5P_DEFAULT); //create dataset
	herr_t status_write = H5Dwrite(dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_in);

	free(maxdims);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);

	if(status_write>=0)
		status= 0;
	else
		status= -1;
	return status;
}

int max(int a, int b){
	if(a>b)
		return a;
	else
		return b;
}

//appends data to dataset, create dataset if it is not yet existing
int h5md_append_dataset(struct h5md_file *file, char* absolute_name_of_dataset, hid_t datatype, void* data_in, int rank_in, hsize_t* dims_in){
	int status;
	//check existence of dataset
	h5md_hide_hdf5_error_messages();
	hid_t dataset_id = H5Dopen2(file->file_id, absolute_name_of_dataset ,H5P_DEFAULT);
	h5md_show_hdf5_error_messages();

	if(dataset_id>0){
		hid_t datatype_read  = H5Dget_type(dataset_id);     // datatype handle
		H5T_class_t type_class_read     = H5Tget_class(datatype_read);
		hid_t dataspace_id = H5Dget_space(dataset_id);    //dataspace handle
		int rank_dataset_read = H5Sget_simple_extent_ndims(dataspace_id);
		unsigned long long int dims_dataset_read[rank_dataset_read];
		H5Sget_simple_extent_dims(dataspace_id, dims_dataset_read, NULL);

		//check consistency // 
		if(rank_in != rank_dataset_read || H5Tget_class(datatype) != type_class_read){
			printf("Data cannot be appended to dataset.\n");
			status=-1;
		}else{
			//append to dataset, compare http://www.hdfgroup.org/HDF5/doc/H5.intro.html#Intro-PMCreateExtendible and linked example there
			unsigned long long int* size=malloc(sizeof(unsigned long long int)*rank_in);
			for(int i=0;i<rank_in;i++){
				if(i==0)
					size[i]=dims_dataset_read[0]+dims_in[0];
				else
					size[i]=max(dims_in[i],dims_dataset_read[i]);
			}
			H5Dextend (dataset_id, size);

			//Select a hyperslab.
			hid_t filespace = H5Dget_space(dataset_id);
			hsize_t* offset=malloc(sizeof(hsize_t)*rank_in);
			for(int i=0;i<rank_in;i++){
				if(i==0)
					offset[i]=dims_dataset_read[0];
				else
					offset[i]=0;
			}
			H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dims_in, NULL);  
			hid_t dataspace = H5Screate_simple(rank_in, dims_in, NULL); //Define memory space
			herr_t status_write = H5Dwrite(dataset_id, datatype_read, dataspace, filespace, H5P_DEFAULT, data_in); //Write the data to the hyperslab.
			if(status_write>=0)
				status=0;
			else
				status=-1;
			H5Dclose(dataset_id);
			H5Sclose(dataspace);
			H5Sclose(filespace);
			free(offset);
			free(size);	
		}


	}else{
		//create dataset if it is not yet existing
		status=h5md_write_dataset(file, absolute_name_of_dataset, datatype, data_in, rank_in, dims_in);
	}
	return status;
}

int get_fill_value(struct h5md_file* file, char* absolute_name_of_dataset, void* filler){
	hid_t dataset_id=H5Dopen(file->file_id, absolute_name_of_dataset, H5P_DEFAULT);
	hid_t plist2 = H5Dget_create_plist(dataset_id);
	herr_t status;
	hid_t datatype  = H5Dget_type(dataset_id);     // datatype handle
	switch(H5Tget_class(datatype)){
		case H5T_INTEGER:
		{	
			status=H5Pget_fill_value(plist2, H5T_NATIVE_INT, (int*) filler);
		break;
		}	
		case H5T_FLOAT:
		{
			status=H5Pget_fill_value(plist2, H5T_NATIVE_FLOAT, (float*) filler);
		break;
		}	
		case H5T_STRING:
		{	
			status=H5Pget_fill_value(plist2, H5T_STRING, (char*) filler);
		break;
		}
		default:
			printf("datatype not implemented\n");
		break;
	}

	
	if(status>=0)
		return 0;
	else
		return -1;
}


int h5md_delete_dataset(struct h5md_file* file, char* absolute_name_of_dataset){ //absolute_name_of_dataset="absolute_path/name_of_dataset"
	//HDF5 does not at this time provide an easy mechanism to remove a dataset from a file or to reclaim the storage space occupied by a deleted object.

	hid_t dataset_id=H5Dopen(file->file_id, absolute_name_of_dataset,H5P_DEFAULT); //get dataset_id

	herr_t status_delete=H5Ldelete(dataset_id, absolute_name_of_dataset ,H5P_DEFAULT);
	if(status_delete>=0)
		return 0;
	else
		return -1;
}


//TODO but use H5LT http://www.hdfgroup.org/HDF5/hdf5_hl/doc/RM_hdf5lt.html#H5LTmake_dataset with easier interface
/*
int h5md_write_attribute(){
}
*/

int h5md_set_author(struct h5md_file* file, char* name, char* email_address){
	H5Gcreate(file->file_id, "/author", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	herr_t status_name =-1;
	int status_username_os =-1;
	if(name == NULL){
		char username[33];	//current (2014) standard username may be 32 characters long + 1 character for null termination
		status_username_os=getlogin_r(username, sizeof(username));
		status_name=H5LTset_attribute_string(file->file_id, "/author", "name", username);
	}else{
		status_name=H5LTset_attribute_string(file->file_id, "/author", "name", name);
	}
	
	int status_email =-1;
	if(email_address==NULL){
		status_email=H5LTset_attribute_string(file->file_id, "/author", "email", "email was not provided");
	}else{
		status_email=H5LTset_attribute_string(file->file_id, "/author", "email", email_address);	
	}

	if( (status_name>0 || status_username_os==0) && status_email>0 )
		return 0;	
	else
		return -1;	
}

//int h5md_set_version(struct h5md_file* file, char* name){

//}

/*
int h5md_delete_attribute(){

}

int h5md_copy_dataset(){

}
*/

