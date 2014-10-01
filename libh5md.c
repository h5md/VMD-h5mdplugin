#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "libh5md.h"

#include <libgen.h>	//for use of basename
#include <unistd.h>	//for use of getlogin_r
#include <math.h>	//for use of acos

#define TRUE	1
#define FALSE	0
#define PI 3.14159265

/*
Return value 0 means everything is fine
*/

//box information
typedef struct h5md_box{
	float A;
	float B;
	float C;
	float alpha;
	float beta;
	float gamma;
} h5md_box;


//group handling
typedef struct h5md_group{
	char* group_path;
	int natoms_group;
	hid_t pos_dataset_id;
	hid_t species_dataset_id;
	int* atom_ids;	//TODO use it
	hid_t datatype;
	H5T_class_t t_class;
	H5T_order_t order;
	size_t size_datatype;
	hid_t dataspace;
	int rank_dataset;
	int nspacedims;
	h5md_box* boxes;
} h5md_group;

//file handling
struct h5md_file{
	hid_t file_id;
	h5md_group* groups;
	int ngroups;	
	int natoms;
	int ntime;
	char *last_error_message;
	int current_time;	//for h5md_seek_timestep()
};

//declaration of "private" functions, cannot be accessed from outside the library
//declaration here in order to be able to sort the content of the file according to it's "importance"
int initialize_h5md_struct(struct h5md_file* file);
int modify_information_about_file_content(struct h5md_file* file, char* group_name);
int check_compatibility(struct h5md_file* file, hid_t new_pos_dataset_id);
herr_t check_for_pos_dataset( hid_t g_id, const char* obj_name, const H5L_info_t* info, void* _file);
h5md_box* get_box_information(struct h5md_file* file, int group_number);

//declaration of boring helper functions
char* concatenate_strings(const char* string1,const char* string2);
int max(int a, int b);
float calculate_length_of_vector(double* vector, int dimensions);
float calculate_angle_between_vectors(double* vector1, double* vector2, int dimensions);

int discover_all_groups(struct h5md_file* file){
	H5Lvisit(file->file_id, H5_INDEX_NAME, H5_ITER_NATIVE, check_for_pos_dataset, (void*) file);	//discover all groups with position datasets
	return 0;
}

//checks whether current object is a position group, then adds the dataset_id to the pos_dataset_id array in the h5md_file file
herr_t check_for_pos_dataset( hid_t g_id, const char* obj_name, const H5L_info_t* info, void* _file){
	int status=0;
	struct h5md_file* file=_file;

	if(! strcmp(basename((char*)obj_name), "position")){
		H5G_stat_t statbuf;
		H5Gget_objinfo(g_id , obj_name , FALSE, &statbuf);	//check, whether the object is a group 
		if(statbuf.type==H5G_GROUP){
			char* full_path_position_dataset=concatenate_strings((const char*) obj_name,(const char*) "/value");	
			hid_t pos_dataset_id=H5Dopen2(file->file_id, full_path_position_dataset ,H5P_DEFAULT);
			free(full_path_position_dataset);

			if(pos_dataset_id>=0)
				status=modify_information_about_file_content(file, dirname((char*)obj_name));
				printf("Position dataset found in group /%s.\n", obj_name);
		}
	}
	return status;	//if status is 0 search for other position datasets continues. If status is negative, search is aborted.
}

int modify_information_about_file_content(struct h5md_file* file, char* group_name){
	int status=-1;
	
	//get pos_dataset_id
	char* full_path_position_dataset=concatenate_strings((const char*) group_name,(const char*) "/position/value");	
	hid_t pos_dataset_id=H5Dopen2(file->file_id, full_path_position_dataset ,H5P_DEFAULT);
	free(full_path_position_dataset);

	//get species_dataset_id
	char* full_path_species_dataset=concatenate_strings((const char*) group_name,(const char*) "/species");	
	hid_t species_dataset_id=H5Dopen2(file->file_id, full_path_species_dataset ,H5P_DEFAULT);
	free(full_path_species_dataset);

	if(check_compatibility(file, pos_dataset_id)==0){
		file->ngroups+=1;
		h5md_group* groups=(h5md_group*) realloc(file->groups, sizeof(h5md_group)*(file->ngroups));	//effectively appends one entry to array

		groups[file->ngroups-1].pos_dataset_id=pos_dataset_id;
		groups[file->ngroups-1].species_dataset_id=species_dataset_id;		
	
		/*
		* Get datatype and dataspace handles and then query
		* dataset class, order, size, rank and dimensions. Since all datasets are checked to be compatible do this only for the first dataset
		*/

		groups[file->ngroups-1].datatype  = H5Dget_type(pos_dataset_id);     // datatype handle
		groups[file->ngroups-1].t_class     = H5Tget_class(groups[file->ngroups-1].datatype);
		groups[file->ngroups-1].order     = H5Tget_order(groups[file->ngroups-1].datatype);

		groups[file->ngroups-1].size_datatype  = H5Tget_size(groups[file->ngroups-1].datatype);
		groups[file->ngroups-1].dataspace = H5Dget_space(pos_dataset_id);	//dataspace handle
		groups[file->ngroups-1].rank_dataset      = H5Sget_simple_extent_ndims(groups[file->ngroups-1].dataspace);
		hsize_t dims_out[groups[file->ngroups-1].rank_dataset];
		H5Sget_simple_extent_dims(groups[file->ngroups-1].dataspace, dims_out, NULL);
		file->ntime = dims_out[0];
		file->natoms += dims_out[1];
		groups[file->ngroups-1].nspacedims = dims_out[2];
		groups[file->ngroups-1].natoms_group=dims_out[1];
		groups[file->ngroups-1].group_path=group_name;


		file->groups=groups; // assign file->groups here since get_box_information() uses all group_path entries
		file->groups[file->ngroups-1].boxes=get_box_information(file, file->ngroups-1);
		//file->groups=groups;

		status=0;
	}else{
		printf("position datasets are not compatible\n");
		status=-1;	
	}
	return status;
}



int check_compatibility(struct h5md_file* file, hid_t new_pos_dataset_id){
	hid_t new_dataspace = H5Dget_space(new_pos_dataset_id);	//dataspace handle
	int new_rank_dataset = H5Sget_simple_extent_ndims(new_dataspace);
	hsize_t dims_out[new_rank_dataset];
	H5Sget_simple_extent_dims(new_dataspace, dims_out, NULL);
	int ntime_new_pos_dataset =(int) dims_out[0];

	if(file->ntime==ntime_new_pos_dataset || (file->ngroups==0 && new_pos_dataset_id>0))	//ngroups=0 from initialization
		return 0;
	else
		return -1;
}

// opens the file, creates the internal structure and goes to the first timestep
// you have to use double pointers in order to be able to change a pointer in a foreing function
int h5md_open(struct h5md_file** _file, const char *filename, int can_write){
	struct h5md_file *file = malloc(sizeof(struct h5md_file));

	if(can_write==TRUE)
		file->file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT); //read&write access
	else
		file->file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	h5md_hide_hdf5_error_messages();

	initialize_h5md_struct(file);
	discover_all_groups(file);

	*_file = file;
	if(file->file_id <0){
		return -1;
	}else{
		return 0;
	}
}

// close file, datasets and frees the internal structure
int h5md_close(struct h5md_file* file){
	if(file!=NULL){
		if(file->ngroups>0){
			for(int i=0; i<file->ngroups; i++){//close datasets
				H5Dclose(file->groups[i].pos_dataset_id);
				H5Dclose(file->groups[i].species_dataset_id);
			}
			free(file->groups);
		}
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
int h5md_get_timestep(struct h5md_file* file, int* natoms, float **coords){

	float* data_out= malloc(sizeof(float)*file->natoms * 3); //allocate space for data in memory, which have the order data_out[3*atom_nr+coord_i] Remark: 3 could be the maximum over i of file->groups[i].nspacedims

	int previous_atoms=0;
	for(int i=0; i<file->ngroups; i++){//go through all groups
		hid_t dataspace_id=H5Dget_space(file->groups[i].pos_dataset_id); //Define dataset dataspace in file.

		/* 
		* Define hyperslab in the dataset. 
		*/
		hsize_t dataset_slab_offset[file->groups[i].rank_dataset];
		dataset_slab_offset[0] = file->current_time;
		dataset_slab_offset[1] = 0;
		dataset_slab_offset[2] = 0;

		hsize_t dataset_slab_count[file->groups[i].rank_dataset];
		dataset_slab_count[0] = 1;
		dataset_slab_count[1] = file->groups[i].natoms_group;
		dataset_slab_count[2] = file->groups[i].nspacedims;
		H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, dataset_slab_offset, NULL, dataset_slab_count, NULL);
	
		/*
		* Define memory dataspace.
		*/
		int rank=1; //linear data representation
		hsize_t dimsm[rank];
		dimsm[0]=file->natoms* file->groups[i].nspacedims;

		hid_t memspace_id = H5Screate_simple(rank,dimsm,NULL);

		/* 
		* Define memory hyperslab. 
		*/
		hsize_t offset_out[rank];
		hsize_t count_out[rank];
		if(i<1){//set offset
			offset_out[0]=0;
		}else{
			offset_out[0]=previous_atoms*file->groups[i].nspacedims;
		}
		count_out[0]=file->groups[i].natoms_group*file->groups[i].nspacedims;	
		previous_atoms+=file->groups[i].natoms_group;
		H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
		

		/*
		* Read data from hyperslab in the file into the hyperslab in memory
		*/
		
		hid_t wanted_memory_datatype = H5T_NATIVE_FLOAT;
		H5Dread (file->groups[i].pos_dataset_id, wanted_memory_datatype, memspace_id, dataspace_id, H5P_DEFAULT, data_out);
		H5Sclose(memspace_id); //close resources
		H5Sclose(dataspace_id);
	}

	*coords=data_out;
	natoms=&(file->natoms); //set natoms to the actual number of atoms
	
	int current_time;
	h5md_get_current_time(file,&current_time);
	int status_seek= h5md_seek_timestep(file, current_time+1); //modify timestep
	
	
	if(status_seek==0){
		return 0;
	}else{
		return -1;
	}
}


// reads the next timestep and writes the data to coords (for which the memory has been allocated externally before the call of this function) iff natoms is the number of atoms in the timestep
int h5md_read_timestep(struct h5md_file* file, int natoms, float* coords){
	float* data_out_temp;
	int status_read =h5md_get_timestep(file, &natoms, (float**) &(data_out_temp));
	memcpy(coords, data_out_temp, sizeof(float)*natoms*3);
	free(data_out_temp);
	return status_read;
}

// internally reads the box information of a given group into the memory
h5md_box* get_box_information(struct h5md_file* file, int group_number){
	h5md_box* boxes=malloc(sizeof(h5md_box)*file->ntime);

	//check whether box_dataset is timedependent, if it is timedependent use it, otherwise copy the box information ntime times
	//try to open time-dependent box dataset, get box_dataset_timedependent_id
	char* full_path_box_dataset_timedependent=concatenate_strings((const char*) file->groups[group_number].group_path,(const char*) "/box/value");	
	hid_t box_timedependent_dataset_id=H5Dopen2(file->file_id, full_path_box_dataset_timedependent ,H5P_DEFAULT);
	//try to open time-independent box dataset, get box_dataset_timeindependent_id
	char* full_path_box_dataset_timeindependent=concatenate_strings((const char*) file->groups[group_number].group_path,(const char*) "/box/edges");	
	hid_t box_timeindependent_dataset_id=H5Dopen2(file->file_id, full_path_box_dataset_timeindependent ,H5P_DEFAULT);

	if(box_timedependent_dataset_id>0){
		//timedependent dataset exists, use it
		//read timedependent dataset 

		float data_box[file->ntime*3*3];
		H5Dread(box_timedependent_dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_box);
		for(int i; i<file->ntime; i++){
			h5md_box box;
			double vector_a[3]={data_box[i*9+0], data_box[i*9+1], data_box[i*9+2]};
			double vector_b[3]={data_box[i*9+3], data_box[i*9+4], data_box[i*9+5]};
			double vector_c[3]={data_box[i*9+6], data_box[i*9+7], data_box[i*9+8]};
			//process to angles and lengths
			box.A=calculate_length_of_vector(vector_a,3);
			box.B=calculate_length_of_vector(vector_b,3);
			box.C=calculate_length_of_vector(vector_c,3);
			box.alpha= calculate_angle_between_vectors(vector_b,vector_c,3);
			box.beta= calculate_angle_between_vectors(vector_a,vector_c,3);
			box.gamma= calculate_angle_between_vectors(vector_a,vector_b,3);
			boxes[i]=box;
		}

	}else{
		if(box_timeindependent_dataset_id>0){
			h5md_box box;
			//read timeindependent dataset
			//decided whether box is cubic (dataset contains a vector) or triclinic (dataset contains a matrix)
			int dims_box;
			int is_cubic=h5md_get_length_of_one_dimensional_dataset(file, full_path_box_dataset_timedependent, &dims_box);
			float* data_box;
			H5T_class_t box_class_out;
			int status=h5md_read_timeindependent_dataset_automatically(file, full_path_box_dataset_timeindependent,(void**) &data_box, &box_class_out);
			
			//process to angles and lengths			
			if(is_cubic==0){
				//box is cubic implies vector
				box.A=data_box[0];
				box.B=data_box[1];
				box.C=data_box[2];
				box.alpha=90;
				box.beta=90;
				box.gamma=90;
			}else{
				//box is triclinic implies matrix
				//VMD expects system to be 3dimensional -> assume 3x3 matrix
				double vector_a[3]={data_box[0],data_box[1],data_box[2]};
				double vector_b[3]={data_box[3],data_box[4],data_box[5]};
				double vector_c[3]={data_box[6],data_box[7],data_box[8]};
				
				//according to VMD's molfile_timestep_t documentation
				box.A=calculate_length_of_vector(vector_a,3);
				box.B=calculate_length_of_vector(vector_b,3);
				box.C=calculate_length_of_vector(vector_c,3);
				box.alpha= calculate_angle_between_vectors(vector_b,vector_c,3);
				box.beta= calculate_angle_between_vectors(vector_a,vector_c,3);
				box.gamma= calculate_angle_between_vectors(vector_a,vector_b,3);
			}


			//copy timeindependend dataset ntimes
			for(int i=0; i<file->ntime; i++){
				boxes[i]=box;
			}
		}else{
			printf("No box information found");	
		}
	}

	
	free(full_path_box_dataset_timedependent);
	free(full_path_box_dataset_timeindependent);

	return boxes; //TODO make correct

}

int free_box_information(struct h5md_file* file){
	for(int i=0;i<file->ngroups;i++){
		free(file->groups[i].boxes);	
	}
	return 0;
}


int h5md_get_box_information(struct h5md_file* file, float* out_box_information){
	//only use information of first group (file->groups[0]...) since VMD only supports one box
	out_box_information[0]=file->groups[0].boxes[file->current_time].A;
	out_box_information[1]=file->groups[0].boxes[file->current_time].B;
	out_box_information[2]=file->groups[0].boxes[file->current_time].C;
	out_box_information[3]=file->groups[0].boxes[file->current_time].alpha;
	out_box_information[4]=file->groups[0].boxes[file->current_time].beta;
	out_box_information[5]=file->groups[0].boxes[file->current_time].gamma;
	if(file->current_time>=file->ntime) //free box information
		free_box_information(file);
	return 0;
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
		//printf("Dataset %s could not be opened.\n", dataset_name);
		status=-1;
	}else{
		switch (*type_class_out) {
		case H5T_INTEGER:{
			hid_t wanted_memory_datatype = H5T_NATIVE_INT;			
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
			hid_t wanted_memory_datatype = H5T_NATIVE_FLOAT;
			//determine needed size
			int needed_size=dims_dataset[0];
			int len_dims_dataset=sizeof(dims_dataset)/sizeof(dims_dataset[0]);
			for(int i=1; i<len_dims_dataset; i++){
				needed_size*=dims_dataset[i];
			}
			float* data_out=(float*) malloc(sizeof(size_datatype)*needed_size);
			int status_read=H5Dread(dataset_id, wanted_memory_datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out);
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

int h5md_free_timeindependent_dataset_automatically(H5T_class_t type_class, void* old_data_out, int length_array_of_strings){
	switch (type_class) {
	case H5T_INTEGER:	//stacking case conditions before break produces OR condition
	case H5T_FLOAT:
		free(old_data_out);
	break;
	case H5T_STRING:{//TODO differentiate between fixed and variable length strings for correct free

//The following is only for variable length strings
/*		char** old_data_string= (char**) old_data_out;*/
/*		for(int i=0;i<length_array_of_strings;i++){*/
/*			free(((char**)old_data_out)[i]);*/
/*		}*/

//This is only for fixed length strings
		free(((char**)old_data_out)[0]);
		free(old_data_out); 
	}
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
		//printf("Dataset %s is not one dimensional.\n", dataset_name);	
		return -1;	
	}
}


void h5md_hide_hdf5_error_messages(){
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
}

void h5md_show_hdf5_error_messages(){
	H5Eset_auto(H5E_DEFAULT, (H5E_auto_t) H5Eprint, stderr);
}

int h5md_get_file_id(struct h5md_file *file, hid_t *file_id){
	*file_id=file->file_id;
	if(file_id !=0)
		return 0;
	else
		return -1;
}

int h5md_get_all_species_infromation(struct h5md_file *file, int** species_infromation_out){
	//derived from h5md_get_timestep() above
	//TODO generalize function to h5md_get_dataset_information_from_all_groups()

	int* data_out= malloc(sizeof(int)*file->natoms * file->ngroups); //allocate space for data in memory, which have the order data_out[atom_nr]

	int previous_atoms=0;
	for(int i=0; i<file->ngroups; i++){//go through all groups
		hid_t dataspace_id=H5Dget_space(file->groups[i].species_dataset_id); //Define dataset dataspace in file.

		/* 
		* Define hyperslab in the dataset. 
		*/
		hsize_t dataset_slab_offset[file->groups[i].rank_dataset];
		dataset_slab_offset[0] = 0;
		hsize_t dataset_slab_count[file->groups[i].rank_dataset];
		dataset_slab_count[0] = file->groups[i].natoms_group;
		H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, dataset_slab_offset, NULL, dataset_slab_count, NULL);
	
		/*
		* Define memory dataspace.
		*/
		int rank=1; //linear data representation
		hsize_t dimsm[rank];
		dimsm[0]=file->natoms;
		hid_t memspace_id = H5Screate_simple(rank,dimsm,NULL);
	
	
		/* 
		* Define memory hyperslab. 
		*/
		hsize_t offset_out[rank];
		hsize_t count_out[rank];
		offset_out[0]=previous_atoms;
		count_out[0]=file->groups[i].natoms_group;
		
		previous_atoms+=file->groups[i].natoms_group;
		H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
		

		/*
		* Read data from hyperslab in the file into the hyperslab in memory
		*/
		hid_t wanted_memory_datatype = H5T_NATIVE_INT;
		herr_t status=H5Dread(file->groups[i].species_dataset_id, wanted_memory_datatype, memspace_id, dataspace_id, H5P_DEFAULT, data_out);
		H5Sclose(memspace_id); //close resources
		H5Sclose(dataspace_id);
	}

	*species_infromation_out=data_out;
	
	return 0;
}

/*int h5md_get_dataset_information_from_all_groups(struct h5md_file *file, char* relative_dataset_name , void** data_out){*/

/*	return 0;*/
/*}*/


/* write operations */
//creates a h5md_file iff it does not exist yet
int h5md_create_file(struct h5md_file **_file, const char* filename){
	struct h5md_file *file = (struct h5md_file*) malloc(sizeof(struct h5md_file));
	initialize_h5md_struct(file);
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


//appends data to dataset, create dataset if it is not yet existing
int h5md_append_dataset(struct h5md_file *file, char* absolute_name_of_dataset, hid_t datatype, void* data_in, int rank_in, hsize_t* dims_in){
	int status;
	//check existence of dataset
	h5md_hide_hdf5_error_messages();
	hid_t dataset_id = H5Dopen2(file->file_id, absolute_name_of_dataset, H5P_DEFAULT);
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
			H5Sclose(dataspace);
			H5Sclose(filespace);
			free(offset);
			free(size);	
		}
		H5Tclose(datatype_read);
		H5Sclose(dataspace_id);
		H5Dclose(dataset_id);
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

int h5md_set_author(struct h5md_file* file, char* name, char* email_address){
	hid_t link_crt_plist = H5Pcreate(H5P_LINK_CREATE);
	H5Pset_create_intermediate_group(link_crt_plist, TRUE);	// Set flag for intermediate group creation

	H5Gcreate(file->file_id, "/h5md/author", link_crt_plist, H5P_DEFAULT, H5P_DEFAULT);

	herr_t status_name =-1;
	int status_username_os =-1;
	if(name == NULL){
		char* username=getlogin();
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


//use H5LT http://www.hdfgroup.org/HDF5/hdf5_hl/doc/RM_hdf5lt.html#H5LTmake_dataset with easy interface


//int h5md_set_version(struct h5md_file* file, char* name){

//}

/*
int h5md_delete_attribute(){

}

int h5md_rename_dataset(){

}
*/






/* boring helper functions */

char* concatenate_strings(const char* string1, const char* string2){
	char* concat_string=malloc(strlen(string1)+strlen(string2)+1);	//+1 for 0 termination
	concat_string=strcpy(concat_string, string1);
	concat_string=strcat(concat_string, string2);
	return concat_string;
}

int initialize_h5md_struct(struct h5md_file* file){
	file->ngroups=0; //initialize the number of position datasets to 0
	file->current_time=0;	//initialize current time to 0
	file->groups=NULL;
	return 0;
}

int max(int a, int b){
	if(a>b)
		return a;
	else
		return b;
}

float calculate_length_of_vector(double* vector, int dimensions){
	float length=0.0;
	for(int i=0;i<dimensions;i++){
	length+=vector[i]*vector[i];
	}
	length=sqrt(length);
	return length;
}

//returns 0 if "vector" is the zero vector, otherwise nonzero
int is_zero_vector(double* vector, int dimensions){
	float absolute_error=0.00001;
	int number_of_nonzero_entries=0;
	for(int i=0;i<dimensions;i++){
		if(fabs(vector[i]-0)>absolute_error){
				number_of_nonzero_entries+=1;
				break;
		}	
	}
	return number_of_nonzero_entries;
}

float calculate_angle_between_vectors(double* vector1, double* vector2, int dimensions){
	if(is_zero_vector(vector1, dimensions)==0 || is_zero_vector(vector2, dimensions)==0){
		printf("Error: Cannot calculate angle to the zero vector which has length 0\n");
		return -1.0;
	}
	float angle=0.0;
	for(int i=0;i<dimensions;i++){
		angle+=vector1[i]*vector2[i];	
	}
	float len_vector1=calculate_length_of_vector(vector1, dimensions);
	float len_vector2=calculate_length_of_vector(vector2, dimensions);
	angle=acos(angle/(len_vector1*len_vector2))*180/PI;	//*180/PI converts radian to degree
	return angle;
}



