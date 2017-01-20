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

//#define DEBUG 0

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
	int nspacedims;
	hid_t pos_dataset_id;
	hid_t species_dataset_id;
	hid_t mass_dataset_id;
	hid_t charge_dataset_id;
	hid_t image_dataset_id;
	hid_t id_dataset_id;
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
int get_box_information(struct h5md_file* file, int group_number, int time_i, h5md_box* out_box);
int get_box_vectors(struct h5md_file* file,  int group_i, int time_i, float* vector_a, float* vector_b, float* vector_c);


//declaration of boring helper functions
char* concatenate_strings(const char* string1,const char* string2);
int max(int a, int b);
float calculate_length_of_vector(float* vector, int dimensions);
float calculate_angle_between_vectors(float* vector1, float* vector2, int dimensions);
char* mystrdup(char* str);

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

	//get species_dataset_id for timeindependent species dataset (a timedependent dataset would be located under /species/value)
	char* full_path_species_dataset=concatenate_strings((const char*) group_name,(const char*) "/species");	
	hid_t species_dataset_id=H5Dopen2(file->file_id, full_path_species_dataset ,H5P_DEFAULT);
	free(full_path_species_dataset);
	
	//get mass_dataset_id for timeindependent mass dataset (a timedependent dataset would be located under /mass/value)
	char* full_path_mass_dataset=concatenate_strings((const char*) group_name,(const char*) "/mass");	
	hid_t mass_dataset_id=H5Dopen2(file->file_id, full_path_mass_dataset ,H5P_DEFAULT);
	free(full_path_mass_dataset);
	
	//get charge_dataset_id for timeindependent charge dataset (a timedependent dataset would be located under /charge/value)
	char* full_path_charge_dataset=concatenate_strings((const char*) group_name,(const char*) "/charge");	
	hid_t charge_dataset_id=H5Dopen2(file->file_id, full_path_charge_dataset ,H5P_DEFAULT);
	free(full_path_charge_dataset);

	//get image_dataset_id
	char* full_path_image_dataset=concatenate_strings((const char*) group_name,(const char*) "/image/value");
	hid_t image_dataset_id=H5Dopen2(file->file_id, full_path_image_dataset ,H5P_DEFAULT);
	free(full_path_image_dataset);

	//get id_dataset_id
	char* full_path_id_dataset=concatenate_strings((const char*) group_name,(const char*) "/id/value");
	hid_t id_dataset_id=H5Dopen2(file->file_id, full_path_id_dataset ,H5P_DEFAULT);
	free(full_path_id_dataset);


	if(check_compatibility(file, pos_dataset_id)==0){
		file->ngroups+=1;
		h5md_group* groups=(h5md_group*) realloc(file->groups, sizeof(h5md_group)*(file->ngroups));	//effectively appends one entry to array

		groups[file->ngroups-1].pos_dataset_id=pos_dataset_id;
		groups[file->ngroups-1].species_dataset_id=species_dataset_id;
		groups[file->ngroups-1].mass_dataset_id=mass_dataset_id;
		groups[file->ngroups-1].charge_dataset_id=charge_dataset_id;
		groups[file->ngroups-1].image_dataset_id=image_dataset_id;
		groups[file->ngroups-1].id_dataset_id=id_dataset_id;

		/*
		* Get dataspace handles and then query
		* dataset rank and dimensions. Since all datasets are checked to be compatible do this only for the first dataset
		*/

		hid_t dataspace_id = H5Dget_space(pos_dataset_id);	//dataspace handle
		int rank_dataset      = H5Sget_simple_extent_ndims(dataspace_id);
		hsize_t dims_out[rank_dataset];

		H5Sget_simple_extent_dims(dataspace_id, dims_out, NULL);
		file->ntime = dims_out[0];
		file->natoms += dims_out[1];
		groups[file->ngroups-1].nspacedims = dims_out[2];
		groups[file->ngroups-1].natoms_group=dims_out[1];
		groups[file->ngroups-1].group_path=mystrdup(group_name);
		H5Sclose(dataspace_id);

		file->groups=groups;
		status=0;
	}else{
		printf("position datasets are not compatible\n");
		status=-1;
	}
	return status;
}



int check_compatibility(struct h5md_file* file, hid_t new_pos_dataset_id){
	if(new_pos_dataset_id<0)
		return -1;
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
				H5Dclose(file->groups[i].image_dataset_id);
				H5Dclose(file->groups[i].mass_dataset_id);
				H5Dclose(file->groups[i].charge_dataset_id);
				H5Dclose(file->groups[i].id_dataset_id);
			}
			free(file->groups);
		}
		H5Fflush(file->file_id,H5F_SCOPE_GLOBAL);
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
	if(i>=0 && i<ntime){
		file->current_time=i;
		return 0;
	}else{
		file->current_time=i;
		return -1;
	}
}


//declare struct for binary tree node in order to assign a particle with an id to the correct index in the current position dataset, assumes that the position dataset and the image dataset have the same ordering
typedef struct _idmapper_node{
	int id;
	int current_index_in_dataset;
	struct _idmapper_node* left;
	struct _idmapper_node* right;
}idmapper_node;

int file_contains_variable_number_of_particles=FALSE;

idmapper_node* insert_id(idmapper_node* root, int id, int current_index_in_dataset){
	if(root==NULL){
		root=(idmapper_node*) malloc(sizeof(idmapper_node));
		root->id=id;
		root->current_index_in_dataset=current_index_in_dataset;
		root->left=NULL;
		root->right=NULL;
	}else{ 
		if(id<root->id){
			root->left=insert_id(root->left,id,current_index_in_dataset);
		}else if(id>root->id){
			root->right=insert_id(root->right,id,current_index_in_dataset);
		}else if(id==root->id && id>=0){
			printf("ERROR: id dataset is not unique in id %d, current index in dataset %d\n", id, current_index_in_dataset);
			return NULL;
		}else if(id<0){
			//NOTE: Found negative id in id dataset. This is typically a hint for a file which contains a variable number of particles from a grand canonical simulation.
			file_contains_variable_number_of_particles=TRUE;
		}
	}
	return root;
}

int search_current_index_of_particle_id(idmapper_node* root, int id){
	int current_index_of_particle_id=-1;
	if(root==NULL){
		if(file_contains_variable_number_of_particles==FALSE)
			printf("ERROR: id not found in tree or no correct root provided. id %d is not unique. \n", id);
		return -1;
	}else if(id<root->id){
		current_index_of_particle_id= search_current_index_of_particle_id(root->left,id);
	}else if(id>root->id){
		current_index_of_particle_id= search_current_index_of_particle_id(root->right,id);
	}else{
		current_index_of_particle_id= root->current_index_in_dataset;
	}
	return current_index_of_particle_id;
}

void free_binary_tree_idmapper(idmapper_node* root){
	if(root!=NULL){
		free_binary_tree_idmapper(root->left);
		free_binary_tree_idmapper(root->right);
		free(root);
	}
}

void sort_data_according_to_id_dataset(struct h5md_file* file, int group_number, float* to_be_sorted_data){
	#if defined DEBUG
	h5md_show_hdf5_error_messages();
	#endif
	int i= group_number;
	/////////////////
	//read in id_data of group
	if(file->groups[i].id_dataset_id>=0){
		int data_out_local_id[file->groups[i].natoms_group];
		hid_t dataspace_id_id=H5Dget_space(file->groups[i].id_dataset_id); //Define dataset dataspace (for id_dataset) in file.
		/* 
		* Define hyperslab in the dataset. 
		*/
		int rank_dataset      = H5Sget_simple_extent_ndims(dataspace_id_id);
		hsize_t dataset_slab_offset[rank_dataset];
		dataset_slab_offset[0]=file->current_time;
		dataset_slab_offset[1]=0;
		if(rank_dataset==3)
			dataset_slab_offset[2]=0;

		hsize_t dataset_slab_count[rank_dataset];
		dataset_slab_count[0] = 1;
		dataset_slab_count[1] = file->groups[i].natoms_group;
		if(rank_dataset==3)
			dataset_slab_count[2]=1;
		H5Sselect_hyperslab(dataspace_id_id, H5S_SELECT_SET, dataset_slab_offset, NULL, dataset_slab_count, NULL);
		/*
		* Define memory dataspace.
		*/
		int rank=1; //linear data representation
		hsize_t dimsm[rank];
		dimsm[0]=file->groups[i].natoms_group;

		hid_t memspace_id = H5Screate_simple(rank,dimsm,NULL);

		/* 
		* Define memory hyperslab. 
		*/
		hsize_t offset_out[rank];
		hsize_t count_out[rank];
		offset_out[0]=0;
		count_out[0]=file->groups[i].natoms_group;	
		H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
		/*
		* Read data from hyperslab in the file into the hyperslab in memory
		*/
		hid_t wanted_memory_datatype = H5T_NATIVE_INT;
		H5Dread (file->groups[i].id_dataset_id, wanted_memory_datatype, memspace_id, dataspace_id_id, H5P_DEFAULT, data_out_local_id); 
		H5Sclose(dataspace_id_id);
		//use id data to sort particle positions (use binary tree)
		//create binary tree
		idmapper_node* root=NULL;
		for(int particle_i=0;particle_i<file->groups[i].natoms_group;particle_i++){
			if(particle_i==0){
				//save root of tree
				root=insert_id(root,data_out_local_id[particle_i],particle_i);
			}else{
				insert_id(root,data_out_local_id[particle_i],particle_i);
			}
		}
		//sort data_out_local_pos using the binary tree	
		float _data_out_local_pos_sorted[3*file->groups[i].natoms_group];
		for(int particle_id=0;particle_id<file->groups[i].natoms_group;particle_id++){
			int current_index_of_particle_id=search_current_index_of_particle_id(root,particle_id);
/*				printf("particle with id %d has current index %d at current time %d at x position %f\n", particle_id, current_index_of_particle_id, file->current_time, to_be_sorted_data[3*current_index_of_particle_id+0]);*/
			if(current_index_of_particle_id>=0){
				_data_out_local_pos_sorted[3*particle_id+0]=to_be_sorted_data[3*current_index_of_particle_id+0];
				_data_out_local_pos_sorted[3*particle_id+1]=to_be_sorted_data[3*current_index_of_particle_id+1];
				_data_out_local_pos_sorted[3*particle_id+2]=to_be_sorted_data[3*current_index_of_particle_id+2];
			}else{
				_data_out_local_pos_sorted[3*particle_id+0]=0;
				_data_out_local_pos_sorted[3*particle_id+1]=0;
				_data_out_local_pos_sorted[3*particle_id+2]=0;
			}
		}
		//write sorted data back to to_be_sorted_data
		memcpy(to_be_sorted_data,_data_out_local_pos_sorted,sizeof(float)*3*file->groups[i].natoms_group);
		free_binary_tree_idmapper(root);
	}
	#if defined DEBUG
	h5md_hide_hdf5_error_messages();
	#endif
}

int h5md_sort_data_according_to_id_datasets(struct h5md_file* file, float* to_be_sorted_data){
	int previous_atoms=0;
	for(int i=0; i<file->ngroups; i++){//go through all groups and sort particles according to local ids
		float* local_to_be_sorted_data=&(to_be_sorted_data[3*previous_atoms]); //times 3 since we deal here with flattened position dataset, currently in the function sort_data_according_to_id_dataset
		sort_data_according_to_id_dataset(file,i,local_to_be_sorted_data);
		previous_atoms+=file->groups[i].natoms_group;
	}
	return 0;
}


// reads the next timestep
int h5md_get_timestep(struct h5md_file* file, float *coords){
	#if defined DEBUG
	h5md_show_hdf5_error_messages();
	#endif
	int previous_atoms=0;
	for(int i=0; i<file->ngroups; i++){//go through all groups
		/////////////////
		//read in positions of group
		float data_out_local_pos[file->groups[i].natoms_group*file->groups[i].nspacedims];
		hid_t dataspace_pos_id=H5Dget_space(file->groups[i].pos_dataset_id); //Define dataset dataspace (for pos_dataset) in file.

		/* 
		* Define hyperslab in the dataset. 
		*/
		if(dataspace_pos_id<0)
			continue;

		int rank_dataset      = file->groups[i].nspacedims;
		hsize_t dataset_slab_offset[rank_dataset];
		dataset_slab_offset[0] = file->current_time;
		dataset_slab_offset[1] = 0;
		dataset_slab_offset[2] = 0;

		hsize_t dataset_slab_count[rank_dataset];
		dataset_slab_count[0] = 1;
		dataset_slab_count[1] = file->groups[i].natoms_group;
		dataset_slab_count[2] = file->groups[i].nspacedims;
		H5Sselect_hyperslab(dataspace_pos_id, H5S_SELECT_SET, dataset_slab_offset, NULL, dataset_slab_count, NULL);

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
		offset_out[0]=0;
		count_out[0]=file->groups[i].natoms_group*file->groups[i].nspacedims;	
		H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);


		/*
		* Read data from hyperslab in the file into the hyperslab in memory
		*/
		hid_t wanted_memory_datatype = H5T_NATIVE_FLOAT;
		H5Dread (file->groups[i].pos_dataset_id, wanted_memory_datatype, memspace_id, dataspace_pos_id, H5P_DEFAULT, data_out_local_pos); 
		H5Sclose(memspace_id); //close resources
		H5Sclose(dataspace_pos_id);

		//memcpy to data_out
		memcpy(&(coords[3*previous_atoms]), data_out_local_pos,sizeof(float)*3*file->groups[i].natoms_group);
		
		previous_atoms+=file->groups[i].natoms_group;

	}
	#if defined DEBUG
	h5md_hide_hdf5_error_messages();
	#endif
	return 0;
}

//this function unfolds the *unsorted* position data
int h5md_unfold_positions(struct h5md_file* file, float* unsorted_folded_pos){
	#if defined DEBUG
	h5md_show_hdf5_error_messages();
	#endif
	int previous_atoms=0;
	for(int i=0; i<file->ngroups; i++){//go through all groups
		/////////////////
		//read in image_data of group
		int data_out_local_image[file->groups[i].natoms_group*file->groups[i].nspacedims];
		if(file->groups[i].image_dataset_id>=0){
			hid_t dataspace_image_id=H5Dget_space(file->groups[i].image_dataset_id); //Define dataset dataspace (for image_dataset) in file.
			if(dataspace_image_id<0)
				continue;
			/* 
			* Define hyperslab in the dataset. 
			*/	
			int rank_dataset_images      = file->groups[i].nspacedims;
			hsize_t dataset_slab_offset_images[rank_dataset_images];
			dataset_slab_offset_images[0] = file->current_time;
			dataset_slab_offset_images[1] = 0;
			dataset_slab_offset_images[2] = 0;

			hsize_t dataset_slab_count_images[rank_dataset_images];
			dataset_slab_count_images[0] = 1;
			dataset_slab_count_images[1] = file->groups[i].natoms_group;
			dataset_slab_count_images[2] = file->groups[i].nspacedims;
			H5Sselect_hyperslab(dataspace_image_id, H5S_SELECT_SET, dataset_slab_offset_images, NULL, dataset_slab_count_images, NULL);

			/*
			* Define memory dataspace.
			*/
			int rank_images=1; //linear data representation
			hsize_t dims_images[rank_images];
			dims_images[0]=file->natoms* file->groups[i].nspacedims;

			hid_t memspace_id_images = H5Screate_simple(rank_images,dims_images,NULL);

			/* 
			* Define memory hyperslab. 
			*/
			hsize_t offset_out_images[rank_images];
			hsize_t count_out_images[rank_images];
			offset_out_images[0]=0;
			count_out_images[0]=file->groups[i].natoms_group*file->groups[i].nspacedims;	
			H5Sselect_hyperslab(memspace_id_images, H5S_SELECT_SET, offset_out_images, NULL, count_out_images, NULL);

			/*
			* Read data from hyperslab in the file into the hyperslab in memory
			*/
			hid_t wanted_memory_datatype_images = H5T_NATIVE_INT;
			H5Dread(file->groups[i].image_dataset_id, wanted_memory_datatype_images, memspace_id_images, dataspace_image_id, H5P_DEFAULT, data_out_local_image); 
			H5Sclose(memspace_id_images); //close resources
			H5Sclose(dataspace_image_id);
			//////////////
			//get box_vectors, not using box_information form get_box, since projection of 3 box vectors (9 variables) to only 3 angles and 3 length has inheritent loss of information (only 6 variables)
			float vector_a[3];
			float vector_b[3];
			float vector_c[3];
			int status_box_vectors=get_box_vectors(file, i, file->current_time, vector_a,vector_b,vector_c);
			//use image_data to calculate absolute positions
			float* data_out_local_pos=&(unsorted_folded_pos[3*previous_atoms]);
			if(status_box_vectors==0) {
				for(int j=0; j< file->groups[i].nspacedims*file->groups[i].natoms_group;j=j+3){
					data_out_local_pos[j]=data_out_local_pos[j]+vector_a[0]*data_out_local_image[j]+vector_b[0]*data_out_local_image[j+1]+vector_c[0]*data_out_local_image[j+2];
					data_out_local_pos[j+1]=data_out_local_pos[j+1]+vector_a[1]*data_out_local_image[j]+vector_b[1]*data_out_local_image[j+1]+vector_c[1]*data_out_local_image[j+2];
					data_out_local_pos[j+2]=data_out_local_pos[j+2]+vector_a[2]*data_out_local_image[j]+vector_b[2]*data_out_local_image[j+1]+vector_c[2]*data_out_local_image[j+2];
				}
			}
		}
	}
	#if defined DEBUG
	h5md_hide_hdf5_error_messages();
	#endif
	return 0;
}

//internally reads box_vectors of group_i at time_i into memory. If the dataset is not timedependent, then time_i is ignored.
//memory for the box vectors has to be allocated in advance, it is implicitly assumed that the box is three dimensional.
int get_box_vectors(struct h5md_file* file,  int group_i, int time_i, float* vector_a, float* vector_b, float* vector_c){
	int status;
	//check whether box_dataset is timedependent, if it is timedependent use it, otherwise copy the box information ntime times
	//try to open time-dependent box dataset, get box_dataset_timedependent_id
	char* full_path_box_dataset_timedependent=concatenate_strings((const char*) file->groups[group_i].group_path,(const char*) "/box/edges/value");	
	hid_t box_timedependent_dataset_id=H5Dopen2(file->file_id, full_path_box_dataset_timedependent ,H5P_DEFAULT);

	//try to open time-independent box dataset, get box_dataset_timeindependent_id
	char* full_path_box_dataset_timeindependent=concatenate_strings((const char*) file->groups[group_i].group_path,(const char*) "/box/edges");
	hid_t box_timeindependent_dataset_id=H5Dopen2(file->file_id, full_path_box_dataset_timeindependent ,H5P_DEFAULT);


	if(box_timedependent_dataset_id>=0){
		//timedependent dataset exists, use it
		//read timedependent dataset 
		float *data_box = malloc(file->ntime*3*3*sizeof(float));
		H5Dread(box_timedependent_dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_box);
		vector_a[0]=data_box[time_i*9+0];
		vector_a[1]=data_box[time_i*9+1];
		vector_a[2]=data_box[time_i*9+2];
		vector_b[0]=data_box[time_i*9+3];
		vector_b[1]=data_box[time_i*9+4];
		vector_b[2]=data_box[time_i*9+5];
		vector_c[0]=data_box[time_i*9+6];
		vector_c[1]=data_box[time_i*9+7];
		vector_c[2]=data_box[time_i*9+8];
		free (data_box);
		status=0;
	}else if(box_timeindependent_dataset_id>=0){
		//read timeindependent dataset
		//decided whether box is cubic (dataset contains a vector) or triclinic (dataset contains a matrix)
		int dims_box;
		int is_cubic=h5md_get_length_of_one_dimensional_dataset(file, full_path_box_dataset_timeindependent, &dims_box);
		float* data_box;
		H5T_class_t box_class_out;
		h5md_read_timeindependent_dataset_automatically(file, full_path_box_dataset_timeindependent,(void**) &data_box, &box_class_out);
		if(is_cubic==0){
			//box is cubic implies vector
			vector_a[0]=data_box[0];
			vector_a[1]=0;
			vector_a[2]=0;
			vector_b[0]=0;
			vector_b[1]=data_box[1];
			vector_b[2]=0;
			vector_c[0]=0;
			vector_c[1]=0;
			vector_c[2]=data_box[2];
			status=0;
		}else{
			//box is triclinic implies matrix
			//VMD expects system to be 3dimensional -> assume 3x3 matrix
			vector_a[0]=data_box[0];
			vector_a[1]=data_box[1];
			vector_a[2]=data_box[2];
			vector_b[0]=data_box[3];
			vector_b[1]=data_box[4];
			vector_b[2]=data_box[5];
			vector_c[0]=data_box[6];
			vector_c[1]=data_box[7];
			vector_c[2]=data_box[8];
			status=0;
		}
	}else{
		//printf("No box information found\n");
		vector_a[0]=0;
		vector_a[1]=0;
		vector_a[2]=0;
		vector_b[0]=0;
		vector_b[1]=0;
		vector_b[2]=0;
		vector_c[0]=0;
		vector_c[1]=0;
		vector_c[2]=0;
		status=-1;
	}
	free(full_path_box_dataset_timeindependent);
	free(full_path_box_dataset_timedependent);
	return status;

}

// internally reads the box information of a given group into the memory
int get_box_information(struct h5md_file* file, int group_number, int time_i, h5md_box* out_box){
	int status;
	h5md_box* box=out_box;

	//check whether box_dataset is timedependent, if it is timedependent use it, otherwise copy the box information ntime times
	//try to open time-dependent box dataset, get box_dataset_timedependent_id
	char* full_path_box_dataset_timedependent=concatenate_strings((const char*) file->groups[group_number].group_path,(const char*) "/box/edges/value");
	hid_t box_timedependent_dataset_id=H5Dopen2(file->file_id, full_path_box_dataset_timedependent ,H5P_DEFAULT);
	free(full_path_box_dataset_timedependent);

	//try to open time-independent box dataset, get box_dataset_timeindependent_id
	char* full_path_box_dataset_timeindependent=concatenate_strings((const char*) file->groups[group_number].group_path,(const char*) "/box/edges");
	hid_t box_timeindependent_dataset_id=H5Dopen2(file->file_id, full_path_box_dataset_timeindependent ,H5P_DEFAULT);
	free(full_path_box_dataset_timeindependent);

	float vector_a[3];
	float vector_b[3];
	float vector_c[3];
	if(box_timedependent_dataset_id>=0){
		//timedependent dataset exists, use it
		//read timedependent dataset 
		status=get_box_vectors(file, group_number, time_i, vector_a,vector_b,vector_c);
		//according to VMD's molfile_timestep_t documentation
		//process to angles and lengths
		box->A=calculate_length_of_vector(vector_a,3);
		box->B=calculate_length_of_vector(vector_b,3);
		box->C=calculate_length_of_vector(vector_c,3);
		box->alpha= calculate_angle_between_vectors(vector_b,vector_c,3);
		box->beta= calculate_angle_between_vectors(vector_a,vector_c,3);
		box->gamma= calculate_angle_between_vectors(vector_a,vector_b,3);
	}else{
		if(box_timeindependent_dataset_id>=0){
			time_i=FALSE;
			//read timeindependent dataset
			//decided whether box is cubic (dataset contains a vector) or triclinic (dataset contains a matrix)
			status=get_box_vectors(file, group_number, time_i, vector_a,vector_b,vector_c);
			//according to VMD's molfile_timestep_t documentation
			//process to angles and lengths
			box->A=calculate_length_of_vector(vector_a,3);
			box->B=calculate_length_of_vector(vector_b,3);
			box->C=calculate_length_of_vector(vector_c,3);
			box->alpha=calculate_angle_between_vectors(vector_b,vector_c,3);
			box->beta=calculate_angle_between_vectors(vector_a,vector_c,3);
			box->gamma=calculate_angle_between_vectors(vector_a,vector_b,3);
		}else{
			status=-1;
			printf("No box information found\n");
		}
	}


	return status;

}


int h5md_get_box_information(struct h5md_file* file, float* out_box_information){
	//only use information of first group (file->groups[0]...) since VMD only supports one box
	int group_number=0;
	h5md_box box;
	int status=get_box_information(file, group_number, file->current_time, &box);
	out_box_information[0]=box.A;
	out_box_information[1]=box.B;
	out_box_information[2]=box.C;
	out_box_information[3]=box.alpha;
	out_box_information[4]=box.beta;
	out_box_information[5]=box.gamma;
	return status;
}



//read timeindependent dataset automatically
int h5md_read_timeindependent_dataset_automatically(struct h5md_file* file, char* dataset_name, void** _data_out, H5T_class_t* type_class_out){
	#if defined DEBUG
	h5md_show_hdf5_error_messages();
	#endif
	int status=-1;
	hid_t dataset_id = H5Dopen2(file->file_id, dataset_name,H5P_DEFAULT);
	if(dataset_id<0){
		//printf("Dataset %s could not be opened.\n", dataset_name);
		return status;
	}
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
		int status_read=H5Dread(dataset_id, wanted_memory_datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out);
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

	//close resources
	H5Dclose (dataset_id);
	H5Sclose (dataspace_id);
	H5Tclose (datatype);
	#if defined DEBUG
	h5md_hide_hdf5_error_messages();
	#endif
	return status;
}

int h5md_get_all_infromation_about_property(struct h5md_file *file, char* property, void** infromation_out){
	#if defined DEBUG
	h5md_show_hdf5_error_messages();
	#endif
	int status=0;
	float* data_out= malloc(sizeof(float)*file->natoms); //allocate space for data in memory, which have the order data_out[atom_nr]

	int previous_atoms=0;
	for(int i=0; i<file->ngroups; i++){//go through all groups
		hid_t dataset_id;
		hid_t dataspace_id;
		if(strcmp(property,"mass")==0){
			dataset_id=file->groups[i].mass_dataset_id;
		}
		if(strcmp(property,"charge")==0){
			dataset_id=file->groups[i].charge_dataset_id;
		}
		if(strcmp(property,"species")==0){
			dataset_id=file->groups[i].species_dataset_id;
		}
		if(dataset_id<0){
			status=-1;
			return status;
		}

		dataspace_id=H5Dget_space(dataset_id); //Define dataset dataspace in file.
		/* 
		* Define hyperslab in the dataset. 
		*/
		int rank_dataset      = H5Sget_simple_extent_ndims(dataspace_id);
		hsize_t dataset_slab_offset[rank_dataset];
		dataset_slab_offset[0] = 0;
		hsize_t dataset_slab_count[rank_dataset];
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
		hid_t wanted_memory_datatype = H5T_NATIVE_FLOAT;
		int current_status=H5Dread(dataset_id, wanted_memory_datatype, memspace_id, dataspace_id, H5P_DEFAULT, data_out);
		if(current_status<0)
			status=current_status;
		H5Sclose(memspace_id); //close resources
		H5Sclose(dataspace_id);
	}

	*infromation_out=(void*) data_out;
	#if defined DEBUG
	h5md_hide_hdf5_error_messages();
	#endif
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

	if( (status_read!=0) )
		return -1;
	else
		return 0;
}

int h5md_get_length_of_one_dimensional_dataset(struct h5md_file *file,char *dataset_name, int *length_of_dataset){
	hid_t dataset_id = H5Dopen2(file->file_id, dataset_name,H5P_DEFAULT);
	if(dataset_id<0)
		return -1;
	hid_t dataspace_id = H5Dget_space(dataset_id);    /* dataspace handle */
	hsize_t rank_dataset      = H5Sget_simple_extent_ndims(dataspace_id);
	unsigned long long int dims_dataset[rank_dataset];
	H5Sget_simple_extent_dims(dataspace_id, dims_dataset, NULL);
	if(rank_dataset==1){

		*length_of_dataset=dims_dataset[0];
		return 0;
	}else{
/*		printf("Dataset %s is not one dimensional.\n", dataset_name);	*/
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
	if(file_id != NULL)
		return 0;
	else
		return -1;
}

/* write operations */
//creates a h5md_file iff it does not exist yet
int h5md_create_file(struct h5md_file **_file, const char* filename){
	struct h5md_file *file = (struct h5md_file*) malloc(sizeof(struct h5md_file));
	initialize_h5md_struct(file);
	/* Create a new file using default properties. */
	hid_t file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT); //H5F_ACC_EXCL <-> fail if file already exists, alternative H5F_ACC_TRUNC <-> overwrite file
	file->file_id=file_id;
	h5md_set_author(file, NULL, NULL); //sets author name by default to the user account name, can be overwritten by another call of h5md_set_author(file,"username");
	h5md_set_creator(file, NULL, NULL);
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
	hid_t dataset_id = H5Dopen2(file->file_id, absolute_name_of_dataset, H5P_DEFAULT);
	if(dataset_id>=0){
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
			status=-1;
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
	if(name == NULL){
		char* username=getlogin();
		status_name=H5LTset_attribute_string(file->file_id, "/h5md/author", "name", username);
	}else{
		status_name=H5LTset_attribute_string(file->file_id, "/h5md/author", "name", name);
	}

	int status_email =-1;
	if(email_address==NULL){
		status_email=H5LTset_attribute_string(file->file_id, "/h5md/author", "email", "email was not provided");
	}else{
		status_email=H5LTset_attribute_string(file->file_id, "/h5md/author", "email", email_address);
	}

	if( status_name>=0 && status_email>=0 )
		return 0;
	else
		return -1;
}

int h5md_set_creator(struct h5md_file* file, char* name, char* version){
	hid_t link_crt_plist = H5Pcreate(H5P_LINK_CREATE);
	H5Pset_create_intermediate_group(link_crt_plist, TRUE);	// Set flag for intermediate group creation

	H5Gcreate(file->file_id, "/h5md/creator", link_crt_plist, H5P_DEFAULT, H5P_DEFAULT);

	herr_t status_name =-1;
	if(name == NULL){
		status_name=H5LTset_attribute_string(file->file_id, "/h5md/creator", "name", "libh5md");
	}else{
		status_name=H5LTset_attribute_string(file->file_id, "/h5md/creator", "name", name);
	}

	int status_version =-1;
	if(version==NULL){
		status_version=H5LTset_attribute_string(file->file_id, "/h5md/creator", "version", "version not provided");
	}else{
		status_version=H5LTset_attribute_string(file->file_id, "/h5md/creator", "version", version);
	}

	if(status_name>=0  && status_version>=0 )
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
	file->ntime=0;
	file->last_error_message="";
	file->natoms=0;
	return 0;
}

int max(int a, int b){
	if(a>b)
		return a;
	else
		return b;
}

float calculate_length_of_vector(float* vector, int dimensions){
	float length=0.0;
	for(int i=0;i<dimensions;i++){
	length+=vector[i]*vector[i];
	}
	length=sqrt(length);
	return length;
}

//returns 0 if "vector" is the zero vector, otherwise nonzero
int is_zero_vector(float* vector, int dimensions){
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

float calculate_angle_between_vectors(float* vector1, float* vector2, int dimensions){
	if(is_zero_vector(vector1, dimensions)==0 || is_zero_vector(vector2, dimensions)==0){
		printf("ERROR: Cannot calculate angle to the zero vector which has length 0\n");
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

char* mystrdup(char* str) {
	char* new_str=(char*) malloc((strlen(str)+1)*sizeof(char));
	strcpy(new_str,str);
	return new_str;
}
