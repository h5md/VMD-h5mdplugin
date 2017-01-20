extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "../libh5md.h"
}//extern "C" is needed, since we include C files into a C++ program


int test_create_file(){
	//h5md_show_hdf5_error_messages();
	char* filename="./samples/test_write.h5";
	int status_del= h5md_delete_file(filename); //delete file if exists, otherwise error when trying to write
	
	struct h5md_file *file;
	int status=h5md_create_file(&file,filename); //creates h5 file if it does not exist yet, if it already exists creation fails
	h5md_close(file);
	return status;
}


int test_write_file_with_datasets(){
	test_create_file();
	struct h5md_file* file;
 	int status_open = h5md_open(&file, "./samples/test_write.h5", 1);
	float data_in1[4]={1.9,1,1,1};
	hid_t datatype1= H5T_NATIVE_FLOAT;
	int rank_in1=1;
	hsize_t dims_in1[1] = {4};
	h5md_write_dataset(file, "/dataset/go/test" , datatype1, data_in1, rank_in1, dims_in1);
	char* data_in2[4]={"bla","bla1","bla2","bla3"};
   	hid_t vls_type_c_id = H5Tcopy(H5T_C_S1);
     	H5Tset_size(vls_type_c_id, H5T_VARIABLE) ;
	hid_t datatype2 = vls_type_c_id;
	int rank_in2=1;
	hsize_t dims_in2[1] = {4};
	int status= h5md_write_dataset(file, "/dataset/go/jjjj" , datatype2, (void*) data_in2, rank_in2, dims_in2);
	h5md_close(file);
	return status;

}

//int test_write_fixed_length_string_dataset(){
//	struct h5md_file* file;
//	int status;
// 	h5md_open(&file, "./samples/test_write.h5", 1);
//	//writing fixed length strings not yet working, probably better to give h5md_write another argument string_length (see http://www.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/C/H5T/h5ex_t_string.c) with the following properties, e.g. like
//	//string_length= 0 iff data_in are not strings
//	//string_length >0 <-> strings have length string_length
//	//string_length <0 <-> strings have length H5T_VARIABLE
//	char* data_in3[4]={"xxxxx","zzz1","ggg2","sss3"};
//	int fixed_length_for_string=16;
//	hid_t filetype3 = H5Tcopy (H5T_C_S1);
//	H5Tset_size (filetype3, fixed_length_for_string+1);	//strings are 0 terminated in C
//	hid_t memtype3 = H5Tcopy (H5T_C_S1);
//	H5Tset_size (memtype3, fixed_length_for_string+1);
//	int rank_in3=1;
//	hsize_t dims_in3[1] = {4};
//	status=h5md_write_dataset(file, "/fixed_length_string" , memtype3, data_in3, rank_in3, dims_in3);
//	h5md_close(file);
//	return status;
//}

int read_string_dataset(){
	test_write_file_with_datasets();
	struct h5md_file* file;
 	int status= h5md_open(&file, "./samples/test_write.h5", -1);
	H5T_class_t type_class_out;
	void* data_out;
	char* dataset_name= "/dataset/go/jjjj";
	int status_read=h5md_read_timeindependent_dataset_automatically(file, dataset_name, &data_out, &type_class_out);

	//display read_data, decide about type at runtime
	if(status_read==0){
		switch(type_class_out){
		case H5T_INTEGER:
		{	
			int* data_out_corr= (int *) data_out;
			printf("read data int test: %d\n",data_out_corr[0]);	
		break;
		}	
		case H5T_FLOAT:
		{
			float* data_out_corr= (float*) data_out;
			printf("read data float test: %f\n",data_out_corr[1]);
		break;
		}	
		case H5T_STRING:
		{	
			char** data_out_corr= (char**) data_out;
			printf("read data string test: %s\n",data_out_corr[2]);
		break;
		}
		default:
			printf("datatype not implemented\n");
		break;
		}
	}
	int length_of_dataset;
	h5md_get_length_of_one_dimensional_dataset(file, dataset_name, &length_of_dataset);
	h5md_free_timeindependent_dataset_automatically(type_class_out,data_out, length_of_dataset); // free stuff
	h5md_close(file);
	return status_read;
}

int read_real_dataset_float(){
	struct h5md_file* file_real;
 	h5md_open(&file_real, "../samples/full_vmd_structure.h5", -1);
	H5T_class_t type_class_out;
	void* data_out;
	char* dataset_name= "particles/atoms/mass";
	int status_read=h5md_read_timeindependent_dataset_automatically(file_real, dataset_name, &data_out, &type_class_out);
	//display read_data, decide about type at runtime
	if(status_read==0){
		switch(type_class_out){
		case H5T_INTEGER:
		{	
			int* data_out_corr= (int *) data_out;
			printf("read data int test: %d\n",data_out_corr[0]);	
		break;
		}	
		case H5T_FLOAT:
		{
			float* data_out_corr= (float*) data_out;
			printf("read data float test: %f\n",data_out_corr[1]);
		break;
		}	
		case H5T_STRING:
		{	
			char** data_out_corr= (char**) data_out;
			printf("read data string test: %s\n",data_out_corr[2]);
		break;
		}
		default:
			printf("datatype not implemented\n");
		break;
		}
	}
	int length_of_dataset;
	h5md_get_length_of_one_dimensional_dataset(file_real, dataset_name, &length_of_dataset);
	h5md_free_timeindependent_dataset_automatically(type_class_out,data_out, length_of_dataset); // free stuff
	h5md_close(file_real);
	return status_read;
}


int read_real_dataset_string(){
	struct h5md_file* file_real;
 	h5md_open(&file_real, "../samples/full_vmd_structure.h5", -1);
	H5T_class_t type_class_out;
	void* data_out;
	char* dataset_name= "/parameters/vmd_structure/name";
	int status_read=h5md_read_timeindependent_dataset_automatically(file_real, dataset_name, &data_out, &type_class_out);
	//display read_data, decide about type at runtime
	if(status_read==0){
		switch(type_class_out){
		case H5T_INTEGER:
		{	
			int* data_out_corr= (int *) data_out;
			printf("read data int test: %d\n",data_out_corr[0]);	
		break;
		}	
		case H5T_FLOAT:
		{
			float* data_out_corr= (float*) data_out;
			printf("read data float test: %f\n",data_out_corr[1]);
		break;
		}	
		case H5T_STRING:
		{	
			char** data_out_corr= (char**) data_out;
			printf("read data string test: %s\n",data_out_corr[2]);
		break;
		}
		default:
			printf("datatype not implemented\n");
		break;
		}
	}
	int length_of_dataset;
	h5md_get_length_of_one_dimensional_dataset(file_real, dataset_name, &length_of_dataset);
	h5md_free_timeindependent_dataset_automatically(type_class_out,data_out, length_of_dataset); // free stuff
	h5md_close(file_real);
	return status_read;
}

int test_delete_dataset(char* absolute_name_of_dataset){//depends on test_write_file_with_datasets() since it tries to access this file

	test_write_file_with_datasets();
	struct h5md_file* file;
 	h5md_open(&file, "./samples/test_write.h5", 1);

	int status=h5md_delete_dataset(file, absolute_name_of_dataset);
	printf("status_del dataset %d \n", status);
	h5md_close(file);
	return status;
}

int test_appending_to_dataset_position_style(){//depends on test_create_file() since it tries to access this file
	test_create_file();
	struct h5md_file* file;
 	h5md_open(&file, "./samples/test_write.h5", 1);
	//Test appending to dataset (position style)
	int rank_in3=3; //should be 3 for positions : time_i, atom_i, x_i
	hsize_t dims_in3[3] = {1, 3, 3};
    	float data_in3[3][3] = {	{1.9, 1, 1},       /* data to write */
					{1, 1.0, 1},
					{1, 1, 1.6} };
	hid_t datatype3= H5T_NATIVE_FLOAT;

	h5md_append_dataset(file, "testappend", datatype3, (void*) data_in3, rank_in3, dims_in3);


	int rank_in4=3; //should be 3 for positions : time_i, atom_i, x_i
	hsize_t dims_in4[3] = {1, 3, 3};
    	float data_in4[3][3] = {	{5.9, 5, 5},       /* data to write */
					{5, 6.0, 2},
					{5, 5, 5.6} };
	hid_t datatype4= H5T_NATIVE_FLOAT;
	int status=h5md_append_dataset(file, "testappend", datatype4, (void*) data_in4, rank_in4, dims_in4);
	h5md_close(file);
	return status;
}

int test_appending_to_dataset_position_style_grandcanonical(){
	test_create_file();
	struct h5md_file* file;
 	h5md_open(&file, "./samples/test_write.h5", 1);
	//Test appending to dataset (position style)
	int rank_in3=3; //should be 3 for positions : time_i, atom_i, x_i
	hsize_t dims_in3[3] = {1, 3, 3};
    	float data_in3[3][3] = {	{1.9, 1, 1},       /* data to write */
					{1, 1.0, 1},
					{1, 1, 1.6} };
	hid_t datatype3= H5T_NATIVE_FLOAT;

	h5md_append_dataset(file, "/particles/atoms/position/value", datatype3, (void*) data_in3, rank_in3, dims_in3);


	int rank_in4=3; //should be 3 for position data: time_i, atom_i, x_i
	hsize_t dims_in4[3] = {1, 4, 3};            /* data1 dimensions */
    	float data_in4[4][3] = {	{5.9, 5, 5},       /* data to write */
					{5, 5.0, 5},
					{5, 5.0, 5},
					{5, 5.0, 5} };
	hid_t datatype4= H5T_NATIVE_FLOAT;
	int status=h5md_append_dataset(file, "/particles/atoms/position/value", datatype4, (void*) data_in4, rank_in4, dims_in4);
	h5md_close(file);
	return status;
}

int test_get_fill_value_float(){
	test_appending_to_dataset_position_style_grandcanonical();
	struct h5md_file* file;
 	h5md_open(&file, "./samples/test_write.h5", 1);	
	
	char* dset_name="/particles/atoms/position/value";
	float filler;
	get_fill_value(file, dset_name, (void*) &filler );
	printf("fill value is %f\n", filler);
	return filler;
}


int test_appending_to_onedimensional_dataset_float(){
	struct h5md_file* file;
 	h5md_open(&file, "./samples/test_write.h5", 1);

	//Test appending to dataset (one dimensional array of floats)
	int rank_in5=1;
	hsize_t dims_in5[1] = {1};
    	float data_in5[1] = {999.8};
	hid_t datatype5= H5T_NATIVE_FLOAT;

	int status=h5md_append_dataset(file, "/dataset/go/test", datatype5, (void*) data_in5, rank_in5, dims_in5);
	h5md_close(file);
	return status;

}


int test_appending_to_onedimensional_string_dataset(){
	struct h5md_file* file;
 	h5md_open(&file, "./samples/test_write.h5", 1);
	//Test appending to dataset (one dimensional array of strings)
	int rank_in6=1;
	hsize_t dims_in6[1] = {1};
    	char* data_in6[1] = {"appending"};
	hid_t vls_type_c_id6 = H5Tcopy(H5T_C_S1);
     	H5Tset_size(vls_type_c_id6, H5T_VARIABLE) ;
	hid_t datatype6 = vls_type_c_id6;

	int status=h5md_append_dataset(file, "/dataset/go/jjjj", datatype6, (void*) data_in6, rank_in6, dims_in6);
	h5md_close(file);
	return status;

}


int test_h5md_get_timestep(){
	struct h5md_file* file;
 	h5md_open(&file, "../samples/full_vmd_structure.h5", -1);
	h5md_seek_timestep(file, 1);
	float* coords= (float*) malloc(3*sizeof(float));
	h5md_get_timestep(file, coords);
	int status;
	if(coords[3*2+0]-0.109156<0.0001)
		status= 0;
	else
		status=  -1;
	return status; 
}
