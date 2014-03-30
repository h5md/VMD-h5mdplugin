#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "libh5md.h"


int h5md_test_creation(){
	h5md_show_hdf5_error_messages();
	int status_del= h5md_delete_file("./samples/test_write.h5"); //delete if exists, otherwise error 
	char* filename="./samples/test_write.h5";
	struct h5md_file *file;
	h5md_create_file(&file,filename); //creates h5 file if it does not exist yet

	// #######test1 ######## int, float

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
	h5md_write_dataset(file, "/dataset/go/jjjj" , datatype2, data_in2, rank_in2, dims_in2);

	/* //writing fixed length strings not yet working, probably better to give h5md_write another argument string_length (see http://www.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/C/H5T/h5ex_t_string.c) with the following properties, e.g. like
	//string_length= 0 iff data_in are not strings
	//string_length >0 <-> strings have length string_length
	//string_length <0 <-> strings have length H5T_VARIABLE
	char* data_in3[4]={"xxxxx","zzz1","ggg2","sss3"};
	int fixed_length_for_string=8;
	hid_t filetype3 = H5Tcopy (H5T_FORTRAN_S1);
	H5Tset_size (filetype3, fixed_length_for_string - 1);
	hid_t memtype3 = H5Tcopy (H5T_C_S1);
	H5Tset_size (memtype3, fixed_length_for_string);
	int rank_in3=1;
	hsize_t dims_in3[1] = {4};
	h5md_write_dataset(file, "/dataset/go/kkkk" , memtype3, data_in3, rank_in3, dims_in3);
	*/


	h5md_close(file);




	//read written dataset
	H5T_class_t type_class_out;
	void* data_out;

	//int status_read=h5md_read_timeindependent_dataset_automatically(file, "/dataset/go/jjjj", &data_out, &type_class_out); //test 1


	//########other tests ##########
	
	struct h5md_file* file_real;
 	h5md_open(&file_real, "./samples/full_vmd_structure.h5");
	//int status_read=h5md_read_timeindependent_dataset_automatically(file_real, "/parameters/vmd_structure/mass", &data_out, &type_class_out);  //test 2
	int status_read=h5md_read_timeindependent_dataset_automatically(file_real, "/parameters/vmd_structure/name", &data_out, &type_class_out);	//test 3

	h5md_close(file_real);

	int test_delete=-1;	//test 4
	if(test_delete>0){
		printf("delete\n");
		//delete dataset again
		h5md_delete_dataset(file, "/dataset/go/jjjj");
	}


	//########################################

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
			printf("read data float test: %f\n",data_out_corr[0]);
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

		// free stuff
		h5md_free_timeindependent_dataset_automatically(type_class_out,data_out);
	}


}

int main(){
	h5md_test_creation();
}
