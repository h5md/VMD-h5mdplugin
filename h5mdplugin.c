/***************************************************************************
 *
 *            (C) Copyright 2013
 *
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: xyzplugin.c,v $
 *      $Author: Jonas Landsgesell, Sascha Ehrhardt S$
 *      $Revision: 1.0 $       $Date: 2013/08/01 14:00:00 $
 *
 ***************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "molfile_plugin.h"
#include "hdf5.h"
#include "hdf5_hl.h"

//element symbols
static const char *element_symbols[] = { 
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P" , "S",  "Cl", "Ar", "K",  "Ca", "Sc",
    "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", 
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
    "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
    "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg"
};

const char default_name[16]="O";
const char default_type[16]="O";
const char default_resname[8]="";
const int default_resid= 1;
const char default_segid[8]= "";
const char default_chain[2]= "";
const char default_altloc[2]= "";
const char default_insertion[2]= "";
const const float default_occupancy= 1.0;
const float default_bfactor= 1.0; 
const float default_mass= 1.0;
const float default_charge= 0.0;
const float default_radius= 0.5;
const int default_atomicnumber= 1;

typedef struct {
	int nspacedims;
	molfile_atom_t *atomlist;
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
} h5mddata;

//global non constant variables
int position_already_read =-1;
int reads = 0;
int current_time = 0;

static void *open_h5md_read(const char *filename, const char *filetype, int *natoms) {
	h5mddata *data;
	data = (h5mddata *) malloc(sizeof(h5mddata));

	hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	hid_t dataset_id = H5Dopen(file_id, "/particles/atoms/position/value",H5P_DEFAULT);

	//get dims
	hid_t dataspace = H5Dget_space(dataset_id);
	unsigned long long int dims_out[3];
	H5Sget_simple_extent_dims(dataspace, dims_out, NULL );

	//save repeatedly needed files in h5mddata struct 
	data->file_name = filename;
	data->ntime = dims_out[0];
	data->natoms = dims_out[1];
	data->file_id = file_id;
	*natoms = data->natoms;
	data->nspacedims = dims_out[2];
	data->dataset_id=dataset_id;
	return data;

}

static int read_h5md_structure(void *mydata, int *optflags,molfile_atom_t *atoms) {
	molfile_atom_t *atom;
	h5mddata *data = (h5mddata *) mydata;
	*optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS;
	int error_status_of_read=-1;
	//check whether "/parameters/vmd_structure/indexOfSpecies" exists. if dataset exists, then read vmd_structure, else do not try to read vmd_structure.
	printf("check whether /parameters/vmd_structure/indexOfSpecies exists\n\n");
	hid_t dataset_id_indexOfSpecies = H5Dopen(data->file_id, "/parameters/vmd_structure/indexOfSpecies",H5P_DEFAULT);
	if(dataset_id_indexOfSpecies <0){
		error_status_of_read = read_h5md_structure_no_vmd_structure(mydata,optflags,atoms);
	}else{
		error_status_of_read = read_h5md_structure_vmd_structure(mydata,optflags,atoms);
	}
	return error_status_of_read;
}

//find index of species "int species" in data_index_Of_Species which has length len_index
int find_index_of_species(int *data_indexOfSpecies, int species,int len_index){
	int ndx=-1;
	for(int i=0;i<len_index;i++){
		if(species== data_indexOfSpecies[i]){
			ndx=i;
			break;
		}
	}

	return ndx;
}

//count the different species that are present in data_species, requires pointer to array and length of array as argument
int find_first_index_in_array(int *data_species, int index_of_current_value, int length){
	int index_at_which_the_value_of_current_index_first_occured=-1;
	for(int i=0;i<length;i++){
		if(data_species[i]==data_species[index_of_current_value]){
			index_at_which_the_value_of_current_index_first_occured=i;
			break;
		}
	}
	return index_at_which_the_value_of_current_index_first_occured;	
}

int count_different_species(int *data_species, int length){
	int count =0;
	for(int i=0;i<length;i++){
		if(find_first_index_in_array(data_species,i,length)==i){
			count+=1;
		}	
	}
	return count;
}
//end count different species


int read_h5md_structure_vmd_structure(void *mydata, int *optflags,molfile_atom_t *atoms) {
	molfile_atom_t *atom;
	h5mddata *data = (h5mddata *) mydata;

	*optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS;

	//load species
	int data_species[data->natoms];
	hid_t dataset_id_species = H5Dopen(data->file_id, "/particles/atoms/species/value",H5P_DEFAULT);
	herr_t status_species= H5Dread(dataset_id_species, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_species);
	status_species = H5Dclose (dataset_id_species);

	//load indexOfSpecies
	hid_t dataset_id_indexOfSpecies = H5Dopen(data->file_id, "/parameters/vmd_structure/indexOfSpecies",H5P_DEFAULT);
	hid_t dataspace_indexOfSpecies = H5Dget_space(dataset_id_indexOfSpecies);
	unsigned long long int dims_indexOfSpecies[1];
	H5Sget_simple_extent_dims(dataspace_indexOfSpecies, dims_indexOfSpecies, NULL );

	//check consistency of indexOfSpecies and species/value
	if(count_different_species(data_species,data->natoms)!=dims_indexOfSpecies[0]){
		printf("ERROR: /parameters/vmd_structure/indexOfSpecies does not contain as much different species as different species are present in /particles/atoms/species/value !\n");
		return MOLFILE_ERROR;
	}

	//load content of indexOfSpecies
	int data_indexOfSpecies[dims_indexOfSpecies[0]];
	herr_t status_indexOfSpecies= H5Dread(dataset_id_indexOfSpecies, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_indexOfSpecies);
	status_indexOfSpecies = H5Dclose (dataset_id_indexOfSpecies);


	//load name
	//open dataset
	hid_t dataset_id_name = H5Dopen(data->file_id, "/parameters/vmd_structure/name",H5P_DEFAULT);
	//get dataspace
	hid_t dataspace_name = H5Dget_space(dataset_id_name);
	//get dims
	unsigned long long int dims_name[1];
	H5Sget_simple_extent_dims(dataspace_name, dims_name, NULL );


	//load content of name
	//get datatype and its size
	hid_t filetype= H5Dget_type (dataset_id_name);
	size_t sdim = H5Tget_size (filetype);
	sdim++;// Make room for null terminator
	/*
    	* Get dataspace and allocate memory for read buffer.  This is a
     	* two dimensional dataset so the dynamic allocation must be done
        * in steps.
        */
	int ndims = H5Sget_simple_extent_dims (dataspace_name, dims_name, NULL);
	// Allocate array of pointers to rows.
	char **data_name = (char **) malloc (dims_name[0] * sizeof (char *));		
	// Allocate space for integer data.
	data_name[0] = (char *) malloc (dims_name[0] * sdim * sizeof (char));
	// Set the rest of the pointers to rows to the correct addresses.
	for (int i=1; i<dims_name[0]; i++)
		data_name[i] = data_name[0] + i * sdim;
	// Create the memory datatype.
	hid_t memtype = H5Tcopy (H5T_C_S1);
	herr_t status_name = H5Tset_size (memtype, sdim);
	// Read the data.
	status_name = H5Dread (dataset_id_name, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_name[0]);
	//close
	status_name = H5Dclose (dataset_id_name);
	status_name = H5Sclose (dataspace_name);
	status_name = H5Tclose (filetype);
	status_name = H5Tclose (memtype);



	//load type
	//open dataset
	hid_t dataset_id_type = H5Dopen(data->file_id, "/parameters/vmd_structure/type",H5P_DEFAULT);
	//get datapscae
	hid_t dataspace_type = H5Dget_space(dataset_id_type);
	//get dims
	unsigned long long int dims_type[1];
	H5Sget_simple_extent_dims(dataspace_type, dims_type, NULL );

	//load content of type
	//get datatype and its size
	hid_t filetype_type= H5Dget_type (dataset_id_type);
	size_t sdim_type = H5Tget_size (filetype_type);
	sdim_type++;// Make room for null terminator
    	// allocate memory for read buffer.  
	ndims = H5Sget_simple_extent_dims (dataspace_type, dims_type, NULL);
	// Allocate array of pointers to rows.
	char **data_type = (char **) malloc (dims_type[0] * sizeof (char *));		
	// Allocate space for integer data.
	data_type[0] = (char *) malloc (dims_type[0] * sdim_type * sizeof (char));
	// Set the rest of the pointers to rows to the correct addresses.
	for (int i=1; i<dims_type[0]; i++)
		data_type[i] = data_type[0] + i * sdim_type;
	// Create the memory datatype.
	hid_t memtype_type = H5Tcopy (H5T_C_S1);
	herr_t status_type = H5Tset_size (memtype_type, sdim_type);
	// Read the data.
	status_type = H5Dread (dataset_id_type, memtype_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_type[0]);
	//close
	status_type = H5Dclose (dataset_id_type);
	status_type = H5Sclose (dataspace_type);
	status_type = H5Tclose (filetype_type);
	status_type = H5Tclose (memtype_type);


	//load mass
	hid_t dataset_id_mass = H5Dopen(data->file_id, "/parameters/vmd_structure/mass",H5P_DEFAULT);
	hid_t dataspace_mass = H5Dget_space(dataset_id_mass);
	unsigned long long int dims_mass[1];
	H5Sget_simple_extent_dims(dataspace_mass, dims_mass, NULL );
	//load content of mass
	float data_mass[dims_indexOfSpecies[0]];
	herr_t status_mass= H5Dread(dataset_id_mass, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_mass);
	status_mass = H5Dclose (dataset_id_mass);

	//load radius
	hid_t dataset_id_radius = H5Dopen(data->file_id, "/parameters/vmd_structure/radius",H5P_DEFAULT);
	hid_t dataspace_radius = H5Dget_space(dataset_id_radius);
	unsigned long long int dims_radius[1];
	H5Sget_simple_extent_dims(dataspace_radius, dims_radius, NULL );
	//load content of radius
	float data_radius[dims_indexOfSpecies[0]];
	herr_t status_radius= H5Dread(dataset_id_radius, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_radius);
	status_radius = H5Dclose (dataset_id_radius);

	//load charge
	hid_t dataset_id_charge = H5Dopen(data->file_id, "/parameters/vmd_structure/charge",H5P_DEFAULT);
	hid_t dataspace_charge = H5Dget_space(dataset_id_charge);
	unsigned long long int dims_charge[1];
	H5Sget_simple_extent_dims(dataspace_charge, dims_charge, NULL );
	//load content of charge
	float data_charge[data->natoms];
	herr_t status_charge= H5Dread(dataset_id_charge, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_charge);
	status_charge = H5Dclose (dataset_id_charge);

	//load atomicnumber
	hid_t dataset_id_atomicnumber = H5Dopen(data->file_id, "/parameters/vmd_structure/atomicnumber",H5P_DEFAULT);
	hid_t dataspace_atomicnumber = H5Dget_space(dataset_id_atomicnumber);
	unsigned long long int dims_atomicnumber[1];
	H5Sget_simple_extent_dims(dataspace_atomicnumber, dims_atomicnumber, NULL );
	//load content of atomicnumber
	int data_atomicnumber[dims_indexOfSpecies[0]];
	herr_t status_atomicnumber= H5Dread(dataset_id_atomicnumber, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_atomicnumber);
	status_atomicnumber = H5Dclose (dataset_id_atomicnumber);


	//load num_segid
	hid_t dataset_id_num_segid = H5Dopen(data->file_id, "/parameters/vmd_structure/num_segid",H5P_DEFAULT);
	hid_t dataspace_num_segid = H5Dget_space(dataset_id_num_segid);
	unsigned long long int dims_num_segid[1];
	H5Sget_simple_extent_dims(dataspace_num_segid, dims_num_segid, NULL );
	//load content of segid
	int data_num_segid[dims_num_segid[0]];
	herr_t status_num_segid= H5Dread(dataset_id_num_segid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_num_segid);
	if(dataspace_num_segid<0 || dims_num_segid[0] != data->natoms){
		printf("ERROR: each atom has to be assigned to a segment\n");
		status_num_segid = H5Dclose (dataset_id_num_segid);
		return MOLFILE_ERROR;
	}
	status_num_segid = H5Dclose (dataset_id_num_segid);


	//load segid
	//open dataset
	hid_t dataset_id_segid = H5Dopen(data->file_id, "/parameters/vmd_structure/segid",H5P_DEFAULT);
	//get datapscae
	hid_t dataspace_segid = H5Dget_space(dataset_id_segid);
	//get dims
	unsigned long long int dims_segid[1];
	H5Sget_simple_extent_dims(dataspace_segid, dims_segid, NULL );

	//load content of segid
	//get datatype and its size
	hid_t filetype_segid= H5Dget_type (dataset_id_segid);
	size_t sdim_segid = H5Tget_size (filetype_segid);
	sdim_segid++;// Make room for null terminator
    	// allocate memory for read buffer.  
	ndims = H5Sget_simple_extent_dims (dataspace_segid, dims_segid, NULL);
	// Allocate array of pointers to rows.
	char **data_segid = (char **) malloc (dims_segid[0] * sizeof (char *));		
	// Allocate space for integer data.
	data_segid[0] = (char *) malloc (dims_segid[0] * sdim_segid * sizeof (char));
	// Set the rest of the pointers to rows to the correct addresses.
	for (int i=1; i<dims_segid[0]; i++)
		data_segid[i] = data_segid[0] + i * sdim_segid;
	// Create the memory datatype.
	hid_t memtype_segid = H5Tcopy (H5T_C_S1);
	herr_t status_segid = H5Tset_size (memtype_segid, sdim_segid);
	// Read the data.
	status_segid = H5Dread (dataset_id_segid, memtype_segid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_segid[0]);
	//close
	status_segid = H5Dclose (dataset_id_segid);
	status_segid = H5Sclose (dataspace_segid);
	status_segid = H5Tclose (filetype_segid);
	status_segid = H5Tclose (memtype_segid);


	//load resid
	hid_t dataset_id_resid = H5Dopen(data->file_id, "/parameters/vmd_structure/resid",H5P_DEFAULT);
	hid_t dataspace_resid = H5Dget_space(dataset_id_resid);
	unsigned long long int dims_resid[1];
	H5Sget_simple_extent_dims(dataspace_resid, dims_resid, NULL );
	//load content of resid
	int data_resid[dims_resid[0]];
	herr_t status_resid= H5Dread(dataset_id_resid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_resid);
	status_resid = H5Dclose (dataset_id_resid);

	//load resname
	//open dataset
	hid_t dataset_id_resname = H5Dopen(data->file_id, "/parameters/vmd_structure/resname",H5P_DEFAULT);
	//get datapscae
	hid_t dataspace_resname = H5Dget_space(dataset_id_resname);
	//get dims
	unsigned long long int dims_resname[1];
	H5Sget_simple_extent_dims(dataspace_resname, dims_resname, NULL );

	//load content of resname
	//get datatype and its size
	hid_t filetype_resname= H5Dget_type (dataset_id_resname);
	size_t sdim_resname = H5Tget_size (filetype_resname);
	sdim_resname++;// Make room for null terminator
    	// allocate memory for read buffer.  
	ndims = H5Sget_simple_extent_dims (dataspace_resname, dims_resname, NULL);
	// Allocate array of pointers to rows.
	char **data_resname = (char **) malloc (dims_resname[0] * sizeof (char *));		
	// Allocate space for integer data.
	data_resname[0] = (char *) malloc (dims_resname[0] * sdim_resname * sizeof (char));
	// Set the rest of the pointers to rows to the correct addresses.
	for (int i=1; i<dims_resname[0]; i++)
		data_resname[i] = data_resname[0] + i * sdim_resname;
	// Create the memory datatype.
	hid_t memtype_resname = H5Tcopy (H5T_C_S1);
	herr_t status_resname = H5Tset_size (memtype_resname, sdim_resname);
	// Read the data.
	status_resname = H5Dread (dataset_id_resname, memtype_resname, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_resname[0]);
	//close
	status_resname = H5Dclose (dataset_id_resname);
	status_resname = H5Sclose (dataspace_resname);
	status_resname = H5Tclose (filetype_resname);
	status_resname = H5Tclose (memtype_resname);

	//load chain
	//open dataset
	hid_t dataset_id_chain = H5Dopen(data->file_id, "/parameters/vmd_structure/chain",H5P_DEFAULT);
	//get datapscae
	hid_t dataspace_chain = H5Dget_space(dataset_id_chain);
	//get dims
	unsigned long long int dims_chain[1];
	H5Sget_simple_extent_dims(dataspace_chain, dims_chain, NULL );

	//load content of chain
	//get datatype and its size
	hid_t filetype_chain= H5Dget_type (dataset_id_chain);
	size_t sdim_chain = H5Tget_size (filetype_chain);
	sdim_chain++;// Make room for null terminator
    	// allocate memory for read buffer.  
	ndims = H5Sget_simple_extent_dims (dataspace_chain, dims_chain, NULL);
	// Allocate array of pointers to rows.
	char **data_chain = (char **) malloc (dims_chain[0] * sizeof (char *));		
	// Allocate space for integer data.
	data_chain[0] = (char *) malloc (dims_chain[0] * sdim_chain * sizeof (char));
	// Set the rest of the pointers to rows to the correct addresses.
	for (int i=1; i<dims_chain[0]; i++)
		data_chain[i] = data_chain[0] + i * sdim_chain;
	// Create the memory datatype.
	hid_t memtype_chain = H5Tcopy (H5T_C_S1);
	herr_t status_chain = H5Tset_size (memtype_chain, sdim_chain);
	// Read the data.
	status_chain = H5Dread (dataset_id_chain, memtype_chain, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_chain[0]);
	//close
	status_chain = H5Dclose (dataset_id_chain);
	status_chain = H5Sclose (dataspace_chain);
	status_chain = H5Tclose (filetype_chain);
	status_chain = H5Tclose (memtype_chain);




	for (int i = 0; i < data->natoms; i++) {
		atom = atoms + i;

		//set species related properties
		//use species information if existing
		if(dataset_id_species>0){
			int index_of_species=find_index_of_species(data_indexOfSpecies,data_species[i],dims_indexOfSpecies[0]);
			
			//set elementname for atom of species
			char* elementname=data_name[index_of_species];
			strncpy(atom->name,elementname,sizeof(elementname));
			//set type for atom of species
			char* type=data_type[index_of_species];
			strncpy(atom->type, type, sizeof(type));
			//set atomic number for atom of species
			if(dataset_id_atomicnumber>0){
				if(data_atomicnumber[index_of_species]<112){
					atom->atomicnumber = data_atomicnumber[index_of_species]%112;
				}else{
					printf("ERROR: Too big atomic number was set in dataset.");
					return MOLFILE_ERROR;
				}
			}else{
				atom->atomicnumber=default_atomicnumber;			
			}
			//set mass for atom of species
			if(dataset_id_mass>0){
				atom->mass = data_mass[index_of_species];
			}else{
				atom->mass = default_mass;		
			}
			//set radius for atom of species
			if(dataset_id_radius>0){
				atom->radius = data_radius[index_of_species];
			}else{
				atom->radius = default_radius;					
			}
		}else{
			//set default values if no species information is provided in /parameters/vmd_structure
			//set default color red (oxygen-> 8)
			strncpy(atom->name,element_symbols[8],sizeof(*element_symbols[8]));
			strncpy(atom->type, atom->name, sizeof(atom->name));
			unsigned int idx = (unsigned int) data_species[i];
			idx=idx%112;
			atom->atomicnumber = idx;
			atom->mass=default_mass;
			atom->radius=default_radius;
		}
		
		//set charge
		if(dataset_id_charge>0){
			atom->charge=data_charge[i];
		}else{
			atom->charge=default_charge;
		}

		//set segid (segments may also consist of pure solvent or lipid)
		if(dataset_id_num_segid>0&&dataset_id_segid>0){
			strncpy(atom->segid,data_segid[data_num_segid[i]],sizeof(data_segid[data_num_segid[i]]));
		}else{
			strncpy(atom->segid,default_segid,sizeof(default_segid));
		}

		if(dataset_id_resid>0){
			//set resid
			atom->resid=data_resid[i];
			//set resname
			if(dataset_id_resname>0){
				strncpy(atom->resname,data_resname[data_resid[i]],sizeof(data_resname[data_resid[i]]));
			}else{
				strncpy(atom->resname,default_resname,sizeof(default_resname));			
			}
			//set chain
			if(dataset_id_chain>0){
				strncpy(atom->chain,data_chain[data_resid[i]],sizeof(data_chain[data_resid[i]]));;
			}else{
				strncpy(atom->chain,default_chain,sizeof(default_chain));
			}


		}else{
			atom->resid = default_resid;
			strncpy(atom->resname,default_resname,sizeof(default_resname));
			strncpy(atom->chain,default_chain,sizeof(default_chain));
		}


	}

	free(*data_name);
	free(data_name);
	free(*data_type);
	free(data_type);
	free(*data_segid);
	free(data_segid);
	free(*data_resname);
	free(data_resname);
	free(*data_chain);
	free(data_chain);

	return MOLFILE_SUCCESS;
}

int read_h5md_structure_no_vmd_structure(void *mydata, int *optflags,molfile_atom_t *atoms) {
	molfile_atom_t *atom;
	h5mddata *data = (h5mddata *) mydata;
	*optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS;
	double data_out[data->natoms];
	hid_t dataset_id_species = H5Dopen(data->file_id, "/particles/atoms/species/value",H5P_DEFAULT);
	herr_t status_species= H5Dread(dataset_id_species, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out);
	H5Tclose(dataset_id_species);
	for (int i = 0; i < data->natoms; i++) {
		atom = atoms + i;
		unsigned int idx = (unsigned int) data_out[i];
		idx=idx%112;
		if(dataset_id_species>0){
			//use species information if existing
			strncpy(atom->name,element_symbols[idx],sizeof(*element_symbols[idx]));
		}else{
			//set default color red (oxygen, 8)
			strncpy(atom->name,element_symbols[8],sizeof(*element_symbols[8]));
		}
		atom->atomicnumber = idx;
		atom->mass = default_mass;
		atom->radius = default_radius;
		strncpy(atom->type, atom->name, sizeof(atom->name));

		atom->resname[0] = default_resname;
		atom->resid = default_resid;
		atom->chain[0] = default_chain;
		atom->segid[0] = default_segid;
	}

	return MOLFILE_SUCCESS;
}


static void get_xyz(void *mydata, int atom_nr, int time_i, double xyz_array[3], double*** data_xyz) {
	//IN: *mydata, atom_nr, time_i
	//OUT: xyz_array[3]
	h5mddata *data = (h5mddata *) mydata;

	xyz_array[0] = data_xyz[time_i][atom_nr][0];
	xyz_array[1] = data_xyz[time_i][atom_nr][1];
	xyz_array[2] = data_xyz[time_i][atom_nr][2];

}


void read_position(h5mddata *data){
	//read position data to data_xyz_read
	double data_xyz_read[data->ntime][data->natoms][data->nspacedims];
	H5Dread(data->dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_xyz_read[0]);
	//allocate memory for data_xyz
	double *** data_xyz = (double ***)malloc(data->ntime*sizeof(double**));

        for (int i = 0; i< data->ntime; i++) {
        	data_xyz[i] = (double **) malloc(data->natoms*sizeof(double *));
        	for (int j = 0; j < data->natoms; j++) {
         		data_xyz[i][j] = (double *)malloc(3*sizeof(double));
        	}
        }
	//copy data of data_xyz_read to data_xyz on heap
        for (int i = 0; i< data->ntime; i++) {
        	for (int j = 0; j < data->natoms; j++) {
         		data_xyz[i][j][0]=data_xyz_read[i][j][0];
			data_xyz[i][j][1]=data_xyz_read[i][j][1];
			data_xyz[i][j][2]=data_xyz_read[i][j][2];
        	}
        }

	data->data_xyz=data_xyz;
}

static int read_h5md_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
	float x, y, z;

	h5mddata *data = (h5mddata *) mydata;
	
	if(position_already_read<0){
		read_position(data);
		position_already_read=1;
	}
	double ***data_xyz=data->data_xyz;
	
	/* read the coordinates */
	unsigned int ntime = data->ntime;
	if (current_time >= ntime - 1) {
		return MOLFILE_ERROR;
	}
	for (int i = 0; i < natoms; i++) {

		current_time = reads / natoms;
		reads += 1;
		double xyz_array[3];
		get_xyz(data, i, current_time, xyz_array,data_xyz);
		x = xyz_array[0];
		y = xyz_array[1];
		z = xyz_array[2];

		if (ts != NULL ) {
			ts->coords[3 * i] = x;
			ts->coords[3 * i + 1] = y;
			ts->coords[3 * i + 2] = z;
		} else {
			break;
		}
	}

	return MOLFILE_SUCCESS;
}

static int h5md_read_bonds_from_file(h5mddata *data){
	printf("Try to load bonds\n\n");
	//load bond_from
	hid_t dataset_id_bond_from = H5Dopen(data->file_id, "/parameters/vmd_structure/bond_from",H5P_DEFAULT);
	hid_t dataspace_bond_from = H5Dget_space(dataset_id_bond_from);
	int dims_bond_from[1];
	int ndims=H5Sget_simple_extent_dims(dataspace_bond_from, dims_bond_from, NULL );
	//load content of bond_from
	int *data_bond_from=(int *)malloc(dims_bond_from[0]*sizeof(int));
	herr_t status_bond_from= H5Dread(dataset_id_bond_from, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_bond_from);
	status_bond_from = H5Dclose (dataset_id_bond_from);
	data->bond_from=data_bond_from;

	//load bond_to
	hid_t dataset_id_bond_to = H5Dopen(data->file_id, "/parameters/vmd_structure/bond_to",H5P_DEFAULT);
	hid_t dataspace_bond_to = H5Dget_space(dataset_id_bond_to);
	int dims_bond_to[1];
	H5Sget_simple_extent_dims(dataspace_bond_to, dims_bond_to, NULL );
	//load content of bond_to
	int *data_bond_to=(int*) malloc(dims_bond_to[0]*sizeof(int));
	herr_t status_bond_to= H5Dread(dataset_id_bond_to, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_bond_to);
	status_bond_to = H5Dclose (dataset_id_bond_to);
	data->bond_to=data_bond_to;

	//save number of bonds
	data->nbonds=dims_bond_to[0];
	printf("nbonds %d\n",dims_bond_to[0]);

	if(dataset_id_bond_from <0 || dataset_id_bond_to <0){
		return -1;
	}else{
		return 0;
	}

}

static int h5md_get_bonds(h5mddata *data, int *nbonds, int **from, int **to, float **bondorder, int **bondtype,  int *nbondtypes, char ***bondtypename){
	int err_read_file = h5md_read_bonds_from_file(data);
	*nbonds=data->nbonds;
	*from=data->bond_from;
	*to=data->bond_to;
	*bondtype=NULL;
	*nbondtypes=0;
	bondtypename=NULL;
	if(err_read_file ==0){
  		return MOLFILE_SUCCESS;
	}else{
		return MOLFILE_ERROR;		
	}
}


static void close_h5md_read(void *mydata) {
	h5mddata *data = (h5mddata *)mydata;
	H5Fclose(data->file_id);
	free(data->bond_from);
	free(data->bond_to);

	//free data_xyz
	for(int i=0;i<data->ntime;i++){
		for(int j=0;j<data->natoms;j++){
			free(data->data_xyz[i][j]);
		}
	free(data->data_xyz[i]);
	}
	free(data->data_xyz);

	free(data);
	
	//reset current_time, reads, position_already_read to start values, so that new molecules are loaded correctly after the first molecule was loaded
	current_time=0;
	reads=0;
	position_already_read =-1;
}

/* registration stuff */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
	memset(&plugin, 0, sizeof(molfile_plugin_t));
	plugin.abiversion = vmdplugin_ABIVERSION;
	plugin.type = MOLFILE_PLUGIN_TYPE;
	plugin.name = "h5md";
	plugin.prettyname = "h5md";
	plugin.author = "Sascha Ehrhardt, Jonas Landsgesell";
	plugin.majorv = 1;
	plugin.minorv = 0;
	plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
	plugin.filename_extension = "h5";
	plugin.open_file_read = open_h5md_read;
	plugin.read_structure = read_h5md_structure;
	plugin.read_next_timestep = read_h5md_timestep;
	plugin.read_bonds = h5md_get_bonds;
	plugin.close_file_read = close_h5md_read;
	return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
	(*cb)(v, (vmdplugin_t *) &plugin);
	return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
	return VMDPLUGIN_SUCCESS;
}
