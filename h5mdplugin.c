/***************************************************************************
 *
 *            (C) Copyright 2013, 2014
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
#include "libh5md.h"

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

//global non constant variables
int reads = 0;
int time_i = 0;

static void *open_h5md_read(const char *filename, const char *filetype, int *natoms){
	h5mddata_lib *data= (h5mddata_lib *) malloc(sizeof(h5mddata_lib));
	read_position_of_file_lib(filename,data);
	*natoms = data->natoms;
	return data;
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
			//add plus 1 if the first occurance of a species happens at the current index i (then a new species is found)
			count+=1;
		}	
	}
	return count;
}


int read_h5md_structure_vmd_structure(h5mddata_lib *data, int *optflags,molfile_atom_t *atoms) {
	//mute errors
 	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	molfile_atom_t *atom;
	*optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS;
	
	h5mdspecies *species=(h5mdspecies*) malloc(sizeof(h5mdspecies));
	int *data_species=(int *) malloc(data->natoms*sizeof(int));
	read_species(data, species, data_species);

	//load indexOfSpecies
	hid_t dataset_id_indexOfSpecies = H5Dopen(data->file_id, "/parameters/vmd_structure/indexOfSpecies",H5P_DEFAULT);
	hid_t dataspace_indexOfSpecies = H5Dget_space(dataset_id_indexOfSpecies);
	unsigned long long int dims_indexOfSpecies[1];
	H5Sget_simple_extent_dims(dataspace_indexOfSpecies, dims_indexOfSpecies, NULL );

	//declaration of variables####################
	int * data_indexOfSpecies;
	hid_t dataset_id_name;
	char **data_name;
	hid_t dataset_id_type;
	char **data_type;
	hid_t dataset_id_mass;
	float *data_mass;
	hid_t dataset_id_radius;
	float *data_radius;
	hid_t dataset_id_atomicnumber;
	int *data_atomicnumber;
	int ndims;

	//check consistency of indexOfSpecies and species/value
	if(count_different_species(data_species,data->natoms)!=dims_indexOfSpecies[0]){
		printf("ERROR: /parameters/vmd_structure/indexOfSpecies does not contain as much different species as different species are present in /particles/atoms/species/value !\n");
		printf("Skipping species related data.\n");
		//return MOLFILE_ERROR;
	}else{
		//load content of indexOfSpecies
		data_indexOfSpecies=(int*) malloc(dims_indexOfSpecies[0]*sizeof(int));
		H5Dread(dataset_id_indexOfSpecies, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_indexOfSpecies);
		H5Dclose (dataset_id_indexOfSpecies);
	

		//load name
		//Beginn der HDF5 fehlermeldungen, aber ohne absturz, kein Problem fuer lauffaehigkeit (eventuell ausschalten)
		//open dataset
		dataset_id_name = H5Dopen(data->file_id, "/parameters/vmd_structure/name",H5P_DEFAULT);
		//get dataspace
		hid_t dataspace_name = H5Dget_space(dataset_id_name);
		//get dims
		unsigned long long int dims_name[1];
		H5Sget_simple_extent_dims(dataspace_name, dims_name, NULL );

		if(dataset_id_name>0){
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
			ndims = H5Sget_simple_extent_dims(dataspace_name, dims_name, NULL);
			// Allocate array of pointers to rows.
			data_name = (char **) malloc (dims_name[0] * sizeof (char *));		
			// Allocate space for integer data.
			data_name[0] = (char *) malloc (dims_name[0] * sdim * sizeof (char));
			// Set the rest of the pointers to rows to the correct addresses.
			for (int i=1; i<dims_name[0]; i++)
				data_name[i] = data_name[0] + i * sdim;
			// Create the memory datatype.
			hid_t memtype = H5Tcopy (H5T_C_S1);
			H5Tset_size (memtype, sdim);
			// Read the data.
			H5Dread (dataset_id_name, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_name[0]);
			//close
			H5Dclose (dataset_id_name);
			H5Sclose (dataspace_name);
			H5Tclose (filetype);
			H5Tclose (memtype);
		}
		//load type
		//open dataset
		dataset_id_type = H5Dopen(data->file_id, "/parameters/vmd_structure/type",H5P_DEFAULT);
		//get datapscae
		hid_t dataspace_type = H5Dget_space(dataset_id_type);
		//get dims
		unsigned long long int dims_type[1];
		H5Sget_simple_extent_dims(dataspace_type, dims_type, NULL );

		if(dataset_id_type>0){
			//load content of type
			//get datatype and its size
			hid_t filetype_type= H5Dget_type (dataset_id_type);
			size_t sdim_type = H5Tget_size (filetype_type);
			sdim_type++;// Make room for null terminator
		    	// allocate memory for read buffer.  
			ndims = H5Sget_simple_extent_dims (dataspace_type, dims_type, NULL);
			// Allocate array of pointers to rows.
			data_type = (char **) malloc (dims_type[0] * sizeof (char *));		
			// Allocate space for integer data.
			data_type[0] = (char *) malloc (dims_type[0] * sdim_type * sizeof (char));
			// Set the rest of the pointers to rows to the correct addresses.
			for (int i=1; i<dims_type[0]; i++)
				data_type[i] = data_type[0] + i * sdim_type;
			// Create the memory datatype.
			hid_t memtype_type = H5Tcopy (H5T_C_S1);
			H5Tset_size (memtype_type, sdim_type);
			// Read the data.
			H5Dread (dataset_id_type, memtype_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_type[0]);
			//close
			H5Dclose (dataset_id_type);
			H5Sclose (dataspace_type);
			H5Tclose (filetype_type);
			H5Tclose (memtype_type);
		}

		//load mass
		dataset_id_mass = H5Dopen(data->file_id, "/parameters/vmd_structure/mass",H5P_DEFAULT);
		hid_t dataspace_mass = H5Dget_space(dataset_id_mass);
		unsigned long long int dims_mass[1];
		H5Sget_simple_extent_dims(dataspace_mass, dims_mass, NULL );
		//load content of mass
		data_mass=(float*) malloc(dims_indexOfSpecies[0]*sizeof(float));
		H5Dread(dataset_id_mass, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_mass);
		H5Dclose (dataset_id_mass);
		//load radius
		dataset_id_radius = H5Dopen(data->file_id, "/parameters/vmd_structure/radius",H5P_DEFAULT);
		hid_t dataspace_radius = H5Dget_space(dataset_id_radius);
		unsigned long long int dims_radius[1];
		H5Sget_simple_extent_dims(dataspace_radius, dims_radius, NULL );
		//load content of radius
		data_radius=(float*) malloc(dims_indexOfSpecies[0]*sizeof(float));
		H5Dread(dataset_id_radius, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_radius);
		H5Dclose (dataset_id_radius);
		//load atomicnumber
		dataset_id_atomicnumber = H5Dopen(data->file_id, "/parameters/vmd_structure/atomicnumber",H5P_DEFAULT);
		hid_t dataspace_atomicnumber = H5Dget_space(dataset_id_atomicnumber);
		unsigned long long int dims_atomicnumber[1];
		H5Sget_simple_extent_dims(dataspace_atomicnumber, dims_atomicnumber, NULL );
		//load content of atomicnumber
		data_atomicnumber=(int*) malloc(dims_indexOfSpecies[0]*sizeof(int));
		H5Dread(dataset_id_atomicnumber, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_atomicnumber);
		H5Dclose (dataset_id_atomicnumber);
	}

	//load charge
	hid_t dataset_id_charge = H5Dopen(data->file_id, "/parameters/vmd_structure/charge",H5P_DEFAULT);
	hid_t dataspace_charge = H5Dget_space(dataset_id_charge);
	unsigned long long int dims_charge[1];
	H5Sget_simple_extent_dims(dataspace_charge, dims_charge, NULL );
	//load content of charge
	float data_charge[data->natoms];
	H5Dread(dataset_id_charge, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_charge);
	H5Dclose (dataset_id_charge);


	//Segid related data handling
	//declaration of variables
	hid_t dataset_id_segid;
	int *data_indexOfSegid;
	char **data_segid;

	//load indexOfSegid
	hid_t dataset_id_indexOfSegid = H5Dopen(data->file_id, "/parameters/vmd_structure/indexOfSegid",H5P_DEFAULT);
	hid_t dataspace_indexOfSegid = H5Dget_space(dataset_id_indexOfSegid);
	unsigned long long int dims_indexOfSegid[1];
	H5Sget_simple_extent_dims(dataspace_indexOfSegid, dims_indexOfSegid, NULL );
	if(dataset_id_indexOfSegid>0){
		//load content of segid
		data_indexOfSegid=(int*) malloc(dims_indexOfSegid[0]*sizeof(int));
		H5Dread(dataset_id_indexOfSegid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_indexOfSegid);
		H5Dclose (dataset_id_indexOfSegid);
		//load segid
		//open dataset
		dataset_id_segid = H5Dopen(data->file_id, "/parameters/vmd_structure/segid",H5P_DEFAULT);
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
		data_segid = (char **) malloc (dims_segid[0] * sizeof (char *));		
		// Allocate space for integer data.
		data_segid[0] = (char *) malloc (dims_segid[0] * sdim_segid * sizeof (char));
		// Set the rest of the pointers to rows to the correct addresses.
		for (int i=1; i<dims_segid[0]; i++)
			data_segid[i] = data_segid[0] + i * sdim_segid;
		// Create the memory datatype.
		hid_t memtype_segid = H5Tcopy (H5T_C_S1);
		H5Tset_size (memtype_segid, sdim_segid);
		// Read the data.
		H5Dread (dataset_id_segid, memtype_segid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_segid[0]);
		//close
		H5Dclose (dataset_id_segid);
		H5Sclose (dataspace_segid);
		H5Tclose (filetype_segid);
		H5Tclose (memtype_segid);
	}


	//load resid
	hid_t dataset_id_resid = H5Dopen(data->file_id, "/parameters/vmd_structure/resid",H5P_DEFAULT);
	hid_t dataspace_resid = H5Dget_space(dataset_id_resid);
	unsigned long long int dims_resid[1];
	H5Sget_simple_extent_dims(dataspace_resid, dims_resid, NULL );

	//resid related data handling
	//declaration of variables
	int* data_resid;
	hid_t dataset_id_resname;
	char **data_resname;
	hid_t dataset_id_chain;
	char **data_chain;

	if(dataset_id_resid>0){
		//load content of resid
		data_resid=(int*) malloc(dims_resid[0]*sizeof(int));
		H5Dread(dataset_id_resid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_resid);
		H5Dclose (dataset_id_resid);
		//load resname
		//open dataset
		dataset_id_resname = H5Dopen(data->file_id, "/parameters/vmd_structure/resname",H5P_DEFAULT);
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
		data_resname = (char **) malloc (dims_resname[0] * sizeof (char *));		
		// Allocate space for integer data.
		data_resname[0] = (char *) malloc (dims_resname[0] * sdim_resname * sizeof (char));
		// Set the rest of the pointers to rows to the correct addresses.
		for (int i=1; i<dims_resname[0]; i++)
			data_resname[i] = data_resname[0] + i * sdim_resname;
		// Create the memory datatype.
		hid_t memtype_resname = H5Tcopy (H5T_C_S1);
		H5Tset_size (memtype_resname, sdim_resname);
		// Read the data.
		H5Dread (dataset_id_resname, memtype_resname, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_resname[0]);
		//close
		H5Dclose (dataset_id_resname);
		H5Sclose (dataspace_resname);
		H5Tclose (filetype_resname);
		H5Tclose (memtype_resname);
		//load chain
		//open dataset
		dataset_id_chain = H5Dopen(data->file_id, "/parameters/vmd_structure/chain",H5P_DEFAULT);
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
		data_chain = (char **) malloc (dims_chain[0] * sizeof (char *));		
		// Allocate space for integer data.
		data_chain[0] = (char *) malloc (dims_chain[0] * sdim_chain * sizeof (char));
		// Set the rest of the pointers to rows to the correct addresses.
		for (int i=1; i<dims_chain[0]; i++)
			data_chain[i] = data_chain[0] + i * sdim_chain;
		// Create the memory datatype.
		hid_t memtype_chain = H5Tcopy (H5T_C_S1);
		H5Tset_size (memtype_chain, sdim_chain);
		// Read the data.
		H5Dread (dataset_id_chain, memtype_chain, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_chain[0]);
		//close
		H5Dclose (dataset_id_chain);
		H5Sclose (dataspace_chain);
		H5Tclose (filetype_chain);
		H5Tclose (memtype_chain);
	}

	//give data to VMD
	for (int i = 0; i < data->natoms; i++) {
		atom = atoms + i;
		//set species related properties
		//use species information if existing
		if(species->dataset_id_species>0 && dataset_id_indexOfSpecies>0 ){
			int index_of_species=find_index_of_species(data_indexOfSpecies,data_species[i],dims_indexOfSpecies[0]);
			if(dataset_id_name>0){
				//set elementname for atom of species
				char* elementname=data_name[index_of_species];
				strncpy(atom->name,elementname,sizeof(elementname));
			}else{
				strncpy(atom->name,default_name,sizeof(default_name));
			}
			//set type for atom of species
			if(dataset_id_type>0){
				char* type=data_type[index_of_species];
				strncpy(atom->type, type, sizeof(type));
			}else{	
				strncpy(atom->type, default_type, sizeof(default_type));			
			}
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
		if(dataset_id_indexOfSegid>0&&dataset_id_segid>0 && (dims_indexOfSegid[0] == data->natoms)){
			strncpy(atom->segid,data_segid[data_indexOfSegid[i]],sizeof(data_segid[data_indexOfSegid[i]]));
		}else{
			//printf("ERROR: segid not present or not every atom is assigned to a segment\n");
			//printf("using default values \n");
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

	//free
	if(dataset_id_indexOfSpecies>0){
		free(data_indexOfSpecies);
		if(dataset_id_name>0){
			free(*data_name);
			free(data_name);
		}
		if(dataset_id_type>0){
			free(*data_type);
			free(data_type);
		}
		if(dataset_id_mass>0){
			free(data_mass);
		}
		if(dataset_id_radius>0){
			free(data_radius);
		}
		if(dataset_id_atomicnumber>0){
			free(data_atomicnumber);		
		}
	}
	if(dataset_id_indexOfSegid>0){
		free(data_indexOfSegid);
		free(*data_segid);
		free(data_segid);
	}
	if(dataset_id_resid>0){
		free(data_resid);
		free(*data_resname);
		free(data_resname);
		free(*data_chain);
		free(data_chain);
	}

	free(species->data_species);
	free(species);
	//show errors
	H5Eset_auto(H5E_DEFAULT, H5Eprint, NULL);
	return MOLFILE_SUCCESS;
}

//Fixed interface due to VMD (don t outsource to lib)
static int read_h5md_timestep(h5mddata_lib *data, int natoms, molfile_timestep_t *ts) {
	/* read the coordinates */
	unsigned int ntime = data->ntime;
	if (time_i >= ntime - 1) {
		return MOLFILE_ERROR;
	}
	for (int i = 0; i < natoms; i++) {

		time_i = reads / natoms;
		reads += 1;
		if (ts != NULL ) {
			ts->coords[3 * i] = data->data_xyz[time_i][i][0];
			ts->coords[3 * i + 1] = data->data_xyz[time_i][i][1];
			ts->coords[3 * i + 2] = data->data_xyz[time_i][i][2];
		} else {
			break;
		}
	}
	return MOLFILE_SUCCESS;
}

static int h5md_read_bonds_from_file(h5mddata_lib *data){
	//mute errors
 	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	//Try to load bonds
	//load bond_from
	hid_t dataset_id_bond_from = H5Dopen(data->file_id, "/parameters/vmd_structure/bond_from",H5P_DEFAULT);
	int dims_bond_from[1];
	//load content of bond_from
	int *data_bond_from=(int *)malloc(dims_bond_from[0]*sizeof(int));
	H5Dread(dataset_id_bond_from, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_bond_from);
	H5Dclose (dataset_id_bond_from);
	data->bond_from=data_bond_from;

	//load bond_to
	hid_t dataset_id_bond_to = H5Dopen(data->file_id, "/parameters/vmd_structure/bond_to",H5P_DEFAULT);
	hid_t dataspace_bond_to = H5Dget_space(dataset_id_bond_to);
	int dims_bond_to[1];
	H5Sget_simple_extent_dims(dataspace_bond_to, dims_bond_to, NULL );
	//load content of bond_to
	int *data_bond_to=(int*) malloc(dims_bond_to[0]*sizeof(int));
	H5Dread(dataset_id_bond_to, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_bond_to);
	H5Dclose (dataset_id_bond_to);
	data->bond_to=data_bond_to;

	//save number of bonds
	data->nbonds=dims_bond_to[0];
	
	//show errors
	H5Eset_auto(H5E_DEFAULT, H5Eprint, NULL);
	
	if(dataset_id_bond_from <0 || dataset_id_bond_to <0){
		return -1;
	}else{
		return 0;
	}
}

static int h5md_get_bonds(h5mddata_lib *data, int *nbonds, int **from, int **to, float **bondorder, int **bondtype,  int *nbondtypes, char ***bondtypename){
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


static void close_h5md_read(h5mddata_lib *data) {
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
	
	//reset time_i, reads, position_already_read to start values, so that new molecules are loaded correctly after the first molecule was loaded
	time_i=0;
	reads=0;
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
	plugin.read_structure = read_h5md_structure_vmd_structure;
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
