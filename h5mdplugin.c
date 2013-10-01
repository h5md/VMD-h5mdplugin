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

typedef struct {
	int ntime;
	int natoms;
	int nspacedims;
	molfile_atom_t *atomlist;
	char *file_name;
	hid_t file_id;
	hid_t dataset_id;
} h5mddata;


int reads = 0;
int current_time = 0;

static void *open_h5md_read(const char *filename, const char *filetype, int *natoms) {
	h5mddata *data;
	data = (h5mddata *) malloc(sizeof(h5mddata));

	hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	hid_t dataset_id = H5Dopen(file_id, "/particles/atoms/position/value",H5P_DEFAULT);
	hid_t dataspace = H5Dget_space(dataset_id);

	unsigned long long int dims_out[3];
	H5Sget_simple_extent_dims(dataspace, dims_out, NULL );

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
	hid_t dataset_id_indexOfSpecies = H5Dopen(data->file_id, "/parameters/vmd_structure/indexOfSpecies",H5P_DEFAULT);
	if(dataset_id_indexOfSpecies<0){
		error_status_of_read = read_h5md_structure_no_vmd_structure(mydata,optflags,atoms);
	}else{
		error_status_of_read = read_h5md_structure_vmd_structure(mydata,optflags,atoms);
	}
	return error_status_of_read;
}

int find_index_of_species(int *data_out_indexOfSpecies, int species,int len_index){
	//find index of species "species" in data_out_index_Of_Species which has length len_index
	int ndx=-1;
	for(int i=0;i<len_index;i++){
		if(species== data_out_indexOfSpecies[2*i+1]){
			ndx=i;
			break;
		}
	}

	return ndx;
}

int read_h5md_structure_vmd_structure(void *mydata, int *optflags,molfile_atom_t *atoms) {
	molfile_atom_t *atom;
	h5mddata *data = (h5mddata *) mydata;

	//check consistency of parameters/vmd_structure
	hid_t dataset_id_indexOfSpecies = H5Dopen(data->file_id, "/parameters/vmd_structure/indexOfSpecies",H5P_DEFAULT);
	hid_t dataspace_indexOfSpecies = H5Dget_space(dataset_id_indexOfSpecies);
	unsigned long long int dims_out_indexOfSpecies[2];
	H5Sget_simple_extent_dims(dataspace_indexOfSpecies, dims_out_indexOfSpecies, NULL );
	int len_index=dims_out_indexOfSpecies[0];

	hid_t dataset_id_name = H5Dopen(data->file_id, "/parameters/vmd_structure/name",H5P_DEFAULT);
	hid_t dataspace_name = H5Dget_space(dataset_id_name);
	unsigned long long int dims_out_name[1];
	H5Sget_simple_extent_dims(dataspace_name, dims_out_name, NULL );

	hid_t dataset_id_type = H5Dopen(data->file_id, "/parameters/vmd_structure/type",H5P_DEFAULT);
	hid_t dataspace_type = H5Dget_space(dataset_id_type);
	unsigned long long int dims_out_type[1];
	H5Sget_simple_extent_dims(dataspace_type, dims_out_type, NULL );


	if(dims_out_indexOfSpecies[0] == dims_out_name[0] && dims_out_name[0] == dims_out_type[0]){
		*optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS;
		//get content

		//*species
		int data_out_species[data->natoms];
		hid_t dataset_species_id = H5Dopen(data->file_id, "/particles/atoms/species/value",H5P_DEFAULT);
		herr_t status_species= H5Dread(dataset_species_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out_species);
		status_species = H5Dclose (dataset_species_id);
		

		//*indexOfSpecies
		int *data_out_indexOfSpecies;
		data_out_indexOfSpecies=(int *) malloc(dims_out_indexOfSpecies[0]*dims_out_indexOfSpecies[1]*sizeof(int));
		hid_t dataset_indexOfSpecies_id = H5Dopen(data->file_id, "/parameters/vmd_structure/indexOfSpecies",H5P_DEFAULT);
		herr_t status_indexOfSpecies= H5Dread(dataset_indexOfSpecies_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out_indexOfSpecies);
		status_indexOfSpecies = H5Dclose (dataset_indexOfSpecies_id);

		//*name
		hsize_t     dims[1] = {dims_out_indexOfSpecies[0]};
		//open dataset
		hid_t dset_name= H5Dopen (data->file_id, "/parameters/vmd_structure/name", H5P_DEFAULT);
		//get datatype and its size
		hid_t filetype= H5Dget_type (dset_name);
		size_t sdim = H5Tget_size (filetype);
		sdim++;// Make room for null terminator
	        /*
    	        * Get dataspace and allocate memory for read buffer.  This is a
     	        * two dimensional dataset so the dynamic allocation must be done
                * in steps.
                */
		hid_t space = H5Dget_space (dset_name);
		int ndims = H5Sget_simple_extent_dims (space, dims, NULL);
	        // Allocate array of pointers to rows.
	        char **rdata = (char **) malloc (dims[0] * sizeof (char *));		
	        // Allocate space for integer data.
	        rdata[0] = (char *) malloc (dims[0] * sdim * sizeof (char));
	        // Set the rest of the pointers to rows to the correct addresses.
	        for (int i=1; i<dims[0]; i++)
			rdata[i] = rdata[0] + i * sdim;
	        // Create the memory datatype.
	        hid_t memtype = H5Tcopy (H5T_C_S1);
	        herr_t status = H5Tset_size (memtype, sdim);
	        // Read the data.
	        status = H5Dread (dset_name, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata[0]);
		//close
		status = H5Dclose (dset_name);
	        status = H5Sclose (space);
	        status = H5Tclose (filetype);
	        status = H5Tclose (memtype);


		//*type
		hsize_t     dims_type[1] = {dims_out_indexOfSpecies[0]};
		//open dataset
		hid_t dset_type= H5Dopen (data->file_id, "/parameters/vmd_structure/type", H5P_DEFAULT);
		//get datatype and its size
		hid_t filetype_type= H5Dget_type (dset_type);
		size_t sdim_type = H5Tget_size (filetype_type);
		sdim_type++;// Make room for null terminator
    	        // Get dataspace and allocate memory for read buffer.  
		space = H5Dget_space (dset_type);
		ndims = H5Sget_simple_extent_dims (space, dims_type, NULL);
	        // Allocate array of pointers to rows.
	        char **rdata_type = (char **) malloc (dims_type[0] * sizeof (char *));		
	        // Allocate space for integer data.
	        rdata_type[0] = (char *) malloc (dims_type[0] * sdim_type * sizeof (char));
	        // Set the rest of the pointers to rows to the correct addresses.
	        for (int i=1; i<dims_type[0]; i++)
			rdata_type[i] = rdata_type[0] + i * sdim_type;
	        // Create the memory datatype.
	        hid_t memtype_type = H5Tcopy (H5T_C_S1);
	        status = H5Tset_size (memtype_type, sdim_type);
	        // Read the data.
	        status = H5Dread (dset_type, memtype_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_type[0]);
		//close
		status = H5Dclose (dset_type);
	        status = H5Sclose (space);
	        status = H5Tclose (filetype_type);
	        status = H5Tclose (memtype_type);

		for (int i = 0; i < data->natoms; i++) {
			atom = atoms + i;
			unsigned int idx = (unsigned int) data_out_species[i];
			if(status_species==0){
				//use species information if existing
				int ndx=find_index_of_species(data_out_indexOfSpecies,data_out_species[i],len_index);
				char* elementname=rdata[ndx];
				strncpy(atom->name,elementname,sizeof(elementname));
				char* type=rdata_type[ndx];
				strncpy(atom->type, type, sizeof(type));
				idx=idx%112;
			}else{
				//set default color red (oxygen-> 8)
				strncpy(atom->name,element_symbols[8],sizeof(*element_symbols[8]));
				strncpy(atom->type, atom->name, sizeof(atom->name));
			}
			atom->atomicnumber = idx;
			atom->mass = 1.0;
			atom->radius = 0.5;
			atom->resname[0] = '\0';
			atom->resid = 1;
			atom->chain[0] = '\0';
			atom->segid[0] = '\0';
		}

		return MOLFILE_SUCCESS;
	}else{
			return MOLFILE_ERROR;
	}

}


int read_h5md_structure_no_vmd_structure(void *mydata, int *optflags,molfile_atom_t *atoms) {
	molfile_atom_t *atom;
	h5mddata *data = (h5mddata *) mydata;

	*optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS;
	double data_out[data->natoms];
	hid_t dataset_species_id = H5Dopen(data->file_id, "/particles/atoms/species/value",H5P_DEFAULT);
	herr_t status= H5Dread(dataset_species_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out);
	status = H5Tclose(dataset_species_id);
	for (int i = 0; i < data->natoms; i++) {
		atom = atoms + i;
		unsigned int idx = (unsigned int) data_out[i];
		idx=idx%112;
		if(status==0){
			//use species information if existing
			strncpy(atom->name,element_symbols[idx],sizeof(*element_symbols[idx]));
		}else{
			//set default color red (oxygen, 8)
			strncpy(atom->name,element_symbols[8],sizeof(*element_symbols[8]));
		}
		atom->atomicnumber = idx;
		atom->mass = 1.0;
		atom->radius = 0.5;
		strncpy(atom->type, atom->name, sizeof(atom->name));

		atom->resname[0] = '\0';
		atom->resid = 1;
		atom->chain[0] = '\0';
		atom->segid[0] = '\0';
	}

	return MOLFILE_SUCCESS;
}

static void get_xyz(void *mydata, int atom_nr, int time_i,
		double xyz_array[3]) {
	h5mddata *data = (h5mddata *) mydata;
	double data_out[data->ntime][data->natoms][data->nspacedims];
	H5Dread(data->dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out);

	xyz_array[0] = data_out[time_i][atom_nr][0];
	xyz_array[1] = data_out[time_i][atom_nr][1];
	xyz_array[2] = data_out[time_i][atom_nr][2];
}

static int read_h5md_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
	float x, y, z;

	h5mddata *data = (h5mddata *) mydata;

	/* read the coordinates */
	unsigned int ntime = data->ntime;
	if (current_time >= ntime - 1) {
		return MOLFILE_ERROR;
	}
	for (int i = 0; i < natoms; i++) {

		current_time = reads / natoms;
		reads += 1;
		double xyz_array[3];
		get_xyz(data, i, current_time, xyz_array);
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

static void close_h5md_read(void *mydata) {
	h5mddata *data = (h5mddata *)mydata;
	H5Fclose(data->file_id);
	free(data);
	//reset current_time and reads to zero, so that new molecules are loaded correctly after the first molecule was loaded
	current_time=0;
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
	plugin.read_structure = read_h5md_structure;
	plugin.read_next_timestep = read_h5md_timestep;
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
