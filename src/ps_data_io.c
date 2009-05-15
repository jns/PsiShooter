#include "ps_data_io.h"
#include "ps_errors.h"
#include <stdlib.h>

int ps_data_write_bin(PS_DATA data, FILE *fptr) {
	
}

PS_DATA ps_data_read_bin(FILE *fptr) {

	unsigned int rows, cols;
	unsigned int i,j;
	double *d;
	double read_val;
	PS_DATA retval;
	
	// Seek to beginning of file
	fseek(fptr, 0, SEEK_SET);
	
	// Read 8byte float for number of columns
	cols = ps_read_int(fptr);
	
	if (PS_OK != ps_errno) 
		return NULL;
	
	rows = ps_read_int(fptr);

	if (PS_OK != ps_errno) 
		return NULL;
	
	// Read the x and y values
	for (i = 0; i < rows && PS_OK == ps_errno; i++) {
		read_val = ps_read_float64(fptr);
	}

	if (PS_OK != ps_errno) 
		return NULL;

	for (j = 0; j < cols && PS_OK == ps_errno; j++) {
		read_val = ps_read_float64(fptr);		
	}
	
	if (PS_OK != ps_errno) 
		return NULL;

	d = (double*)malloc(sizeof(double)*rows*cols);
	for (i=0; i < rows; i++) {
		for (j=0; j < cols; j++) {
			*(d+i*cols+j) = ps_read_float64(fptr);
			if (PS_OK != ps_errno) 
				break;
		}
	}

	if (PS_OK != ps_errno) 
		return NULL;
	
	// Data formats don't match up right now. 
	// PS_DATA stores an xsize,ysize.  but file format
	// allows rectilinear with varying xsize,ysize
	retval = ps_create_data(cols, rows, 1, 1); 
	ps_data_init_with_array(retval, d, rows, cols);
	free(d);
	
	return retval;
    
}

double ps_read_float64(FILE *fptr) {
	int i;
	FLOAT64 value;
	int nbytes = 1;
	for (i=7; i>=0 && 1 == nbytes; i--) {
		nbytes = fread(&(value.c[i]), 1, 1, fptr);
	}
	
	if (1 == nbytes)
		ps_errno = PS_OK;
	else	
		ps_errno = PS_ERROR_UNEXPECTED_EOF;

	return value.d;
}

unsigned int ps_read_int(FILE *fptr) {
	return (unsigned int)ps_read_float64(fptr);
}