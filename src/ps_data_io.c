#include "ps_data_io.h"
#include "ps_errors.h"
#include <stdlib.h>

PS_DATA ps_data_read_ascii(FILE *fptr) {
	unsigned int rows, cols;
	unsigned int i,j;
	double *d;
	PS_DATA retval;
	
	// Seek to beginning of file
	fseek(fptr, 0, SEEK_SET);
	
	// Read the rows and columns
	fscanf(fptr, "%i\n%i\n", &cols, &rows);
	printf("File at %i\n", ftell(fptr));
	
	// Create the data
	retval = ps_data_create(rows, cols);
	
	// Read the X values
	d = (double*)malloc(sizeof(double)*cols);
	double dval;
	for(i=0; i<cols; i++) {
		fscanf(fptr,"%e ", &dval);
		printf("File at %i\n", ftell(fptr));
//		fseek(fptr, 1, SEEK_CUR);
		printf("%g ", dval);
	}
	ps_data_set_x_values(retval, d);
	free(d);
	
	// Read the Y values
	d = (double*)malloc(sizeof(double)*rows);
	for(i=0; i<rows; i++) {
		fscanf(fptr, "%e", (d+i));
	}
	ps_data_set_y_values(retval, d);
	free(d);
	
	return retval;
}

int ps_data_write_ascii(PS_DATA data, FILE *fptr) {
	int rows = ps_data_rows(data);
	int cols = ps_data_columns(data);
	int i,j;
	
	fprintf(fptr, "%i\n", cols);
	fprintf(fptr, "%i\n", rows);
	for(i=0; i<cols;i++) {
		fprintf(fptr, "%e", ps_data_xvalue_at(data, i));
		if (cols-1 == i) {
			fprintf(fptr, "\n");
		} else {
			fprintf(fptr, " ");
		}
	}

	for(j=0; j<rows;j++) {
		fprintf(fptr, "%e", ps_data_yvalue_at(data,j));
		if (rows-1 == j) {
			fprintf(fptr, "\n");
		} else {
			fprintf(fptr, " ");
		}		
	}
	
	for (i=0; i<rows;i++) {
		for(j=0; j<cols;j++) {
			fprintf(fptr, "%e", ps_data_value(data, i, j));
			if (cols-1 == j) {
				fprintf(fptr, "\n");
			} else {
				fprintf(fptr, " ");
			}					
		}
	}
	fprintf(fptr, "\n");	
	return PS_OK;
}

int ps_data_write_bin(PS_DATA data, FILE *fptr) {
	int rows = ps_data_rows(data);
	int cols = ps_data_columns(data);
	int i, j;
	
	// Write the number of cols
	ps_write_float64(fptr, (double)cols);
	
	// Write the number of rows
	ps_write_float64(fptr, (double)rows);
	
	// Write out x-values
	for(j=0; j < cols; j++) {
		ps_write_float64(fptr, ps_data_xvalue_at(data, j));
	}
	
	// Write out y-values
	for(i=0; i < rows; i++) {
		ps_write_float64(fptr, ps_data_yvalue_at(data, i));
	}
	
	// Write out data
	for(i = 0; i < rows; i++) {
		for (j=0; j < cols; j++) {
			ps_write_float64(fptr, ps_data_value(data, i, j));
		}
	}
	
}

PS_DATA ps_data_read_bin(FILE *fptr) {

	unsigned int rows, cols;
	unsigned int i,j;
	double *d;
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
	
	// Create the structure
	retval = ps_data_create(cols, rows); 
	
	// Read the X values
	d = (double*)malloc(sizeof(double)*cols);
	for (i = 0; i < cols && PS_OK == ps_errno; i++) {
		*(d+i) = ps_read_float64(fptr);
	}

	if (PS_OK != ps_errno)  {
		free(d);
		ps_data_destroy(retval);
		return NULL;
	}

	ps_data_set_x_values(retval, d);
	free(d);

	// Read in Y Values
	d = (double*)malloc(sizeof(double)*rows);
	for (j = 0; j < rows && PS_OK == ps_errno; j++) {
		*(d+j) = ps_read_float64(fptr);		
	}
	
	if (PS_OK != ps_errno) {
		free(d);
		ps_data_destroy(retval);
		return NULL;		
	}
	
	ps_data_set_y_values(retval, d);
	free(d);

	d = (double*)malloc(sizeof(double)*rows*cols);
	for (i=0; i < rows; i++) {
		for (j=0; j < cols; j++) {
			d[i*cols+j] = ps_read_float64(fptr);
			if (PS_OK != ps_errno) {
				fprintf(stderr, "Error in ps_data_io#read_bin\n");
				break;
			}
		}
	}

	if (PS_OK != ps_errno) {
		free(d);
		ps_data_destroy(retval);
		return NULL;		
	}
	
	// Data formats don't match up right now. 
	// PS_DATA stores an xsize,ysize.  but file format
	// allows rectilinear with varying xsize,ysize
	ps_data_init_with_array(retval, d, rows, cols);
	free(d);
	
	return retval;
    
}

double ps_read_float64(FILE *fptr) {
	int i;
	FLOAT64 value;
	int nbytes = 1;
	for (i=0; i<8 && 1 == nbytes; i++) {
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

void ps_write_float64(FILE *fptr, double d) {
	FLOAT64 value;
	value.d = d;
	int i;
	for (i=0; i<8; i++) {
		fwrite(&(value.c[i]), 1, 1, fptr);
	}
}