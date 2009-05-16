#ifndef PS_DATA_IO_H
#define PS_DATA_IO_H

#include <stdio.h>
#include "ps_data.h"

typedef union float64 {
	char c[8];
	double d;
} FLOAT64;

/**
 * Read 8bytes from the file and store double in ptr
 */
double ps_read_float64(FILE *fptr);
void ps_write_float64(FILE *fprt, double d);
unsigned int ps_read_int(FILE *fptr);

/**
 * Read and write a binary file
 */
int ps_data_write_bin(PS_DATA data, FILE *fptr);
PS_DATA ps_data_read_bin(FILE *fptr);

/**
 * Read and write an ASCII file
 */
int ps_data_write_ascii(PS_DATA data, FILE *fptr);
PS_DATA ps_data_read_ascii(FILE *fptr);

/** 
 * Read and write a FITS file
 */
int ps_data_write_fits(PS_DATA data, FILE *fptr);
int ps_data_read_fits(PS_DATA data, FILE *fptr);	

#endif