#ifndef PS_DATA_IO_H
#define PS_DATA_IO_H

/**
 * Read and write a binary file
 */
int ps_data_write_bin(PS_DATA data, FILE *fptr);
int ps_data_read_bin(PS_DATA data, FILE *fptr);

/**
 * Read and write an ASCII file
 */
int ps_data_write_ascii(PS_DATA data, FILE *fptr);
int ps_data_read_ascii(PS_DATA data, FILE *fptr);

/** 
 * Read and write a FITS file
 */
int ps_data_write_fits(PS_DATA data, FILE *fptr);
int ps_data_read_fits(PS_DATA data, FILE *fptr);	

#endif