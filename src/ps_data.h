#ifndef PS_DATA_H
#define PS_DATA_H

/** 
 * Definition of 2-dimensional numerical data structure.
 **/
typedef struct {
	unsigned int xsize;
	unsigned int ysize;
	double xstep;
	double ystep;
	double **data;
} PS_DATA_T;

typedef PS_DATA_T* PS_DATA;

/** Create and destroy functions **/
PS_DATA ps_create_data(unsigned int xsize, unsigned int ysize, double xstep, double ystep);
void ps_destroy_data(PS_DATA data);
	
/** Query the grid size **/
int ps_data_rows(PS_DATA data);
int ps_data_columns(PS_DATA data);

/** 
 * Initialize the data from a two-dimensional array of the correct size. 
 * Returns 0 upon success, error otherwise.
*/
int ps_data_init_with_array(PS_DATA data, double **values);

/**
 * Set the value at a particular row/column
 * Returns 0 upon success, error otherwise.
 */
int ps_data_set_value_at_row_column(PS_DATA data, double value, unsigned int row, unsigned int col);

/**
 * Set the value at a particular x/y
 * Returns 0 upon success, error otherwise.
 */
int ps_data_set_value_at_x_y(PS_DATA data, double value, double x, double y);

/**
 * Get the value at a particular row/column
 * row - the row of interest
 * col - the column of interest
 * val - Pointer to a double where data will be stored
 * Returns 0 upon success, error otherwise
 */
int ps_data_value_at_row_column(PS_DATA data, unsigned int row, unsigned int col, double *val);

/**
 * Get the value at a particular x/y
 * row - the row of interest
 * col - the column of interest
 * val - Pointer to a double where data will be stored
 * Returns 0 upon success, error otherwise
 */
int ps_data_value_at_x_y(PS_DATA data, double x, double y, double *val);

#endif