#ifndef PS_DATA_H
#define PS_DATA_H

#include "ps_errors.h"

/** 
 * Definition of 2-dimensional numerical data structure.
 **/
typedef struct {
	unsigned int xsize;
	unsigned int ysize;
	double xstep; // temporary for the sake of not upsetting psi-shooter
	double ystep; // temporary for the sake of not upsetting psi-shooter
	double *x_values;
	double *y_values;
	double *data;
} PS_DATA_T;

typedef PS_DATA_T* PS_DATA;

#define ps_data_value(ps_data,r,c) (ps_data->data[ps_data->xsize*r + c])

/** Create and destroy functions **/
PS_DATA ps_create_data(unsigned int xsize, unsigned int ysize, double xstep, double ystep); // DEPRECATED use ps_data_create
PS_DATA ps_data_create(unsigned int rows, unsigned int cols);
PS_DATA ps_data_copy(PS_DATA data);

void ps_destroy_data(PS_DATA data); // DEPRECATED use ps_data_destroy
void ps_data_destroy(PS_DATA data);
	
/** Query the grid size **/
int ps_data_rows(PS_DATA data);
int ps_data_columns(PS_DATA data);

/** Set the x and y values */
int ps_data_set_x_values(PS_DATA data, double *x_values);
int ps_data_set_y_values(PS_DATA data, double *y_values);
int ps_data_set_x_value_at(PS_DATA data, int x, double x_value);
int ps_data_set_y_value_at(PS_DATA data, int y, double y_value);

/**
 * Query the x and y value
 */
double ps_data_xvalue_at(PS_DATA data, int col);
double ps_data_yvalue_at(PS_DATA data, int row);
	
/** 
 * Initialize the data from a two-dimensional array of the correct size. 
 * Returns 0 upon success, error otherwise.
*/
int ps_data_init_with_array(PS_DATA data, double *values, int nrows, int ncols); //DEPRECATED, use ps_data_set_data
int ps_data_set_data(PS_DATA, double *values);

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


/**
 * Find the extreme values in the data 
 */
double ps_data_max_value(PS_DATA data);
double ps_data_min_value(PS_DATA data);
#endif