
#include "ps_data.h"
#include "ps_errors.h"
#include <stdlib.h>

// Internal function to check index
int ps_data_check_index(PS_DATA data, int row, int col) {
	if (row < 0 || row >= data->ysize) {
		return PS_ERROR_INDEX_OUT_OF_BOUNDS;
	}

	if (col < 0 || col >= data->xsize) {
		return PS_ERROR_INDEX_OUT_OF_BOUNDS;
	}	
	return PS_OK;
}

PS_DATA ps_create_data(unsigned int xsize, unsigned int ysize, double xstep, double ystep) {
	PS_DATA new_data = (PS_DATA)malloc(sizeof(PS_DATA_T));
	new_data->xsize = xsize;
	new_data->ysize = ysize;
	new_data->xstep = xstep;
	new_data->ystep = ystep;
	new_data->data = (double*)malloc(xsize*ysize*sizeof(double));
}

void ps_destroy_data(PS_DATA data) {
	free(data->data);
	free(data);
}

int ps_data_rows(PS_DATA data) {
	return data->xsize;
}

int ps_data_columns(PS_DATA data) {
	return data->ysize;
}

int ps_data_init_with_array(PS_DATA data, double *values, int nrows, int ncols) {
	int i,j;
	if (nrows == data->xsize && ncols == data->ysize) {
		for (i = 0; i < nrows; i++) {
			for (j=0; j< ncols; j++) {
				data->data[i*ncols + j] = values[i*ncols + j];
			}
		}
		return PS_OK;
	} else {
		return PS_ERROR_DIMENSION_MISMATCH;
	}
	
}

int ps_data_set_value_at_row_column(PS_DATA data, double value, unsigned int row, unsigned int col) {
	int err = ps_data_check_index(data, row, col);
	if (PS_OK != err) {
		return err;
	}
	data->data[row*data->xsize + col] = value;
	return PS_OK;
}

int ps_data_value_at_row_column(PS_DATA data, unsigned int row, unsigned int col, double *val) {
	int err = ps_data_check_index(data, row, col);
	if (PS_OK != err) {
		return err;
	}
	*val = data->data[row*data->xsize + col];
	return PS_OK;
}

