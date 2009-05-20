
#include "ps_data.h"
#include "ps_errors.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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
	PS_DATA retval = ps_data_create(xsize, ysize);
}

PS_DATA ps_data_create(unsigned int xsize, unsigned int ysize) {
	PS_DATA new_data = (PS_DATA)malloc(sizeof(PS_DATA_T));
	new_data->xsize = xsize;
	new_data->ysize = ysize;
	new_data->x_values = (double*)malloc(sizeof(double)*xsize);
	new_data->y_values = (double*)malloc(sizeof(double)*ysize);
	new_data->data = (double*)malloc(xsize*ysize*sizeof(double));
	return new_data;
}

PS_DATA ps_data_copy(PS_DATA data) {
	PS_DATA newdata = ps_data_create(data->xsize, data->ysize);
	ps_data_set_x_values(newdata, data->x_values);
	ps_data_set_y_values(newdata, data->y_values);
	ps_data_init_with_array(newdata, data->data, data->ysize, data->xsize);
	return newdata;
}

void ps_data_destroy(PS_DATA data) {
	free(data->data);
	free(data->x_values);
	free(data->y_values);
	free(data);	
}

void ps_destroy_data(PS_DATA data) {
	ps_data_destroy(data);
}

int ps_data_rows(PS_DATA data) {
	return data->ysize;
}

int ps_data_columns(PS_DATA data) {
	return data->xsize;
}

int ps_data_init_with_array(PS_DATA data, double *values, int nrows, int ncols) {
	ps_data_set_data(data, values);
}

int ps_data_set_data(PS_DATA data, double *values) {
	memcpy(data->data, values, (data->xsize)*(data->ysize)*sizeof(double));	
}

int ps_data_set_value_at_row_column(PS_DATA data, double value, unsigned int row, unsigned int col) {
	int err = ps_data_check_index(data, row, col);
	if (PS_OK != err) {
		return err;
	}
	data->data[row*data->xsize + col] = value;
	return PS_OK;
}

int ps_data_set_x_values(PS_DATA data, double *x_values) {
	memcpy(data->x_values, x_values, (data->xsize)*sizeof(double));
	return PS_OK;
}

int ps_data_set_y_values(PS_DATA data, double *y_values) {
	memcpy(data->y_values, y_values, (data->ysize)*sizeof(double));
	return PS_OK;
}

int ps_data_set_x_value_at(PS_DATA data, int x, double x_value) {
	data->x_values[x] = x_value;
}

int ps_data_set_y_value_at(PS_DATA data, int y, double y_value) {
	data->y_values[y] = y_value;
}

double ps_data_xvalue_at(PS_DATA data, int col) {
	return data->x_values[col];
}

double ps_data_yvalue_at(PS_DATA data, int row) {
	return data->y_values[row];
}

int ps_data_value_at_row_column(PS_DATA data, unsigned int row, unsigned int col, double *val) {
	int err = ps_data_check_index(data, row, col);
	if (PS_OK != err) {
		return err;
	}
	*val = data->data[row*data->xsize + col];
	return PS_OK;
}

double ps_data_max_value(PS_DATA data) {
	int n = data->xsize*data->ysize;
	int i = 0;
	double max = data->data[i];
	for(i=1; i<n; i++) {
		if (data->data[i] > max) { max=data->data[i]; }
	}
	return max;
}

double ps_data_min_value(PS_DATA data) {
	int n = data->xsize*data->ysize;
	int i = 0;
	double min = data->data[i];
	for(i=1; i<n; i++) {
		if (data->data[i] < min) { min=data->data[i]; }
	}
	return min;
}