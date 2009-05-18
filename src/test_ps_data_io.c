#include "ps_data_io.h"
#include <stdio.h>

int main(int argc, char **argv) {
	return read_ascii(argc, argv);
}

int create_write(int argc, char **argv) {

	if (2 != argc) {
		printf("test_ps_data_io outfile");
		return 1;
	}

	char *outfname = argv[1];
	
	int rows = 10; 
	int cols = 10;
	double xvalues[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
	double yvalues[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
	int i,j;
	
	PS_DATA data = ps_data_create(rows, cols);
	ps_data_set_x_values(data, xvalues);
	ps_data_set_y_values(data, yvalues);
	for (i=0; i < rows; i++) {
		for (j=0; j< cols; j++) {
			ps_data_set_value_at_row_column(data, 5.0, i, j);
		}
	}
	
	FILE *outfile = fopen(outfname, "w");
	ps_data_write_bin(data, outfile);
	fclose(outfile);
	ps_data_destroy(data);
	printf("Wrote %s\n", outfname);
	return 0;
}

int read_print(int argc, char **argv) {

	if (2 != argc) {
		printf("test_ps_data_io infile");
		return 1;
	}

	char *infname = argv[1];
	FILE *infile = fopen(infname, "r");
	PS_DATA data = ps_data_read_bin(infile);
	fclose(infile);
	
	ps_data_write_ascii(data, stdout);
	ps_data_destroy(data);
	return 0;
}

int read_ascii(int argc, char **argv) {
	int i,j;
	if (2 != argc) {
		printf("read_ascii infile");
		return 1;
	}
	
	char *infname = argv[1];
	FILE *infile = fopen(infname, "r");
	PS_DATA data = ps_data_read_ascii(infile);
	fclose(infile);
	
	ps_data_write_ascii(data, stdout);
	ps_data_destroy(data);
	return 0;
}

int read_write(int argc, char **argv) {
	int i,j;
	
	if (3 != argc) {
		printf("test_ps_data infile outfile\n");
		return 1;
	}
	
	char *infname = argv[1];
	char *outfname = argv[2];
	FILE *infile = fopen(infname, "r");
	PS_DATA data = ps_data_read_bin(infile);
	fclose(infile);
	
	if (NULL == data) {
		printf("Error!\n");
		return 1;
	}

	int rows = ps_data_rows(data);
	int cols = ps_data_columns(data);
	double v;
	
	printf("Read %s\n", infname);
	printf("rows=%i , cols=%i\n", rows, cols);
	printf("X Values:\n");
	for (i = 0; i < rows; i++) {
		printf("%e\n", ps_data_xvalue_at(data, i));
	}
	
	
	FILE *outfile = fopen(outfname, "w");
	ps_data_write_bin(data, outfile);
	fclose(outfile);
	
	printf("Wrote %s\n", outfname);
	ps_data_destroy(data);
	
	// for (i=0; i < rows; i++) {
	// 	for (j=0; j< cols; j++) {
	// 		v = ps_data_value(data, i, j);
	// 		printf("%f", i, j, v );				
	// 		if (cols-1 == j) {
	// 			printf("\n");
	// 		} else {
	// 			printf(" ");
	// 		}
	// 	}
	// }
	return 0;
}	