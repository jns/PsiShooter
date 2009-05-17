#include "ps_data.h"
#include "ps_data_io.h"
#include "ps_list.h"
#include "psi_shooter.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {

	int i;
	
	if (2 != argc) {
		printf("test_ls_list infile\n");
		return 1;
	}

	char* infname = argv[1];
	FILE *infile = fopen(infname, "r");
	PS_DATA data = ps_data_read_bin(infile);
	fclose(infile);
	
	PS_SOLUTION *solution;
	
	PS_LIST list = ps_list_create();
	for(i=0; i<5; i++) {
		solution = (PS_SOLUTION*)malloc(sizeof(PS_SOLUTION));
		solution->energy = i;
		solution->wavefunction = ps_data_copy(data);
		ps_list_add(list, solution);
		printf("List Size=%i\n", ps_list_size(list));
	}
	free(data);

	PS_SOLUTION *listnode = ps_list_front(list);
	i = 0;
	while (listnode != NULL) {
		printf("Node %i\n", ++i);
		printf("Energy=%f\n", listnode->energy);
		ps_data_write_ascii(listnode->wavefunction, stdout);
		listnode = ps_list_next(list);
	}
	
	ps_list_destroy_all(list);
	return 0;
}	