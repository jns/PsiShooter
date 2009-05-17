#include "ps_data.h"
#include "ps_data_io.h"
#include "ps_list.h"
#include <stdio.h>

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
	
	PS_LIST list = ps_list_create();
	for(i=0; i<5; i++) {
		ps_list_add(list, data);
		printf("List Size=%i\n", ps_list_size(list));
	}

	PS_DATA node = ps_list_front(list);
	i = 0;
	while (node != NULL) {
		printf("Node %i\n", ++i);
		ps_data_write_ascii(node, stdout);
		node = ps_list_next(list);
	}
	return 0;
}	