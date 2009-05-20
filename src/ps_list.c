#include "ps_list.h"
#include "ps_data.h"
#include <stdlib.h>

PS_LIST ps_list_create(void) {
	PS_LIST newlist = (PS_LIST)malloc(sizeof(PS_LIST_T));
	newlist->head = NULL;
	newlist->tail = NULL;
	newlist->current = NULL;
	return newlist;
}

/** Release the list */
void ps_list_destroy(PS_LIST list) {
	PS_LIST_NODE *node = list->head;
	if (node != NULL) {	
		PS_LIST_NODE *next = node->next;
		while (next != NULL) {
			free(node);
			node = next;
			next = node->next;
		}
		// At end of loop we still have to free the last node
		free(node);
	}
	free(list);
}

/** Release the list and the data */
void ps_list_destroy_all(PS_LIST list) {
	PS_LIST_NODE *node = list->head;
	if (node != NULL) {	
		PS_LIST_NODE *next = node->next;
		while (next != NULL) {
			ps_data_destroy(node->node_data->wavefunction);
			free(node->node_data);
			free(node);
			node = next;
			next = node->next;
		}
		// At end of loop we still have to free the last node
		ps_data_destroy(node->node_data->wavefunction);
		free(node->node_data);
		free(node);
	}
	free(list);	
}

/** Number of elements in list */
int ps_list_size(PS_LIST list) {
	int count = 0;
	PS_LIST_NODE *node = list->head;
	while (node != NULL) {
		count++;
		node = node->next;
	}
	return count;
}

/**
 * Set the current node to the head of the list
 */
PS_SOLUTION* ps_list_front(PS_LIST list) {
	list->current = list->head;
	if (list->current == NULL) {
		return NULL;
	} else {
		return list->current->node_data;		
	}
}

/**
 * Move to the next position in the list and return the PS_DATA object htere.
 */
PS_SOLUTION* ps_list_next(PS_LIST list) {
	if(list->current->next == NULL) {
		return NULL;
	} else {
		list->current = list->current->next;
		return list->current->node_data;
	}
}

/**
 * Add a new node onto the tail and make it current
 */
void ps_list_add(PS_LIST list, PS_SOLUTION *solution) {
	
	// Make a new node
	PS_LIST_NODE *newnode = (PS_LIST_NODE*)malloc(sizeof(PS_LIST_NODE));
	newnode->node_data = solution;
	newnode->next = NULL;
	
	if (list->head == NULL) {
		// Empty lists go onto the head
		list->head = newnode;
	} else {
		// Put onto the tail
		list->tail->next = newnode;	
	}
	list->tail = newnode;		
	list->current = newnode;
		
}

void ps_list_add_all(PS_LIST dest_list, PS_LIST src_list) {
	PS_SOLUTION *s = ps_list_front(src_list);
	while (s != NULL) { 
		ps_list_add(dest_list, s);
		s = ps_list_next(src_list);
	}
}