#ifndef PS_LIST_H
#define PS_LIST_H

#include "ps_data.h"
#include "psi_shooter.h"

typedef struct list_node {
	PS_SOLUTION *node_data;
	struct list_node *next;
} PS_LIST_NODE;

typedef struct {
	PS_LIST_NODE *head;
	PS_LIST_NODE *tail;
	PS_LIST_NODE *current;
} PS_LIST_T;

typedef PS_LIST_T* PS_LIST;

PS_LIST ps_list_create(void);
/** Release the list */
void ps_list_destroy(PS_LIST list);
/** Release the list and the data */
void ps_list_destroy_all(PS_LIST list);

/** List size */
int ps_list_size(PS_LIST list);

/** Add an element */
void ps_list_add(PS_LIST list, PS_SOLUTION *data);

/** Navigation returns NULL when empty or at end.*/
PS_SOLUTION* ps_list_front(PS_LIST list);
PS_SOLUTION* ps_list_next(PS_LIST list);

#endif