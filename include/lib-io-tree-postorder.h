/*
 *  lib-io-tree-minmem.h
 *  lib-io-tree
 *
 *  Created by defbond on 8/10/10.
 *  Copyright 2010 LIP/ENS-Lyon. All rights reserved.
 *
 */
#ifndef LIB_IO_TREE_POSTORDER_H
#define LIB_IO_TREE_POSTORDER_H

#ifdef __cplusplus


#include <list>
#include "lib-io-tree-utils.h"

struct po_node_t {
	int index;
	double ew;
	double w;
};

bool POOptimalSort(OrdoLiu_t a, OrdoLiu_t b);
int ascend_po_comp(const void * a, const void * b);

void PostOrderRecur(Cnode * current_root, OrdoLiu_t & SubSchedule, int quiet, int & count);
double PostOrderIter(const int N,const int *prnts,const double *nwghts,const double *ewghts,const int * chstart, const int * chend, int * children, const int root, int *schedule);
	
#endif

#endif
