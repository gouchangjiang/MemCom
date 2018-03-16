//
//  cut.hpp
//  MinCom
//
//  Created by Changjiang on 30/12/2016.
//  Copyright © 2016 Changjiang. All rights reserved.
//

#ifndef heuristics_hpp
#define heuristics_hpp

#include <iostream>
#include <forward_list>
#include <unordered_map>
#include "lib-io-tree.h"
#include "lib-io-tree-minmem.h"

vector<unsigned int> Greedy(Ctree *treeobj, double m_avail);
vector<unsigned int> combine_sim(Ctree *treeobj,double m_avail, list<Cnode*> &L,schedule_t &schedule);
vector<unsigned int> combine(Ctree *treeobj,double m_avail, list<Cnode*> &L,schedule_t &schedule);

double splitSubtrees(Ctree *tree, unsigned int num_p, double beta);
double improvedSplit_v1(Ctree* tree,unsigned int num_p, double beta, unsigned int & num_subtrees);
double improvedSplit_v2(Ctree* tree,unsigned int num_p, double beta, unsigned int & num_subtrees);
double ASAP(Ctree* tree, unsigned int num_p, double beta, unsigned int &num_subtrees, bool consider_commu,unsigned int deepth);
double AvoidChain(Ctree* tree, unsigned int &num_subtrees,double beta, double & after_merged, Qtree* & out_qtree);
double Upper(Ctree* tree, Qtree* qtree, vector<Ctree*>* subtrees, unsigned int num_subtrees, double memory, double beta);
double LarSav(Ctree* tree, Qtree* qtree, unsigned int & num_subtrees, unsigned int num_p, double beta);
void Immediately(Ctree* tree, unsigned long N, double * nwghts, double * ewghts, int * chstart, int * children, int * schedule, double m_availble, unsigned int & num_para_subtrees);

bool maxoutEqualMinmem(Ctree* tree,double &maxout,double &minMem);
#endif /* heuristics_hpp */
