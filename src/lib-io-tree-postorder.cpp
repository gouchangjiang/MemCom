/*
 *  lib-io-tree-minmem.cpp
 *  lib-io-tree
 *
 *  Created by defbond on 8/10/10.
 *  Copyright 2010 LIP/ENS-Lyon. All rights reserved.
 *
 */
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <limits>


#include "lib-io-tree.h"
#include "lib-io-tree-postorder.h"

/* Postorder traversal*/
bool POOptimalSort(OrdoLiu_t a, OrdoLiu_t b){
  return ((a.max_pebble_cost-a.fi)>(b.max_pebble_cost-b.fi));
}

void PostOrderRecur(Cnode * current_root, OrdoLiu_t & SubSchedule, int quiet, int & count){
  if(!quiet){cerr<<"Treating T"<<current_root->GetId()<<endl;}

  count++;

  SubSchedule.schedule.clear();
  if (current_root->IsLeaf()){
    if(!quiet){cerr<<"Leaf found : "<<current_root->GetId()<<endl;}
    SubSchedule.schedule.push_back(current_root->GetId());
    SubSchedule.max_pebble_cost = current_root->GetCost();
    SubSchedule.fi = current_root->GetEW();
  }
  else{

    /* explore subtree to get their mem consumption*/
    list<OrdoLiu_t> * children_sched = new list<OrdoLiu_t>(current_root->GetChildren()->size());
    int i=0;
    for(list<OrdoLiu_t>::iterator cur_child = children_sched->begin();cur_child!=children_sched->end();cur_child++,i++){
      PostOrderRecur(current_root->GetChild(i), *cur_child, quiet, count);
    }


    /*sort the subtrees accordingly*/
    children_sched->sort(POOptimalSort);

    /*merge and add the root*/
    SubSchedule.max_pebble_cost = children_sched->front().max_pebble_cost;
    SubSchedule.fi = current_root->GetEW();


    double output_processed = 0;
    for(list<OrdoLiu_t>::iterator iter = children_sched->begin();iter!=children_sched->end();iter++){
      if(!quiet){cerr<<"This peak  "<< iter->max_pebble_cost<<" other consumption "<< output_processed <<endl;}
      /* merge two lists in constant time*/
      SubSchedule.schedule.splice(SubSchedule.schedule.end(),iter->schedule);
      SubSchedule.max_pebble_cost = max(SubSchedule.max_pebble_cost, iter->max_pebble_cost + output_processed);
      output_processed += iter->fi;
    }
    SubSchedule.max_pebble_cost = max(SubSchedule.max_pebble_cost,current_root->GetCost());

    if(!quiet){cerr<<"Max Pebble cost "<< SubSchedule.max_pebble_cost<<endl;}
    SubSchedule.schedule.push_back(current_root->GetId());

    delete children_sched;	
  }


  if(!quiet){cerr<<"End of treatment of T"<<current_root->GetId()<<endl;}
}

int ascend_po_comp(const void * a, const void * b){
  const po_node_t * pa = (po_node_t *)a;
  const po_node_t * pb = (po_node_t *)b;
  double diff = ((pb->w - pb->ew) - (pa->w - pa->ew));

  if(diff==0.0){
    return 0;
  }
  else if (diff>0.0){
    return 1;
  }
  else{
    return -1;
  }
}




double PostOrderIter(const int N,const int *prnts,const double *nwghts,const double *ewghts,const int * chstart, const int * chend, int * children, const int root, int *schedule){

  /*chstart and chend are indexed from 1 to N inclusive*/

  int * po = new int[N+1];
  memset ( (void *) po, 0, (N+1)*sizeof(*po) );

  int label =1;
  poaux(chstart, children, N, root, po, &label);



  double * w = new double[N+1];


  double Mr = 0;

#if VERBOSE	
  cerr<<"children list"<<endl;
  for(int ii = 1; ii<N+1;ii++){
    cerr<<children[ii]<<" ";
  }
  cerr<<endl;
#endif

  //TODO :change the size of scratchpad to MaxOutDegree
  po_node_t * scratchpad = new po_node_t[N+1];
  /*Post order is done, compute weight then sort children*/
  for(int i = 1;i<N+1;i++){
    int nd = po[i];
#if VERBOSE	
    cerr<<"processing node "<<nd<<endl;
#endif
    /*if nd is a leaf*/
    if(chstart[nd]==chstart[nd+1]){
      w[nd] = nwghts[nd] + ewghts[nd];
#if VERBOSE	
      cerr<<"node "<<nd<<" is a leaf of w="<<w[nd]<<endl;
#endif
    }
    else{
      int offset = 0;
      double chwghts = 0;
      double maxwghts = 0;
      int ch;
      for(int j = chstart[nd]; j<chstart[nd+1];j++){
        ch = children[j];
        chwghts += ewghts[ch];
        scratchpad[offset].w = w[ch];
        scratchpad[offset].ew = ewghts[ch];
        scratchpad[offset].index = ch;
        offset++;
      }
#if VERBOSE	
      cerr<<"unsorted list"<<endl;
      for(int ii = 0; ii<(chstart[nd+1]-chstart[nd]);ii++){
        cerr<<scratchpad[ii].index<<" ";
      }
      cerr<<endl;
#endif	


      qsort ( scratchpad, (chstart[nd+1]-chstart[nd]), sizeof(po_node_t), ascend_po_comp );
#if VERBOSE	
      cerr<<"sorted list"<<endl;
      for(int ii = 0; ii<(chstart[nd+1]-chstart[nd]);ii++){
        cerr<<scratchpad[ii].index<<" ";
      }
      cerr<<endl;
#endif


      double out_processed = 0;
      for(offset = 0;offset<(chstart[nd+1]-chstart[nd]);offset++){
#if VERBOSE	
        cerr<<"\tlook children "<< scratchpad[offset].index<<" ew = "<< scratchpad[offset].ew<<" peak = "<<scratchpad[offset].w<<endl;	
#endif
        children[chstart[nd]+offset] = scratchpad[offset].index;
        maxwghts = max(maxwghts,scratchpad[offset].w+out_processed);
        out_processed += /*ewghts[scratchpad[offset].index];*/scratchpad[offset].ew;
      }


      w[nd] = max(maxwghts,ewghts[nd]+nwghts[nd]+chwghts);
#if VERBOSE	
      cerr<<"node "<<nd<<" is not a leaf of w="<<w[nd]<<endl;
#endif
    }
#if VERBOSE	
    cerr<<endl;
#endif
    Mr = max(Mr,w[nd]);	
  }

#if VERBOSE	
  cerr<<"children list"<<endl;
  for(int ii = 1; ii<N+1;ii++){
    cerr<<children[ii]<<" ";
  }
  cerr<<endl;
#endif

  //	Mr = w[root];

  delete[] scratchpad;
  delete[] w;

  /*run PO again*/

  memset ( (void *) po, 0, (N+1)*sizeof(*po) );


  label =1;
  poaux(chstart, children, N, root, po, &label);

#if VERBOSE	
  cerr<<"children list"<<endl;
  for(int ii = 1; ii<N+1;ii++){
    cerr<<children[ii]<<" ";
  }
  cerr<<endl;
#endif

  for(int ii=1;ii<N+1;ii++){
    schedule[ii-1]= po[ii];
  }

  delete[] po;
  return Mr;

}



double PostOrderIterAlgorithm(int N, int *prnts, double *nwghts, double *ewghts, int *schedule){
  double Mr;
  //int count = 0;
  int * chstart,*chend,*children;
  int root;


  po_construct(N, prnts, &chstart,&chend,&children, &root);

  Mr = PostOrderIter(N, prnts, nwghts, ewghts, chstart,chend,children,root,schedule);

  delete[] chstart;
  delete[] chend;
  delete[] children;

  return Mr;
}

double PostOrderIterAlgorithm_timed(int N, int *prnts, double *nwghts, double *ewghts, int *schedule, double * usec,int quiet){
  double Mr;
  //int count = 0;
  int * chstart,*chend,*children;
  int root;


  po_construct(N, prnts, &chstart,&chend,&children, &root);

  *usec = -u_wseconds();
  Mr = PostOrderIter(N, prnts, nwghts, ewghts, chstart,chend,children,root,schedule);
  *usec += u_wseconds();

  delete[] chstart;
  delete[] chend;
  delete[] children;

  return Mr;
}


