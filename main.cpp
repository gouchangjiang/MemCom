//
//  main.cpp
//  Heuristics
//
//  Created by changjiang GOU on 24/05/2017.
//  Copyright © 2017 Changjiang GOU. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "heuristics.hpp"
#include "lib-io-tree.h"
#include <math.h>

int main(int argc, const char * argv[]) {
    int tree_size=0;
    int *prnts;
    double *ewghts, *nwghts, *MSwghts;
    string dir=argv[1];
    string buffer;
    string filename;
    double memory;
    std::vector<double> memory_vec;
    std::vector<double> beta_vec;
    char cur_char;
    list<Cnode*> L;
    int * chstart,*chend,*children;
    int root;
    schedule_t schedule;
    double maxoutd;
    double minMem;
    uint64_t count=0;
    double beta;
    double ratio[] = {1/16,0.5,1,2,16};
    double MakeSpan;
    unsigned int com_freq=0;
    list<int>::iterator ite_sche;
    vector<Ctree*> subtrees;
    double comm_cost,makespan_chains_merged;
    unsigned int num_subtrees=0;
    double pro_per [] = {0.1,0.2,0.4};
    unsigned int num_p;
    
    cout.precision(20);
    
    cout<<"Tree_name"<<" "<<"Available_memory"<<" "<<"processors_perc"<<" "<<"CCR"<<" "<<"Communication_freq"<<" "<<"Makespan"<<" "<<"algorithm"<<endl;
    
    ifstream OpenFile(dir+argv[2]);
    do{
        OpenFile>>filename;
        parse_tree((dir+filename).c_str(), &tree_size, &prnts, &nwghts, &ewghts, &MSwghts);
        
        memory_vec.clear();
        OpenFile>>buffer;
        memory=stod(buffer);
        memory_vec.push_back(memory);
        OpenFile>>buffer;
        memory=stod(buffer);
        memory_vec.push_back(memory);
        
        do{
            OpenFile>>buffer;
            beta=stod(buffer);
            beta_vec.push_back(beta);
            cur_char = OpenFile.get();
        }while(cur_char != '\n' && OpenFile.good());
        
        po_construct(tree_size, prnts, &chstart,&chend,&children, &root);
        
        Ctree *treeobj = new Ctree(tree_size,prnts,nwghts,ewghts,MSwghts);
        
        schedule_t * schedule_f = new schedule_t();
        maxoutd = MaxOutDegree(treeobj, true);
        MinMem(treeobj, maxoutd, minMem, *schedule_f, true, count);
        
        /*cout<<"Schedule : {";
        for (list<int>::iterator ite=schedule_f->begin();ite!=schedule_f->end();ite++){ cout<<*ite<<" "; }
        cout<<"}"<<endl;*/
        
        int * schedule_copy = new int[tree_size+1];
        
        for (unsigned int memory_index=0; memory_index<2; ++memory_index) {
            for (unsigned int beta_index=0; beta_index<5; ++beta_index) {
                beta=beta_vec[beta_index];
                memory=memory_vec[memory_index];
                if (minMem<=memory) {//abundant memory
                    //MakeSpan=ASAP(treeobj, num_p, beta, num_subtrees,false,1);
                    //cout<<filename<<" "<<memory<<" "<<ratio[beta_vec.size()]<<" "<<num_subtrees-1<<" "<<MakeSpan<<" ASAP"<<endl;
                    
                    //MakeSpan=ASAP(treeobj, num_p, beta, num_subtrees,true,1);
                    //cout<<filename<<" "<<memory<<" "<<ratio[beta_vec.size()]<<" "<<num_subtrees-1<<" "<<MakeSpan<<" ASAPc"<<endl;
                    
                    //MakeSpan=splitSubtrees(treeobj, num_p, beta);
                    //cout<<filename<<" "<<memory<<" "<<ratio[beta_vec.size()]<<" "<<-1<<" "<<MakeSpan<<" split"<<endl;
                    
                    for (unsigned int index_p=0; index_p<3; ++index_p) {
                        num_p=ceil(pro_per[index_p]*tree_size);
                        MakeSpan=ASAP(treeobj, num_p, beta, num_subtrees,true,10);
                        cout<<filename<<" "<<memory<<" "<<pro_per[index_p]<<" "<<ratio[beta_index]<<" "<<num_subtrees-1<<" "<<MakeSpan<<" ASAPc10"<<endl;
                        
                        Qtree* qtree;
                        MakeSpan=AvoidChain(treeobj, num_subtrees, beta, makespan_chains_merged,qtree);
                        cout<<filename<<" "<<memory<<" "<<pro_per[index_p]<<" "<<ratio[beta_index]<<" "<<num_subtrees-1<<" "<<makespan_chains_merged<<" optASAPc10"<<endl;
                        
                        MakeSpan=LarSav(treeobj, qtree, num_subtrees, num_p, beta);//qtree will be deleted in  LarSav
                        cout<<filename<<" "<<memory<<" "<<pro_per[index_p]<<" "<<ratio[beta_index]<<" "<<num_subtrees-1<<" "<<MakeSpan<<" LarSav"<<endl;
                    }
                    
                    //MakeSpan=improvedSplit_v1(treeobj, num_p, beta, num_subtrees);
                    //cout<<filename<<" "<<memory<<" "<<1<<" "<<ratio[beta_vec.size()]<<" "<<-1<<" "<<MakeSpan<<" improve1"<<endl;
                    
                    //MakeSpan=improvedSplit_v2(treeobj, num_p, beta, num_subtrees);
                    //cout<<filename<<" "<<memory<<" "<<1<<" "<<ratio[beta_vec.size()]<<" "<<num_subtrees-1<<" "<<MakeSpan<<" improve2"<<endl;
                    
                }else{//memory constraint case
                    ite_sche=schedule_f->begin();
                    for (unsigned int i=tree_size; i>=1; --i) {
                        schedule_copy[i]=*ite_sche;
                        advance(ite_sche, 1);
                    }
                    schedule_copy[0]=tree_size+1;
                    
                    com_freq=0;
                    subtrees.clear();
                    comm_cost = IOCounter(treeobj,tree_size+1, nwghts, ewghts,chstart,children, schedule_copy, memory, false, true,com_freq,FIRST_FIT);
                    Qtree* qtree = new Qtree(treeobj->GetRoot(),com_freq+1,beta,subtrees);
                    cout<<filename<<" "<<memory<<" "<<"NA"<<" "<<ratio[beta_index]<<" "<<com_freq<<" "<<qtree->GetMakeSpan(com_freq+1)<<" firstfit"<<endl;
                    MakeSpan = Upper(treeobj, qtree, &subtrees, com_freq+1, memory, beta);
                    cout<<filename<<" "<<memory<<" "<<"NA"<<" "<<ratio[beta_index]<<" "<<com_freq<<" "<<MakeSpan<<" "<<" firstfit+upper"<<endl;
                    delete qtree;
                    
                    ite_sche=schedule_f->begin();
                    for (unsigned int i=tree_size; i>=1; --i) {
                        schedule_copy[i]=*ite_sche;
                        advance(ite_sche, 1);
                    }
                    schedule_copy[0]=tree_size+1;
                    com_freq=0;
                    subtrees.clear();
                    comm_cost = IOCounter(treeobj,tree_size+1, nwghts, ewghts,chstart,children, schedule_copy, memory, false, true,com_freq,LARGEST_FIT);
                    qtree = new Qtree(treeobj->GetRoot(),com_freq+1,beta,subtrees);
                    cout<<filename<<" "<<memory<<" "<<"NA"<<" "<<ratio[beta_index]<<" "<<com_freq<<" "<<qtree->GetMakeSpan(com_freq+1)<<" largestfirst"<<endl;
                    MakeSpan = Upper(treeobj, qtree, &subtrees, com_freq+1, memory, beta);
                    cout<<filename<<" "<<memory<<" "<<"NA"<<" "<<ratio[beta_index]<<" "<<com_freq<<" "<<MakeSpan<<" "<<" largestfirst+upper"<<endl;
                    delete qtree;
                    
                    ite_sche=schedule_f->begin();
                    for (unsigned int i=tree_size; i>=1; --i) {
                        schedule_copy[i]=*ite_sche;
                        advance(ite_sche, 1);
                    }
                    schedule_copy[0]=tree_size+1;
                    com_freq=0;
                    subtrees.clear();
                    Immediately(treeobj, tree_size+1, nwghts, ewghts, chstart, children, schedule_copy, memory, com_freq);
                    qtree = new Qtree(treeobj->GetRoot(),com_freq+1,beta,subtrees);
                    cout<<filename<<" "<<memory<<" "<<"NA"<<" "<<ratio[beta_index]<<" "<<com_freq<<" "<<qtree->GetMakeSpan(com_freq+1)<<" immediately"<<endl;
                    MakeSpan = Upper(treeobj, qtree, &subtrees, com_freq+1, memory, beta);
                    cout<<filename<<" "<<memory<<" "<<"NA"<<" "<<ratio[beta_index]<<" "<<com_freq<<" "<<MakeSpan<<" "<<" immediately+upper"<<endl;
                    delete qtree;
                    
                    //com_freq++;
                    //MakeSpan = LarSav(treeobj, qtree, com_freq, num_p, beta);//qtree will be deleted in LarSav
                    //cout<<filename<<" "<<memory<<" "<<com_freq<<" "<<MakeSpan<<" "<<" larsav"<<endl;
                }
            }
        }
        
        delete treeobj;
        delete schedule_f;
        delete[] schedule_copy;
        delete[] prnts;
        delete[] ewghts;
        delete[] nwghts;
        delete[] MSwghts;
        delete[] chstart;
        delete[] chend;
        delete[] children;
    }while (OpenFile.good());
    OpenFile.close();
    
    return 0;
}

