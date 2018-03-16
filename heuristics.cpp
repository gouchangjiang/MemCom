//
//  cut.cpp
//  MinCom
//
//  Created by Changjiang on 30/12/2016.
//  Copyright © 2016 Changjiang. All rights reserved.
//

#include "heuristics.hpp"
#include <time.h>

//comparing function used for finding the minimum edge in a tree
bool cmp(const Cnode *a,const Cnode *b){
    return a->GetEW()<b->GetEW();
}

//comparing function used for sorting subtrees
bool cmp_s(Ctree *a, Ctree *b){
    return (a->GetRoot()->GetEW())>(b->GetRoot()->GetEW());
}

//finding the node with minimum sum of outputfiles and execution file
bool cmp_oe(const Cnode *a, const Cnode *b){
    return a->GetEO()<b->GetEO();
}

vector<unsigned int> greedyC(Ctree *tree,vector<Ctree*> &subtrees,double m_avail, list<Cnode*> &L, schedule_t &schedule, bool skipCheckM){
    //cout<<"-----------------In GreedyC--------------------\n";
    //cout<<"   "<<tree->GetNodes()->size()-1<<" nodes in tree "<<tree->GetTreeId()<<".\n";
    uint64_t count = 0;
    vector<unsigned int> cut_set;
    const vector<Cnode*> *nodes;
    double memory_peak=m_avail;
    double M = numeric_limits<double>::infinity();
    
    //double Mp = MaxOutDegree(tree, 1);//1 means quiet
    
    /*if(Mp>m_avail){
        //cerr<<"Error: available memory is too small. "<<endl;
        //cout<<"Set it to 1.013 times "<<Mp<<'\n'<<endl;
        m_avail=1.0125*Mp;
        //exit(EXIT_FAILURE);
    }*/
    
    //clock_t t;
    if(skipCheckM==false){
        //t=clock();
        //MinMem(tree, Mp, min_memory, *schedule, true, count);
        L.clear();
        schedule.clear();
        explore(tree->GetRoot(), m_avail, &L, &schedule, M, L, schedule, memory_peak, 1, 0,count);
        //printf ("   Memory_peak is %g , it took me %lu clicks.\n",memory_peak,clock()-t);
    }//else{cout<<"   skip\n";}
    
    if (skipCheckM==true||!L.empty()) {
        nodes = tree->GetNodes();
        //t=clock();
        vector<Cnode*>::const_iterator min_edge = min_element(nodes->begin()+1,nodes->end()-1,cmp);
        
        //check
        list<int>::iterator find_in_Sche;
        find_in_Sche=find(schedule.begin(),schedule.end(),(*min_edge)->GetId());
        
        if (find_in_Sche!=schedule.end()) {
            if ((*min_edge)->IsLeaf()) {
                skipCheckM=true;
            }else{
                skipCheckM=false;
            }
        }else{
            bool roll_back=true;
            for (list<Cnode*>::iterator ite=L.begin(); ite!=L.end(); ++ite) {
                if ((*ite)->GetEO()+M==memory_peak) {
                    roll_back=false;
                    break;
                }
            }
            if (roll_back==true) {
                skipCheckM=false;
            }else{
                list<Cnode*>::iterator find_in_L;
                find_in_L=find(L.begin(), L.end(), *min_edge);
                if (find_in_L!=L.end()) {
                    skipCheckM=false;
                }else{
                    find_in_L=find(L.begin(),L.end(), (*min_edge)->GetParent());
                    if (find_in_L!=L.end()) {
                        skipCheckM=false;
                    }else{skipCheckM=true;}
                }
            }
        }
        
        //Cnode *node_t=(*min_edge);
        //printf("   Min_edge is %u, edge weight is %g\n",(*min_edge)->GetId(),node_t->GetEW());
        //printf ("   Search min_edge took me %lu clicks.\n",clock()-t);
        //cout<<"-cut---";
        cut_set.push_back((*min_edge)->GetId());
        //cout<<(*min_edge)->GetId()<<'\n'<<endl;
        Ctree *subtree = new Ctree();
        //t=clock();
        subtree = tree->Cut(*min_edge);
        //printf ("   Cut it took me %lu clicks.\n",clock()-t);
        
        
        //skipCheckM=false;
        
        /*if (L.size()==1) {
            Cnode *root_node = L.front();
            bool only_child = root_node->OnlyOneChild();
            while (only_child) {
                root_node=root_node->GetChildren()->front();
                only_child = root_node->OnlyOneChild();
            }
            
            vector<Cnode*>::const_iterator ite=tree->GetNodes()->begin();
            while ((*ite)!=root_node) {
                advance(ite,1);
            }
            tree->SetRootId(ite+1-tree->GetNodes()->begin());//change the root id
        }*/
        
        schedule_t schedule_sub;
        list<Cnode*> L_sub;
        L_sub.clear();
        schedule_sub.clear();
        vector<unsigned int> cut_i = greedyC(subtree, subtrees, m_avail,L_sub,schedule_sub,false);
    
        vector<unsigned int> cut_r = greedyC(tree, subtrees, m_avail,L,schedule,skipCheckM);//remained part
        
        cut_set.insert(cut_set.begin(), cut_i.begin(), cut_i.end());
        cut_set.insert(cut_set.begin(), cut_r.begin(), cut_r.end());

        return cut_set;
    } else {
        subtrees.push_back(tree);
        return cut_set;
    }
}


void greedyR(vector<unsigned int> &cut,vector<Ctree*> &subtrees,double m_avail){
    uint64_t count = 0;
    schedule_t schedule;
    list<Cnode*> L;
    double memory_peak=m_avail;
    double M = numeric_limits<double>::infinity();
    
    sort(subtrees.begin(),subtrees.end(),cmp_s);
    vector<Ctree*>::iterator father_tree;
    for (vector<Ctree*>::iterator iter=subtrees.begin(); iter!=subtrees.end(); ++iter) {
        Cnode *root_node = (*iter)->GetNodeByPos(0);
        Cnode *parent_tree_root = root_node->GetParent();
        if (parent_tree_root!=NULL) {
            while (!parent_tree_root->IsRoot()) {
                parent_tree_root=parent_tree_root->GetParent();
            }
            
            father_tree = subtrees.begin();
            while ((*father_tree)->GetTreeId()!=parent_tree_root->GetId()) {
                ++father_tree;
            }// find father tree
            
            (*father_tree)->Merge(*iter, root_node->GetParent());
            //double Mp = MaxOutDegree(*father_tree, 1);
            //MinMem(*father_tree, Mp, min_memory, *schedule, true, count);
            L.clear();
            schedule.clear();
            explore((*father_tree)->GetRoot(), m_avail, &L, &schedule, M, L, schedule, memory_peak, 1, 0,count);
            
            
            if (L.empty()) {
                vector<unsigned int>::const_iterator position=find(cut.begin(), cut.end(), root_node->GetId());
                cut.erase(position);
                (*iter)->Clear();
                //(*iter)->GetRoot()->SetParentId(parent_tree_root->GetId());//equal to deleting tree in subtrees set
                //cout<<"-restore-"<<endl;
            } else {
                (*father_tree)->Cut(root_node);
            }
        }
    }
}

vector<unsigned int> Greedy(Ctree *treeobj, double m_avail){
    vector<unsigned int> cut;
    vector<Ctree*> subtrees;
    schedule_t schedule;
    list<Cnode*> L;
    schedule.clear();
    L.clear();
    
    cut = greedyC(treeobj, subtrees, m_avail, L, schedule, false);
    
    if (cut.size()>1) {
        greedyR(cut, subtrees, m_avail);
    }
    for(vector<Ctree*>::iterator iter = subtrees.begin();iter!=subtrees.end();iter++){
        delete *iter;
    }
    return cut;
}

void exact_subset_sum(list<Cnode*> *set,unsigned long target,list<Cnode*> &keep, list<Cnode*> &cut){
    map<unsigned long int, Cnode*> path;
    list<unsigned long int> set_first(1,0);
    list<unsigned long int> set_second;
    list<Cnode*>::iterator first_ite=set->begin();
    for (list<Cnode*>::iterator ite_out=first_ite; ite_out!=set->end(); ++ite_out) {
        for (list<unsigned long int>::iterator ite_in=set_first.begin(); ite_in!=set_first.end(); ++ite_in) {
            if ((*ite_in+(*ite_out)->GetEW())<=target) {
                set_second.push_back(*ite_in+(*ite_out)->GetEW());
                path.insert(pair<unsigned long int, Cnode*>((*ite_in+(unsigned long int)(*ite_out)->GetEW()),*ite_out));
            }
        }
        set_first.sort();
        set_second.sort();
        set_first.merge(set_second);
        set_first.unique();
    }
    
    unsigned long int subset_sum=set_first.back();
    
    map<unsigned long int, Cnode*>::iterator it;
    while (subset_sum!=0) {
        it = path.find(subset_sum);
        keep.push_back(it->second);
        subset_sum = subset_sum - it->second->GetEW();
        path.erase(it);
    }
    
    cut=*set;
    for (list<Cnode*>::iterator it=keep.begin(); it!=keep.end(); ++it) {
        cut.remove(*it);
    }

    if (cut.empty()) {//if cut is empty, just return false, redo it with a new target
        /*list<Cnode*>::iterator ite = keep.begin();
        advance(ite,1);
        list<Cnode*>::const_iterator min_edge = min_element(ite,keep.end(),cmp_oe);
        cut.push_back(*min_edge);
        keep.erase(min_edge);*/
    }
}

vector<unsigned int> combine_sim(Ctree *treeobj,double m_avail, list<Cnode*> &L,schedule_t &schedule){
    double memory_peak=m_avail;
    uint64_t count=0;
    vector<unsigned int> cut_set;
    vector<unsigned int> cut_set_sub;
    list<Cnode*> keep_set;
    vector<Ctree*> subtrees;
    list<Cnode*> cut;
    unsigned long target;
    double M = numeric_limits<double>::infinity();
    
    if (!L.empty()) {
        M=0;
        for (list<Cnode*>::iterator it=L.begin(); it!=L.end(); ++it) {
            M=M+(*it)->GetEW();
        }
    }
    
    //cout<<"m_avail in explore is: "<<m_avail<<endl;
    //cout<<"The root of the tree is "<<treeobj->GetRoot()->GetId()<<" , it's nw is "<<treeobj->GetRoot()->GetNW()<<endl;
    explore(treeobj->GetRoot(), m_avail, &L, &schedule, M, L, schedule, memory_peak, 1, 0,count);
    //cout<<"The size of L is: "<<L.size()<<". They are: "<<endl;
    /*for (list<Cnode*>::iterator ite=L.begin(); ite!=L.end(); ++ite) {
        cout<<(*ite)->GetId()<<" ";
    }*/
    //cout<<"\n";
    //cout<<"The size of Schedule is: "<<schedule.size()<<". They are:"<<endl;
    /*for (list<int>::iterator ite=schedule.begin(); ite!=schedule.end(); ++ite) {
        cout<<*ite;
        if (!treeobj->GetNode(*ite)->IsLeaf()) {
            cout<<"*";
        }
        cout<<" ";
    }*/
    //cout<<"M is "<<M<<endl;
    //cout<<"memory_peak is "<<memory_peak<<endl;
    
    if (L.empty()) {
        delete treeobj;
        return cut_set;
    }else if (L.size()>1){
        if (L.size()==2) {
            //cout<<"Size of L: 2"<<endl;
            if(L.front()->GetEW() > L.back()->GetEW()){
                keep_set.push_back(L.front());
                cut.push_back(L.back());
            }else{
                keep_set.push_back(L.back());
                cut.push_back(L.front());
            };
        }else{
            //cout<<"Size of L: "<<L.size()<<endl;
            list<Cnode*>::const_iterator min_edge = min_element(L.begin(),L.end(),cmp_oe);
            keep_set.push_back(*min_edge);
            target=m_avail-(*min_edge)->GetCost();
            L.erase(min_edge);
            exact_subset_sum(&L, target, keep_set,cut);
            }
        if (cut.empty()) {//use the peak memory as target
            keep_set.clear();
            exact_subset_sum(&L, memory_peak, keep_set, cut);
        }
        
        for (list<Cnode*>::iterator ite=cut.begin(); ite!=cut.end(); ++ite) {
            cut_set.push_back((*ite)->GetId());
            subtrees.push_back(treeobj->Cut(*ite));
        }
        
        cut_set_sub=combine_sim(treeobj, m_avail, keep_set, schedule);
        cut_set.insert(cut_set.end(), cut_set_sub.begin(), cut_set_sub.end());
        
        for (vector<Ctree*>::iterator ite=subtrees.begin(); ite!=subtrees.end(); ++ite) {
            L.clear();
            schedule.clear();
            cut_set_sub=combine_sim(*ite, m_avail, L,schedule);
            cut_set.insert(cut_set.end(), cut_set_sub.begin(), cut_set_sub.end());
        }
        return cut_set;

    }else{
        Cnode *root_node = L.front();
        bool only_child = root_node->OnlyOneChild();
        while (only_child) {
            root_node=root_node->GetChildren()->front();
            only_child = root_node->OnlyOneChild();
        }
        
        vector<Cnode*>::const_iterator ite=treeobj->GetNodes()->begin();
        while (*ite!=root_node) {
            ++ite;
        }
        treeobj->SetRootId(ite+1-treeobj->GetNodes()->begin());//change the root id
        L.clear();
        schedule.clear();
        //cout<<"the root of treeobj is: "<<treeobj->GetRootId()<<endl;
        //double maxout1,minMem1;
        //bool temp = maxoutEqualMinmem(treeobj,maxout1,minMem1);
        //cout<<"maxout: "<<maxout1<<"  minMem: "<<minMem1<<endl;
        cut_set_sub=combine_sim(treeobj, m_avail, L, schedule);
        cut_set.insert(cut_set.end(), cut_set_sub.begin(), cut_set_sub.end());
        return cut_set;
    }
}

vector<unsigned int> combine(Ctree *treeobj,double m_avail, list<Cnode*> &L,schedule_t &schedule){
	double memory_peak=m_avail;
    uint64_t count=0;
    vector<unsigned int> cut_set;
    vector<unsigned int> cut_set_sub;
    list<Cnode*> keep_set;
    vector<Ctree*> subtrees;
    list<Cnode*> cut;
    unsigned long target;
    double M = numeric_limits<double>::infinity();
    
    if (!L.empty()) {
        M=0;
        for (list<Cnode*>::iterator it=L.begin(); it!=L.end(); ++it) {
            M=M+(*it)->GetEW();
        }
    }
    
    //cout<<"Tree "<<treeobj->GetTreeId()<<endl;
    //cout<<"   Enter Explore.\n";
    explore(treeobj->GetRoot(), m_avail, &L, &schedule, M, L, schedule, memory_peak, 1, 0,count);

    if (L.empty()) {
        //cout<<"   Out Explore, L is empty.\n";
        delete treeobj;
        return cut_set;
    }else if (L.size()>1){
        /*cout<<"   Out Explore, The size of L is "<<L.size()<<".\n";
        cout<<"   ";
        for (list<Cnode*>::iterator ite=L.begin(); ite!=L.end(); ++ite) {
            cout<<(*ite)->GetId()<<"|"<<(*ite)->GetEW()<<"|"<<(*ite)->GetNW()<<"|"<<(*ite)->GetCost()<<"   ";
        }*/
        //cout<<endl;
        list<Cnode*> copy_keep_set;
        list<Cnode*> keep;
        L.sort(cmp_oe);
        unsigned long edges_sum;
        unsigned int len_copy;
        
        do {
            if (L.empty()) {//combine fail, use combine simple
                L.splice(L.begin(), keep);
                target=L.front()->GetCost();
                keep_set.clear();
                cut.clear();
                exact_subset_sum(&L, target, keep_set,cut);
                break;
            }
            keep.push_back(L.front());
            keep_set.clear();
            keep_set=keep;
            
            target=m_avail-L.front()->GetCost();
            L.pop_front();
            cut.clear();
            exact_subset_sum(&L, target, keep_set,cut);
            if (cut.empty()) {//use peak memory as target
                keep_set.clear();
                exact_subset_sum(&L, memory_peak, keep_set, cut);
            }
            edges_sum=0;
            len_copy=0;
            for (list<Cnode*>::iterator ite=keep_set.begin(); ite!=keep_set.end(); ++ite) {
                edges_sum=edges_sum+(*ite)->GetEW();
            }
            
            copy_keep_set=keep_set;
            while (len_copy!=copy_keep_set.size()) {
                len_copy=copy_keep_set.size();
                for (list<Cnode*>::iterator ite=copy_keep_set.begin(); ite!=copy_keep_set.end(); /*Empty on purpose*/) {
                    if ((*ite)->GetEO()<=m_avail-edges_sum) {
                        edges_sum=edges_sum-(*ite)->GetEW();
                        ite=copy_keep_set.erase(ite);
                    }else{
                        ++ite;
                    }
                }
            }
        } while (edges_sum!=0);

        //cout<<"   Edges cut are: ";
        for (list<Cnode*>::iterator ite=cut.begin(); ite!=cut.end(); ++ite) {
            //cout<<(*ite)->GetId()<<" ";
            cut_set.push_back((*ite)->GetId());
            subtrees.push_back(treeobj->Cut(*ite));
        }
        //cout<<"\n";

        cut_set_sub=combine(treeobj, m_avail, keep_set, schedule);
        cut_set.insert(cut_set.end(), cut_set_sub.begin(), cut_set_sub.end());
        
        for (vector<Ctree*>::iterator ite=subtrees.begin(); ite!=subtrees.end(); ++ite) {
            L.clear();
            schedule.clear();
            cut_set_sub=combine(*ite, m_avail, L,schedule);
            cut_set.insert(cut_set.end(), cut_set_sub.begin(), cut_set_sub.end());
        }
        return cut_set;
    }else{
        //cout<<"   Out Explore, 1 node in L.\n";
        Cnode *root_node = L.front();
        bool only_child = root_node->OnlyOneChild();
        while (only_child) {
            root_node=root_node->GetChildren()->front();
            only_child = root_node->OnlyOneChild();
        }
        
        vector<Cnode*>::const_iterator ite=treeobj->GetNodes()->begin();
        while (*ite!=root_node) {
            ++ite;
        }
        treeobj->SetRootId(ite+1-treeobj->GetNodes()->begin());
        L.clear();
        schedule.clear();
        cut_set_sub=combine(treeobj, m_avail, L, schedule);
        cut_set.insert(cut_set.end(), cut_set_sub.begin(), cut_set_sub.end());
        return cut_set;
    }
}

bool maxoutEqualMinmem(Ctree* tree,double &maxout,double &minMem){
    uint64_t count = 0;
    schedule_t * schedule = new schedule_t();
    maxout = MaxOutDegree(tree, 1);
    MinMem(tree, maxout, minMem, *schedule, true, count);
    delete schedule;
    if (minMem==maxout)
    {return true;}else{return false;}
}


//comparing function used for sorting by non-decreasing Makespan with communication
bool cmp_PQ(Cnode *a,Cnode *b){
    return a->GetMSCCost()<b->GetMSCCost();
}

//non-decreasing without conmmunication
bool cmp_more(Cnode *a, Cnode *b){
    return a->GetMSCost()<b->GetMSCost();
}

double splitSubtrees(Ctree *tree, unsigned int num_p, double beta){
    if (beta<=0) {
        return -1; //beta should be positive
    }
    
    vector<Cnode*> PQ(1,tree->GetRoot());
    vector<double> MS(1,tree->GetRoot()->GetMSCost(beta));
    vector<Cnode*> SeqSet;
    Cnode* node;
    
    double MS_seqSet=0;
    double WMore,W_PQ;
    while (!PQ.back()->IsLeaf()) {//only when head(PQ) is a leaf node, W equals to w
        node = PQ.back();
        SeqSet.push_back(node);
        MS_seqSet = MS_seqSet+node->GetMSW();
        PQ.pop_back();
        for (vector<Cnode*>::iterator iter=node->GetChildren()->begin();iter!=node->GetChildren()->end();iter++) {
            PQ.push_back(*iter);
        }
        
        WMore=0;
        if (PQ.size()>=num_p) {
            sort(PQ.begin(), PQ.end(), cmp_more);//no communication
            for (unsigned int i=0;i<(PQ.size()-num_p+1);++i) {
                WMore+=PQ[i]->GetMSCost();// no communication
            }
        }
        
        sort(PQ.begin(), PQ.end(),cmp_PQ);//non-decreasing sort
        W_PQ=PQ.back()->GetMSCCost();
        //cout<<"seqSet: "<<MS_seqSet<<" ---"<<" WMore: "<<WMore<<" ---"<<"para: "<<W_PQ<<endl;
        MS.push_back(MS_seqSet+WMore+W_PQ); // communication added
    }

    return *min_element(MS.begin(),MS.end());
}

vector<Cnode*> splitSubtreesVector(Ctree *tree, unsigned int num_p, double beta, double & ms_seq, double & ms_all){
    vector<Cnode*> PQ(1,tree->GetRoot());
    vector<double> MS(1,tree->GetRoot()->GetMSCost(beta));
    vector<Cnode*> SeqSet;
    Cnode* node;
    
    if (beta<=0) {
        return SeqSet; //beta should be positive
    }
    
    double MS_seqSet=0;
    double WMore,W_PQ;
    while (!PQ.back()->IsLeaf()) {//only when head(PQ) is a leaf node, W equals to w
        node = PQ.back();
        SeqSet.push_back(node);
        MS_seqSet = MS_seqSet+node->GetMSW();
        node->SetCut();
        PQ.pop_back();
        for (vector<Cnode*>::iterator iter=node->GetChildren()->begin();iter!=node->GetChildren()->end();iter++) {
            (*iter)->SetCut();
            PQ.push_back(*iter);
        }
        
        WMore=0;
        if (PQ.size()>=num_p) {
            sort(PQ.begin(), PQ.end(), cmp_more);
            for (unsigned int i=0;i<(PQ.size()-num_p+1);++i) {
                WMore+=PQ[i]->GetMSCost(beta);// no communication
            }
        }
        
        sort(PQ.begin(), PQ.end(),cmp_PQ);//non-decreasing sort, communication counted
        W_PQ=PQ.back()->GetMSCCost();
        
        //cout<<"seqSet: "<<MS_seqSet<<" ---"<<" WMore: "<<WMore<<" ---"<<"para: "<<MS_para.back()->GetMSCCost()<<endl;
        MS.push_back(MS_seqSet+WMore+W_PQ); // communication added
    }
    
    //return edges we cut instead of makespan
    unsigned long back_step = MS.end()-min_element(MS.begin(),MS.end())-1;
    vector<Cnode*>::iterator it;
    for (unsigned int i=0; i<back_step; ++i) {
        node = SeqSet.back();
        SeqSet.pop_back();
        for (vector<Cnode*>::iterator child=node->GetChildren()->begin(); child!=node->GetChildren()->end(); ++child) {
            it = find(PQ.begin(), PQ.end(),*child);
            (*it)->SetCut();
            PQ.erase(it);
        }
        node->SetCut();
        PQ.push_back(node);
    }
    
    ms_all=*min_element(MS.begin(),MS.end());
    ms_seq=ms_all-PQ.back()->GetMSCost();
    
    if (num_p<=PQ.size()) {
        vector<Cnode*> PQ_return;
        it=PQ.begin();
        PQ_return.assign(PQ.end()-num_p+1,PQ.end());
        for (vector<Cnode*>::iterator ite=PQ.begin(); ite!=PQ.end()-num_p+1; ++ite) {
            (*ite)->SetCut();//set it as non-cut
        }
        return PQ_return;
    }else{
        return PQ;
    }
}

vector<Cnode*> splitSubtreesVector(Ctree *tree, unsigned int num_p, double beta, double & ms_seq, double & ms_all, Qtree* qtree){
    vector<Cnode*> PQ(1,tree->GetRoot());
    qtree->AddNodes(tree->GetRoot(), 1, beta);
    vector<double> MS(1,qtree->GetMakeSpan());
    vector<Cnode*> SeqSet;
    Cnode* node;
    
    if (beta<=0) {
        return SeqSet; //beta should be positive
    }
    
    if (!PQ.back()->IsLeaf()) {
        node = PQ.back();
        SeqSet.push_back(node);
        PQ.pop_back();
        vector<Cnode*>* children=node->GetChildren();
        for (vector<Cnode*>::iterator ite=children->begin(); ite!=children->end(); ++ite) {
            (*ite)->SetCut();//set it as cut, to initilazie qtree
            PQ.push_back(*ite);
        }

        qtree->Clear();
        qtree->AddNodes(tree->GetRoot(), children->size()+1, beta);
        MS.push_back(qtree->GetMakeSpan(num_p));
        
        sort(PQ.begin(), PQ.end(),cmp_PQ);//non-decreasing sort, communication counted
    }
    
    while(!PQ.back()->IsLeaf()){//only when head(PQ) is a leaf node, W equals to w
        qtree->popupPQ(PQ.back(), num_p, beta);
        MS.push_back(qtree->GetMakeSpan(num_p));
        
        node = PQ.back();
        SeqSet.push_back(node);
        PQ.pop_back();
        for (vector<Cnode*>::iterator iter=node->GetChildren()->begin();iter!=node->GetChildren()->end();iter++) {
            (*iter)->SetCut();//set it as cut
            PQ.push_back(*iter);
        }
        
        sort(PQ.begin(), PQ.end(),cmp_PQ);//non-decreasing sort, communication counted
    }
    
    //return edges we cut instead of makespan
    unsigned long back_step = MS.end()-min_element(MS.begin(),MS.end())-1;
    vector<Cnode*>::iterator it;
    for (unsigned int i=0; i<back_step; ++i) {
        node = SeqSet.back();
        SeqSet.pop_back();
        for (vector<Cnode*>::iterator child=node->GetChildren()->begin(); child!=node->GetChildren()->end(); ++child) {
            it = find(PQ.begin(), PQ.end(),*child);
            (*it)->SetCut();//set it as non-cut
            PQ.erase(it);
        }
        node->SetCut();//set it as cut
        PQ.push_back(node);
    }
    
    ms_all=*min_element(MS.begin(),MS.end());
    
    if (num_p<=PQ.size()) {
        vector<Cnode*> PQ_return;
        it=PQ.begin();
        PQ_return.assign(PQ.end()-num_p+1,PQ.end());
        for (vector<Cnode*>::iterator ite=PQ.begin(); ite!=PQ.end()-num_p+1; ++ite) {
            (*ite)->SetCut();//set it as non-cut
        }
        return PQ_return;
    }else{
        return PQ;
    }
}


//comparing function used for sorting subtrees by non-decreasing node weight, communication counted
bool cmp_Qnode(Qnode *a,Qnode *b){
    return (a->GetNW()+a->GetEW())<(b->GetNW()+b->GetEW());
}

bool cmp_nonincreMS(Cnode *a, Cnode *b){//non-increasing without conmmunication
    return a->GetMSCost()>b->GetMSCost();
}

bool cmp_nondecreMS(Cnode *a, Cnode *b){//non-decreasing without conmmunication
    return a->GetMSCost()<b->GetMSCost();
}

//used for sorting subtrees by non-decreasing Makespan
bool cmp_MS(Qnode *a,Qnode *b){
    return (a->GetMS())<(b->GetMS());
}

void CriticalPath(Qtree* qtree, list<Qnode*> & path){
    Qnode* parent_on_cri = qtree->GetNode(1);// root of qtree
    path.push_back(parent_on_cri);
    while (!parent_on_cri->IsLeaf()) {
        sort(parent_on_cri->GetChildren()->begin(), parent_on_cri->GetChildren()->end(), cmp_MS);
        parent_on_cri=parent_on_cri->GetChildren()->back();
        path.push_back(parent_on_cri);
    }
}

double improvedSplit_v1(Ctree* tree,unsigned int num_p, double beta, unsigned int & num_subtrees){
    double ms_seq,msn;
    vector<double> ms;
    list<Qnode*> cpath;
    vector<Cnode*> Cz=splitSubtreesVector(tree, num_p, beta, ms_seq, msn);
    ms.push_back(msn);
    num_subtrees=1;//needs modification
    
    if (Cz.size()<=1) {//communication is too heavy or num_p is too less
        return msn;
    }
    
    unsigned long q=Cz.size()-1;
    unsigned long add_q=1;
    unsigned long pidle=num_p-Cz.size()-1;
    Qtree* qtree=new Qtree(tree->GetRoot(),Cz.size()+1,beta);
    
    sort(Cz.begin(), Cz.end(), cmp_PQ);//sorting by non-decreasing Makespan with communication
    //Ctree* Tlargest = tree->Cut(Cz.back());//leaf subtree on critical path
    Ctree* Tlargest = SubtreeRooted(Cz.back());//leaf subtree on critical path
    cpath.push_back(qtree->GetNodebyCorID(Cz.back()->GetId()));
    
    vector<Cnode*> Ci;
    while (q+pidle>=2) {
        Ci=splitSubtreesVector(Tlargest, 3, beta, ms_seq, msn);
        
        if (Ci.size()<=1) {//commu is too heavy or Tlargest is too small
            delete qtree;
            return *min_element(ms.begin(), ms.end());
        }
        
        qtree->Split(cpath.back(), &Ci, beta);
        
        if (pidle<2) {
            sort(Cz.begin(), Cz.end(), cmp_nonincreMS);//non-increasing without conmmunication
            for (unsigned int removei=0; removei<(2-pidle); ++removei) {
                Cz.back()->SetCut();//set it non-cut
                Cz.pop_back();
            }
            q=q-(2-pidle);
        }
        
        ms.push_back(qtree->GetMakeSpan(q+add_q+1));
        
        cpath.clear();
        CriticalPath(qtree, cpath);
        Cnode* rootTlargest=cpath.back()->GetCorresNode();
        if (pidle+q<2) {
            break;
        }else{
            Tlargest->Clear();//free the pointer to nodes
            delete Tlargest;
        }
        
        Tlargest = SubtreeRooted(rootTlargest);//leaf subtree on critical path
        
        if (cpath.size()==2) {//parent(T') is root of T
            q=q-1;
            ++add_q;
        }
        
        if (pidle>0) {
            if (pidle<2) {
                pidle=0;
            }else{
                pidle=pidle-2;
            }
        }
    }
    Tlargest->Clear();
    delete Tlargest;
    
    delete qtree;
    msn=*min_element(ms.begin(), ms.end());
    
    return msn;
}

double improvedSplit_v2(Ctree* tree,unsigned int num_p, double beta, unsigned int & num_subtrees){
    double ms_seq,msb,msn;
    double ms_largest;
    vector<Cnode*> Ci;
    vector<double> makespanVec;
    unsigned int s=0;
    ms_largest=tree->GetRoot()->GetMSCost(beta);//initialize MS
    
    vector<Cnode*> Cz=splitSubtreesVector(tree, num_p, beta, ms_seq, msb);
    
    if (Cz.size()<=1) {//communication is too heavy or num_p is too less
        return msb;
    }
    
    makespanVec.push_back(msb);
    Qtree* qtree=new Qtree(tree->GetRoot(), Cz.size()+1, beta);
    num_subtrees=Cz.size()+1;
    sort(Cz.begin(), Cz.end(), cmp_PQ);//non-decreasing Makespan with communication
    ms_largest=Cz.back()->GetMSCCost();
    while (!Cz.empty()&&ms_largest==Cz.back()->GetMSCCost()) {
        ++s;
        tree->Cut(Cz.back(),true);
        qtree->MoveNode(Cz.back());
        Cz.pop_back();
    }
    for (vector<Cnode*>::iterator ite=Cz.begin(); ite!=Cz.end(); ++ite) {
        (*ite)->SetCut();//set it as non-cut;
    }
    
    num_p=num_p-s;
    while (num_p>1){
        qtree->Clear();
        Ci=splitSubtreesVector(tree, num_p, beta, ms_seq, msn, qtree);
        qtree->Clear();
        qtree->AddNodes(tree->GetRoot(), Ci.size()+1, beta);
        makespanVec.push_back(msn);
        
        if (Ci.size()<=1) {//communication is too heavy or num_p is too less
            break;
        }
        
        msb=msn;
        num_subtrees=qtree->GetSize();
        s=0;
        sort(Ci.begin(), Ci.end(), cmp_PQ);
        ms_largest=Ci.back()->GetMSCCost();
        while (!Ci.empty()&&ms_largest==Ci.back()->GetMSCCost()) {
            ++s;
            tree->Cut(Ci.back(),true);
            qtree->MoveNode(Ci.back());
            Ci.pop_back();
        }
        num_p=num_p-s;
        
        for (vector<Cnode*>::iterator ite=Ci.begin(); ite!=Ci.end(); ++ite) {
            (*ite)->SetCut();//set it as non-cut;
        }
    }
    
    msn=*min_element(makespanVec.begin(), makespanVec.end());
    
    qtree->deleteCorCnode();
    delete qtree;
    return msn;
}

//comparing function used for sorting by non-decreasing W-f/beta
bool cmp_W_f(Cnode *a, Cnode *b){
    return (2*a->GetMSCost()-a->GetMSCCost())<(2*b->GetMSCost()-b->GetMSCCost());
}

//does not really cut any node, but label it as cut by SecCut()
double ASAP(Ctree* tree, unsigned int num_p, double beta, unsigned int &num_subtrees, bool consider_commu,unsigned int deepth){
    list<Cnode*> PQ;
    Qtree *quotient_tree;
    double makespan;
    Cnode* largest_node;
    vector<double> makespan_vec;
    vector<unsigned int> cutCnodeID;
    list<Cnode*> children_buffer;
    
    makespan_vec.push_back(tree->GetRoot()->GetMSCost(beta));
    
    vector<Cnode*>* children=tree->GetRoot()->GetChildren();
    for (vector<Cnode*>::iterator ite=children->begin(); ite!=children->end(); ++ite) {
        PQ.push_back(*ite);
        children_buffer.push_back(*ite);
    }
    
    unsigned long buffer_size;
    unsigned int i=1;
    while (i<deepth) {
        buffer_size=children_buffer.size();
        while (buffer_size>0) {
            children=children_buffer.front()->GetChildren();
            children_buffer.pop_front();
            for (vector<Cnode*>::iterator child=children->begin(); child!=children->end(); ++child) {
                children_buffer.push_back(*child);
                PQ.push_back(*child);
            }
            --buffer_size;
        }
        ++i;
    }
    
    Cnode* node_remove;
    list<Cnode*>::iterator node_find;
    unsigned int add_child_deepth;
    while (num_p>1) {
        if (consider_commu==false) {
            PQ.sort(cmp_more);//non-decreasing makespan
        }else{
            PQ.sort(cmp_W_f);//non-decreasing makespan minus communication:W-f/beta
            
            /*cout<<"mscost"<<"            "<<"f/beta"<<"            "<<"substraction"<<endl;
            for (list<Cnode*>::iterator iter=PQ.begin(); iter!=PQ.end(); ++iter) {
                cout<<(*iter)->GetMSCost()<<" "<<(*iter)->GetEW()/beta<<" "<<(*iter)->GetMSCost()-(*iter)->GetEW()/beta<<endl;
                cout<<" "<<(*iter)->GetMSCCost()-(*iter)->GetMSCost()<<endl;
            }*/
        }
        
        if (PQ.empty()) {
            break;//when having more processors than nodes
        }
        
        PQ.back()->SetCut();//set it as cut
        cutCnodeID.push_back(PQ.back()->GetId());
        largest_node=PQ.back();
        
        node_remove=PQ.back()->GetParent();
        add_child_deepth=1;
        node_find=find(PQ.begin(), PQ.end(), node_remove);
        
        while (node_find!=PQ.end()) {
            ++add_child_deepth;
            PQ.erase(node_find);
            node_remove=node_remove->GetParent();
            node_find=find(PQ.begin(), PQ.end(), node_remove);
        }
        PQ.pop_back();
        
        children_buffer.clear();
        children_buffer.push_back(largest_node);
        i=0;
        while (i<deepth-add_child_deepth) {
            ++i;
            buffer_size=children_buffer.size();
            while (buffer_size>0) {
                children=children_buffer.front()->GetChildren();
                children_buffer.pop_front();
                for (vector<Cnode*>::iterator child=children->begin(); child!=children->end(); ++child) {
                    children_buffer.push_back(*child);
                }
                --buffer_size;
            }
        }
        
        while (add_child_deepth>0) {
            buffer_size=children_buffer.size();
            while (buffer_size>0) {
                children=children_buffer.front()->GetChildren();
                children_buffer.pop_front();
                for (vector<Cnode*>::iterator child=children->begin(); child!=children->end(); ++child) {
                    children_buffer.push_back(*child);
                    PQ.push_back(*child);
                }
                --buffer_size;
            }
            --add_child_deepth;
        }
        
        --num_p;
    }
    
    quotient_tree=new Qtree(tree->GetRoot(),cutCnodeID.size()+1,beta);
    
    Qnode* rootQtree=quotient_tree->GetNode(1);
    
    //clock_t time;
    //time=clock();
    unordered_map<unsigned int, Qnode*> cid_qnode_MAP;
    vector<Qnode*>* qnodes=quotient_tree->GetNodes();
    for (vector<Qnode*>::iterator qnode=qnodes->begin(); qnode!=qnodes->end(); ++qnode) {
        cid_qnode_MAP.emplace((*qnode)->GetCorresId(),*qnode);
    }
    //time=clock()-time;
    //printf("Insert elements into MAP took me %d clicks. \n",time);
    
    unordered_map<unsigned int, Qnode*>::iterator qnode;
    double sequential,parallel;
    for (vector<unsigned int>::iterator cnodeid=cutCnodeID.begin(); cnodeid!=cutCnodeID.end(); ++cnodeid) {
        //time=clock();
        qnode=cid_qnode_MAP.find(*cnodeid);
        //time=clock()-time;
        //printf("Find a element took me %d clicks. \n",time);
        
        qnode->second->SetSolid();//set is as a subtree in parallel part
        //time=clock();
        rootQtree->GetMSasap(sequential,parallel);//input file of root doesn't count
        //time=clock()-time;
        //printf("Calculate makespan took me %d clicks. \n",time);

        makespan_vec.push_back(sequential+parallel);
    }

    vector<double>::iterator min_element;
    min_element=std::min_element(makespan_vec.begin(), makespan_vec.end());
    makespan=*min_element;
    
    unsigned long index=min_element-makespan_vec.begin();
    if (index==0) {
        num_subtrees=1;
    }else{
        num_subtrees=index+1;
    }
    
    /*vector<unsigned int>::iterator iter=cutCnodeID.begin()+index;
    for (; iter!=cutCnodeID.end(); ++iter) {
       tree->GetNode((*iter))->SetCut();//set it as non-cut
    }*///for avoidChain
    
    vector<unsigned int>::iterator iter=cutCnodeID.begin();
    for (; iter!=cutCnodeID.end(); ++iter) {
     tree->GetNode((*iter))->SetCut();//set it as non-cut
    }//for comparing ASAP and ASAPc
    
    delete quotient_tree;
    
    return makespan;
}

double AvoidChain(Ctree* tree, unsigned int &num_subtrees,double beta, double & after_merged, Qtree* & out_qtree){
    Qtree *qtree = new Qtree(tree->GetRoot(),num_subtrees,beta);
    
    forward_list<Qnode*> que;
    vector<Qnode*> chain;
    vector<Qnode*> *children;
    vector<Cnode*> children_set;
    vector<unsigned int> index;
    unsigned int label=1;
    double makespan;
    
    Qnode *cur_node;
    que.insert_after(que.before_begin(),qtree->GetNode(1));//initialize with root
    bool add_tail=false;
    bool add_chain=false;
    while (!que.empty()) {
        cur_node=que.front();
        que.pop_front();
        add_tail=false;
        
        if (add_chain==true) {
            index.push_back(label);
            add_chain=false;
        }
        
        if (cur_node->GetChildren()->size()==1) {
            add_chain=true;
            chain.push_back(cur_node);
            cur_node=cur_node->GetChildren()->front();
            ++label;
            add_tail=true;
        }
        
        while(cur_node->GetChildren()->size()==1) {
            ++label;
            chain.push_back(cur_node);
            cur_node=cur_node->GetChildren()->front();
        }
        
        if (add_tail==true) {
            chain.push_back(cur_node);
            ++label;
        }
        
        children=cur_node->GetChildren();
        for (vector<Qnode*>::iterator ite=children->begin(); ite!=children->end(); ++ite) {
            que.push_front(*ite);
        }
    }
    
    //after all chains have been found
    if (chain.empty()) {
        makespan=qtree->GetMakeSpan();
        after_merged=makespan;
        out_qtree=qtree;
        //delete qtree;
        
        return makespan;
    }
    
    //cout<<"chains found"<<endl;
    
    Qnode* fat_node=chain.front();
    unsigned int j=1;
    for (unsigned int i=0; i<index.size(); i++) {
        for (; j<index[i]-1; j++) {
            qtree->MergeChain(fat_node, chain[j]);
        }
        fat_node=chain[j];
        j++;
    }
    
    for (; j<chain.size(); j++) {
        qtree->MergeChain(fat_node, chain[j]);
    }
    
    num_subtrees=qtree->GetSize();
    after_merged=qtree->GetMakeSpan();
    
    /*
    unsigned int num_p_idle = chain.size()-index.size()-1;
    
    //after chains merged
    que.insert_after(que.before_begin(),qtree->GetNode(1));//initialize with root
    unsigned long num_children=0;
    Qnode* front;
    vector<Qnode*>* children_vector;
    Cnode* cur_Cnode;
    while (!que.empty()) {
        cur_Cnode=que.front()->GetCorresNode();
        front=que.front();
        que.pop_front();
        num_children=cur_Cnode->GetChildren()->size();
        if (!front->IsLeaf()) {
            children_vector=front->GetChildren();
            for (vector<Qnode*>::iterator ite=children_vector->begin(); ite!=children_vector->end(); ++ite) {
                que.push_front(*ite);
            }
            
            if (num_children>1) {
                for (vector<Cnode*>::iterator ite=cur_Cnode->GetChildren()->begin(); ite!=cur_Cnode->GetChildren()->end(); ++ite) {
                    if (!(*ite)->GetCut()) {//problem here
                        (*ite)->SetLabel(front->GetId());
                        children_set.push_back(*ite);
                        //bug here!
                    }
                }
            }
        }
    }
    
    sort(children_set.begin(), children_set.end(), cmp_PQ);//makespan with commu, non-decreasing way
    
    unsigned int p_remain;
    if (children_set.size()>num_p_idle) {
        p_remain=num_p_idle;
    }else{p_remain=children_set.size();}
    
    for (vector<Cnode*>::iterator ite=children_set.begin(); ite<children_set.begin()+p_remain; ite++) {
        Qnode* child = new Qnode(qtree->GetSize()+1,(*ite)->GetMSCost(),(*ite)->GetEW()/beta);
        qtree->AddNode(child);
        child->SetCorresNode((*ite));
        child->SetParent_id((*ite)->GetLabel());
        child->SetParent(qtree->GetNode((*ite)->GetLabel()));
        qtree->GetNode((*ite)->GetLabel())->AddChild(child);
        qtree->GetNode((*ite)->GetLabel())->SetNW(qtree->GetNode((*ite)->GetLabel())->GetNW()-(*ite)->GetMSCost());
    }

    makespan=qtree->GetMakeSpan();*/
    
    /*vector<Qnode*>::iterator iter=qtree->GetNodes()->begin();
    for (; iter!=qtree->GetNodes()->end(); ++iter) {
        (*iter)->GetCorresNode()->SetCut();//set it as non-cut
    }*/
    
    out_qtree=qtree;
    //delete qtree;
    
    return after_merged;
}

bool up(Ctree* child_tree,Ctree* parent_tree, Qnode* child_Qnode,double memory,double beta){
    Cnode* old_root;
    Cnode* parent;
    bool cut_in_brothers=false;
    vector<Cnode*> children;
    vector<Cnode*> store_children;
    Cnode* back_node;
    list<Cnode*>* l_init = new list<Cnode*>;
    schedule_t* s_init = new schedule_t;
    double cut_val;
    double m_peak;
    int count;
    bool stop=false;
    vector<double> ms_vec;
    
    double ms_0=child_Qnode->GetParent()->GetMS();
    //cout<<"ms_0: "<<ms_0<<endl;
    ms_vec.push_back(ms_0);//ms0
    Cnode* original_root=child_tree->GetRoot();
    unsigned long original_root_poi=child_tree->GetRootId();
    unsigned long original_size=child_tree->GetNodes()->size();
    vector<double>::iterator min_ite;
    unsigned long step=0;
    
    do {
        old_root=child_tree->GetRoot();
        parent = old_root->GetParent();
        //cout<<" now root is "<<old_root->GetId()<<" "<<"it's partent "<<parent->GetId()<<endl;
        
        if (parent->GetParentId()==0) {
            stop = true;// when its parent node is a root
            break;
        }
        
        children.clear();
        store_children.clear();
        
        for (vector<Cnode*>::iterator iter=parent->GetChildren()->begin(); iter!=parent->GetChildren()->end(); ++iter) {
            if ((*iter)->GetId()!=old_root->GetId()) {
                children.push_back(*iter);
                store_children.push_back(*iter);
            }
        }
        
        while (!children.empty()) {
            if (children.back()->GetCut()) {
                cut_in_brothers=true;
                stop = true;
                break;
            }
            back_node=children.back();
            children.pop_back();
            for (vector<Cnode*>::iterator iter=back_node->GetChildren()->begin(); iter!=back_node->GetChildren()->end(); ++iter) {
                children.push_back(*iter);
                store_children.push_back(*iter);
            }
        }
        
        if (cut_in_brothers==false) {
            for (vector<Cnode*>::iterator iter=store_children.begin(); iter!=store_children.end(); ++iter) {
                child_tree->AddNode(*iter);
            }
            child_tree->AddNode(parent);
            //do not change tree id
            child_tree->SetRootId(child_tree->GetNodes()->size());
            
            explore(child_tree->GetRoot(), memory, l_init, s_init, cut_val, *l_init, *s_init, m_peak, 1, 0, count);
            //cout<<"cut_val: "<<m_peak<<endl;
            
            if (l_init->empty()) {
                old_root->SetCut();//restore this edge
                old_root->SetParentId(parent->GetId());
                parent->SetCut();//cut this edge
                parent->SetParentId(0);
                //cout<<" remove node "<<parent->GetId()<<endl;
                parent_tree->RemoveNode(parent);
                unsigned int node_weight_exchan=parent->GetMSW();
                for (vector<Cnode*>::iterator ite=store_children.begin(); ite!=store_children.end(); ++ite) {
                    node_weight_exchan=node_weight_exchan+(*ite)->GetMSW();
                    parent_tree->RemoveNode(*ite);
                }
                
                if (parent_tree->GetRootId()!=1) {
                    parent_tree->SetRootId(parent_tree->GetNodes()->size());//update the RootId when the root is in the end
                }
                
                //update corresponding Qnodes
                child_Qnode->SetCorresId(parent->GetId());
                child_Qnode->SetCorresNode(parent);
                child_Qnode->SetNW(child_Qnode->GetNW()+node_weight_exchan);
                child_Qnode->SetEW(parent->GetEW()/beta);
                child_Qnode->SetMS(0);//set ms equals to 0, it will be recomputed. so as to parent
                
                child_Qnode->GetParent()->SetNW(child_Qnode->GetParent()->GetNW()-node_weight_exchan);
                //cout<<"ms_i: "<<child_Qnode->GetParent()->GetMS()<<endl;
                ms_vec.push_back(child_Qnode->GetParent()->GetMS());//ms_i
            }else{
                for (vector<Cnode*>::iterator iter=store_children.begin(); iter!=store_children.end(); ++iter) {
                    child_tree->RemoveNode(*iter);
                }
                child_tree->RemoveNode(parent);
                if (child_tree->GetNodeByPos(0)->GetParentId()==0) {
                    child_tree->SetRootId(1);
                }else{child_tree->SetRootId(child_tree->GetNodes()->size());}
                stop=true;
            }
        }else{
            stop=true;
        }
    } while (stop==false);
    
    min_ite=min_element(ms_vec.begin(), ms_vec.end());
    step=min_ite-ms_vec.begin();
    
    if (ms_vec.size()==1) {
        return false;
    }
    
    if (step==ms_vec.size()-1) {
        return true;
    }else{// if the min_element is not the last element
        Cnode* new_root=original_root;
        for (unsigned int i=0; i<step; ++i) {
            new_root=new_root->GetParent();
        }
        
        child_tree->GetRoot()->SetCut();//set it as non-cut
        child_tree->GetRoot()->GetParent()->SetLabel(child_tree->GetRoot()->GetParent()->GetId());
        child_tree->GetRoot()->SetParentId(child_tree->GetRoot()->GetParent()->GetId());
        
        new_root->SetCut();//set it as cut
        new_root->GetParent()->SetLabel(-1);
        new_root->SetParentId(0);
        
        child_Qnode->SetCorresId(new_root->GetId());
        child_Qnode->SetCorresNode(new_root);
        child_Qnode->SetEW(new_root->GetEW()/beta);
        
        const vector<Cnode*>* nodes_childtree = child_tree->GetNodes();
        vector<Cnode*>::const_iterator index;
        
        if (step!=0) {
            index=find(nodes_childtree->begin(), nodes_childtree->end(), new_root);
        }else{
            index=nodes_childtree->begin()+original_size-1;
        }
        
        double work_trans = 0;
        ++index;
        vector<Cnode*>::const_iterator index_copy=index;
        for (; index!=nodes_childtree->end(); ++index) {
            work_trans=work_trans+(*index)->GetMSW();
            parent_tree->AddNode(*index);
        }
        child_tree->RemoveNodes(index_copy, nodes_childtree->end());
        
        if (step!=0) {
            child_tree->SetRootId(child_tree->GetNodes()->size());
        }else{
            if (original_root_poi==1) {
                child_tree->SetRootId(1);
            }else{
                child_tree->SetRootId(original_size);
            }
        }
        
        child_Qnode->SetNW(child_Qnode->GetNW()-work_trans);
        child_Qnode->GetParent()->SetNW(child_Qnode->GetParent()->GetNW()+work_trans);
        
        if (step==0) {
            return false;
        }else{
            return true;
        }
    }
    
    return false;
}

double Upper(Ctree* tree, Qtree* qtree, vector<Ctree*>* subtrees, unsigned int num_subtrees, double memory, double beta){
    double msn=qtree->GetMakeSpan(qtree->GetSize());
    double gain=1,msb;
    vector<Qnode*> par_que;
    list<Qnode*> chi_que;
    Qnode* Q_T;
    Ctree* T;
    Ctree* parent_T;
    bool up_sucess = false;
    vector<Qnode*>* children;
    double ms_temp;
    
    unsigned int index=1;
    while (gain>0) {
        //cout<<"----gain>0, enter round "<<index<<endl;
        ++index;
        msb=msn;
        par_que.push_back(qtree->GetNode(1));//initialize with root of Q
        while (!par_que.empty()) {
            children=par_que.back()->GetChildren();
            for (vector<Qnode*>::iterator ite=children->begin(); ite!=children->end(); ++ite) {
                ms_temp=(*ite)->GetMS();
                chi_que.push_back(*ite);
            }
            
            par_que.pop_back();
            chi_que.sort(cmp_Qnode);
            while (!chi_que.empty()) {
                Q_T=chi_que.front();
                chi_que.pop_front();
                par_que.push_back(Q_T);
                T=subtrees->at(Q_T->GetId()-1);
                parent_T=subtrees->at(Q_T->GetParent_id()-1);
                //cout<<"-try to up tree "<<T->GetTreeId()<<endl;
                up_sucess=up(T,parent_T,Q_T,memory,beta);
                /*if (up_sucess) {
                    cout<<"up success, tree's root now is "<<T->GetRoot()->GetId()<<endl;
                }*/
            }
        }
        msn=qtree->GetMakeSpan(qtree->GetSize());
        gain=msb-msn;
    }
    
    /*vector<Qnode*>* qnodes= qtree->GetNodes();
    vector<Qnode*>::iterator iter=qnodes->begin();
    ++iter;
    for (; iter!=qnodes->end(); ++iter) {
        (*iter)->GetCorresNode()->SetParentId((*iter)->GetCorresNode()->GetParent()->GetId());
        (*iter)->GetCorresNode()->SetCut();//set it as non-cut
    }*/// larsav will turn nodes to non-cut
    
    return msb;
}

//sort Cnode by id
bool sort_id(Cnode *a, Cnode *b){
    return (a->GetId() < b->GetId());
}

double NLargest(Qtree* qtree, Qnode* QT, bool really_cut, double beta){
    list<Cnode*> children_nodes;
    list<Cnode*> nodes_on_path;
    Cnode* croot=QT->GetCorresNode();
    Cnode* cnode;
    Cnode* cnode_pre;
    double max_para_computation = 0;
    for (vector<Qnode*>::iterator q_child=QT->GetChildren()->begin(); q_child!=QT->GetChildren()->end(); ++q_child) {
        if (((*q_child)->GetNW()+(*q_child)->GetEW())>max_para_computation) {
            max_para_computation=(*q_child)->GetNW()+(*q_child)->GetEW();
        }
        
        cnode_pre=(*q_child)->GetCorresNode();
        cnode=cnode_pre->GetParent();
        nodes_on_path.push_back(cnode);
        while (cnode!=croot) {
            for (vector<Cnode*>::iterator c_child=cnode->GetChildren()->begin(); c_child!=cnode->GetChildren()->end(); ++c_child) {
                if ((*c_child)->GetCut()==false) {
                    children_nodes.push_back(*c_child);
                }
            }
            cnode_pre=cnode;
            cnode=cnode->GetParent();
            nodes_on_path.push_back(cnode);
        }
    }
    
    for (vector<Cnode*>::iterator c_child=croot->GetChildren()->begin(); c_child!=croot->GetChildren()->end(); ++c_child) {
        if ((*c_child)->GetCut()==false) {
            children_nodes.push_back(*c_child);
        }
    }
    
    children_nodes.sort(sort_id);
    children_nodes.unique();
    nodes_on_path.sort(sort_id);
    nodes_on_path.unique();
    list<Cnode*>::iterator iter=set_difference(children_nodes.begin(), children_nodes.end(), nodes_on_path.begin(), nodes_on_path.end(), children_nodes.begin(),sort_id);
    
    unsigned long dist=std::distance(children_nodes.begin(),iter);
    children_nodes.resize(dist);
    
    children_nodes.sort(cmp_nondecreMS);
    
    double time_saving=0;

    if (!children_nodes.empty()) {
        cnode=children_nodes.back();
        if (max_para_computation>(cnode->GetMSCCost())) {
            time_saving=cnode->GetMSCost();
        }else{
            time_saving=max_para_computation-(cnode->GetMSCCost()-cnode->GetMSCost());
        }
    }
    
    if (really_cut&&time_saving>0) {
        list<Cnode*>::reverse_iterator rit=children_nodes.rbegin();
        (*rit)->SetCut();
        Qnode* child = new Qnode(qtree->GetSize()+1,(*rit)->GetMSCost(),(*rit)->GetEW()/beta);
        qtree->AddNode(child);
        child->SetParent_id(QT->GetId());
        child->SetParent(QT);
        child->SetCorresId((*rit)->GetId());
        child->SetCorresNode((*rit));
        QT->AddChild(child);
        QT->SetNW(QT->GetNW()-child->GetNW());
    }
    
    return time_saving;
}

double cut_two_subtrees(Qtree* qtree, Qnode* QT, bool really_cut, double beta){
    double time_save_i=0;
    Cnode* cnode=QT->GetCorresNode();
    
    while (cnode->GetChildren()->size()==1) {
        cnode=cnode->GetChildren()->front();
    }
    
    if (cnode->GetChildren()->size()>1) {
        sort(cnode->GetChildren()->begin(), cnode->GetChildren()->end(),cmp_nondecreMS);
        vector<Cnode*>::reverse_iterator rit=cnode->GetChildren()->rbegin();
        cnode=cnode->GetChildren()->back();
        double ms_node_a=cnode->GetMSCost(beta);
        double msc_node_a=cnode->GetMSCCost();
        advance(rit,1);
        cnode=(*rit);
        time_save_i=ms_node_a+cnode->GetMSCost(beta)-max(msc_node_a, cnode->GetMSCCost());
    }
    
    if (really_cut==true) {
        vector<Cnode*>::reverse_iterator rit=cnode->GetParent()->GetChildren()->rbegin();
        for (unsigned int i=0; i<2; ++i,++rit) {
            (*rit)->SetCut();//set it as cut
            Qnode* child = new Qnode(qtree->GetSize()+1,(*rit)->GetMSCost(),(*rit)->GetEW()/beta);
            qtree->AddNode(child);
            child->SetParent_id(QT->GetId());
            child->SetParent(QT);
            child->SetCorresId((*rit)->GetId());
            child->SetCorresNode(*rit);
            QT->AddChild(child);
            QT->SetNW(QT->GetNW()-child->GetNW());
        }
    }
    return time_save_i;
};

double LarSav(Ctree* tree, Qtree* qtree, unsigned int & num_subtrees, unsigned int num_p, double beta){
    list<Qnode*> CriPat;
    double cost = tree->GetRoot()->GetMSCost(beta);
    
    CriticalPath(qtree,CriPat);
    
    if (num_p<=num_subtrees) {
        double ms_return=qtree->GetMakeSpan(qtree->GetSize());
        vector<Qnode*>* qnodes= qtree->GetNodes();
        vector<Qnode*>::iterator iter=qnodes->begin();
        ++iter;
        for (; iter!=qnodes->end(); ++iter) {
            (*iter)->GetCorresNode()->SetCut();//set it as non-cut
            (*iter)->GetCorresNode()->SetParentId((*iter)->GetCorresNode()->GetParent()->GetId());//recover its parent id
        }
        delete qtree;
        return ms_return;
    }
    
    long n=num_p-num_subtrees;
    vector<double> time_saving;
    double time_save_i;
    list<Qnode*>::reverse_iterator rit;
    double ms_before=qtree->GetMakeSpan(qtree->GetSize());
    double ms_now;
    while (n>0) {
        time_saving.clear();
        rit=CriPat.rbegin();
        time_save_i=cut_two_subtrees(qtree,*rit,false,beta);
        time_saving.push_back(time_save_i);
        
        advance(rit, 1);
        for (; rit!=CriPat.rend(); advance(rit, 1)) {
            time_save_i=NLargest(qtree,*rit,false,beta);
            time_saving.push_back(time_save_i);
        }
        
        vector<double>::iterator max_ele = max_element(time_saving.begin(), time_saving.end());
        if ((*max_ele)<=0) {// if we can not reduce the makespan, break the while loop
            break;
        }
        
        rit=CriPat.rbegin();
        advance(rit,max_ele-time_saving.begin());
        if (rit==CriPat.rbegin()&&n>=2) {
            cut_two_subtrees(qtree, *rit, true, beta);//cut it really
            n=n-2;
            num_subtrees=num_subtrees+2;
        }else{
            NLargest(qtree,*rit, true,beta);//cut it really
            n=n-1;
            ++num_subtrees;
        }
        
        ms_now=qtree->GetMakeSpan(qtree->GetSize());
        if (ms_now>=ms_before) {
            break;
        }else{
            ms_before=ms_now;
        }
        
        CriPat.clear();
        CriticalPath(qtree, CriPat);
    }
    
    vector<Qnode*>* qnodes= qtree->GetNodes();
    vector<Qnode*>::iterator iter=qnodes->begin();
    ++iter;
    for (; iter!=qnodes->end(); ++iter) {
        (*iter)->GetCorresNode()->SetCut();//set it as non-cut
        (*iter)->GetCorresNode()->SetParentId((*iter)->GetCorresNode()->GetParent()->GetId());//recover its parent id
    }
    
    delete qtree;
    return ms_before;
}

void Immediately(Ctree* tree, unsigned long N, double * nwghts, double * ewghts, int * chstart, int * children, int * schedule, double m_availble, unsigned int & num_para_subtrees){
    double memory_occupied = ewghts[schedule[N-1]];
    
    //cout<<"current task:";
    for (unsigned long rank = N-1; rank>=1; rank--) {
        int cur_task_id = schedule[rank];
        //cout<<" "<<cur_task_id;
        if (cur_task_id !=0) { //=0 means this node has already been moved to another processor
            double node_cost = ewghts[cur_task_id] + nwghts[cur_task_id];
            for(int j = chstart[cur_task_id];j<chstart[cur_task_id+1];j++){
                node_cost += ewghts[children[j]];
            }
            
            double data_to_unload = memory_occupied + node_cost - ewghts[cur_task_id] - m_availble;
            if (data_to_unload>0) {// schedule the subtree that is rooted at this node onto another processor
                //cout<<cur_task_id<<"(move)"<<endl;
                tree->GetNode(cur_task_id)->SetCut();// set it cut, used for building a quotient tree later
                ++num_para_subtrees;
                vector<unsigned int> allNodes;
                list<unsigned int> queue;
                queue.push_back(cur_task_id);
                
                do {
                    unsigned int child_start=*(chstart+queue.front());
                    unsigned int child_end=*(chstart+queue.front()+1);
                    allNodes.push_back(queue.front());
                    queue.pop_front();
                    for (unsigned int i=child_start; i<child_end; ++i) {
                        queue.push_back(*(children+i));
                    }
                } while (!queue.empty());
                
                unsigned long subtree_size=allNodes.size();
                unsigned long subtree_size_backup = subtree_size;
                int* schedule_sub = new int[subtree_size+1];
                schedule_sub[subtree_size]=cur_task_id;
                for (long i=rank-1; i>=0;i-- ) {
                    vector<unsigned int>::iterator iter=find(allNodes.begin(), allNodes.end(), schedule[i]);
                    if (iter!=allNodes.end()) {
                        --subtree_size;
                        schedule_sub[subtree_size]=schedule[i];
                        schedule[i]=0;//IO counter will pass 0;
                        allNodes.erase(iter);
                    }
                    if (allNodes.size()==1) {
                        break;
                    }
                }
                
                Immediately(tree, subtree_size_backup+1, nwghts, ewghts, chstart, children, schedule_sub, m_availble,num_para_subtrees);
                
                memory_occupied -= ewghts[cur_task_id];
                memory_occupied = max(0.0,memory_occupied);
            }else{//memory is enough for executing this node
                memory_occupied += node_cost - 2*ewghts[cur_task_id] - nwghts[cur_task_id];
                memory_occupied = max(0.0,memory_occupied);
            }
    }
}
    //cout<<endl;
}

















