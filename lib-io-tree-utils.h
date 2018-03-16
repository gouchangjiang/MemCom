 /*
 *  lib-io-tree-utils.h
 *  lib-io-tree
 *
 *  Created by defbond on 8/10/10.
 *  Copyright 2010 LIP/ENS-Lyon. All rights reserved.
 *
 */
#ifndef LIB_IO_TREE_UTILS_H
#define LIB_IO_TREE_UTILS_H

#ifdef __cplusplus

#include <ostream>
#include <iostream>
#include <stdio.h>
#include <list>
#include <algorithm>

#include <vector>
#include <assert.h>
#include <map>

#ifndef DEBUG_MEMUSAGE
#define DEBUG_MEMUSAGE 0 
#endif
#ifndef VERBOSE
#define VERBOSE 0
#endif
#ifndef STRONG_ASSERT
#define STRONG_ASSERT 0
#endif


using namespace std;

#ifndef MAX_COMBI_SIZE
#define MAX_COMBI_SIZE 5
#endif

typedef enum {FURTHEST_NODE=1, BEST_K_COMBI, BEST_FIT_ABS, FIRST_FIT_ABS, BEST_FIT, FIRST_FIT,BEST_INC_COMBI, BEST_COMBI,LARGEST_FIT} io_method_t;


double u_wseconds(void);

class Cnode{
  protected:
    bool cost_computed;
    double cost;
    double MS_cost;
    double MSC_cost; //add communication
    double edge_weight;
    double node_weight;
    double MS_weight;//assume execution time for any node is larger than 0
    vector<Cnode*> * children;
    Cnode * parent;
    unsigned int parent_id;
    unsigned int id;
    bool cut; // label if this edge is cut
    int label;

  public :
    Cnode(){
      id = 0;
      Mpeak = 0;
      parent_id = 0;
      Mavail = 0;
      parent = 0;
      cost_computed = false;
        MS_cost=0;
      children = new vector<Cnode*>();
        cut=false;
    }
    Cnode(double nw, double ew, double mw){
      id = 0;
      Mpeak = 0;
      parent_id = 0;
      Mavail = 0;
      parent=0;
      cost_computed = false;
        MS_cost=0;
        cut=false;
      children = new vector<Cnode*>(); 

      edge_weight =ew;
      node_weight = nw;
      MS_weight = mw;
    }
    Cnode(unsigned int pparent_id,double nw, double ew, double mw){
      id = 0;
      Mpeak = 0;
      Mavail = 0;
      parent=0;
      cost_computed = false;
        MS_cost=0;
        cut=false;
      children = new vector<Cnode*>(); 

      edge_weight =ew;
      node_weight = nw;
        MS_weight = mw;
      parent_id = pparent_id;
    }

    ~Cnode(){

      for(vector<Cnode*>::iterator iter = children->begin();iter!=children->end();iter++){
          //cout<<"delete Cnode "<<(*iter)->GetId()<<endl;
        delete *iter;
      }	
      delete children;
    }


    void SetParent(Cnode * pparent){
      this->parent =pparent;
    }


    void AddChild(Cnode * pchild){
      this->children->push_back(pchild);
      cost_computed= false;
    }

    vector<Cnode*> * GetChildren(){
      return children;
    }

    bool OnlyOneChild(){
        if (children->size()==1) {
            return true;
        } else {
            return false;
        }
    }
    
    Cnode * GetChild(unsigned int node_id){
      return children->at(node_id);
    }
    
    void RemoveChild(Cnode *pchild){
        vector<Cnode*>:: const_iterator ite= find(children->begin(), children->end(), pchild);
        children->erase(ite);
        this->cost_computed=false;
    }

    Cnode * GetParent(){
      return parent;
    }


    bool IsLeaf() const{
      return children->size()==0;
    }

    bool IsRoot() const{
      return parent_id==0;
    }

    double GetCost(){
      if(!cost_computed){
        cost = edge_weight + node_weight;
        for(vector<Cnode*>::iterator iter=children->begin();iter!=children->end();iter++){
          cost += (*iter)->GetEW();
        }
        cost_computed = true;
      }
        return cost;
    }

    double GetMSCost(double beta) {
        //if(MS_cost<=0){
            MS_cost = MS_weight;
            for(vector<Cnode*>::iterator iter=children->begin();iter!=children->end();iter++){
                MS_cost += (*iter)->GetMSCost(beta);
            }
        //}
        MSC_cost= MS_cost+edge_weight/beta;
        
        return MS_cost;
    }
    
    double GetMSCost() const{
        return MS_cost;
    }
    
    double GetMSCCost() {//add communicaiton cost
        return MSC_cost;
    }
    
    void ResetMS_MSCcost(double reduce){
        MS_cost=MS_cost-reduce;
        MSC_cost=MSC_cost-reduce;
    }
    
    void SetParentId(unsigned int pparent_id){
      parent_id = pparent_id;
    }

    void SetId(unsigned int pid){
      id = pid;
    }

    void SetLabel(int pid){
      label = pid;
    }

    void SetEW(double ew ){
      edge_weight = ew;
    }

    void SetNW(double nw ){
      node_weight = nw;
    }
    
    void SetMSW(double mw){
        MS_weight = mw;
    }

    unsigned int GetParentId() const{
      return parent_id;
    }

    double GetEW() const{
      return edge_weight;
    }

    double GetNW() const{
      return node_weight;
    }

    double GetMSW() const{
        return MS_weight;
    }
    
    void SetCut() {
        cut=!cut;
    }
    
    bool GetCut() const{
        return cut;
    }
    
    unsigned int GetId() const{
      return id;
    }

    int GetLabel() const{
      return label;
    }

    void Print(ostream & out) const{
      out<<max((unsigned int)0,GetParentId())<<" "<<GetNW()<<" "<<GetEW()<<endl;
      for(vector<Cnode*>::iterator iter=children->begin();iter!=children->end();iter++){
        (*iter)->Print(out);
      }

    }
    
    double GetEO() const{
        return cost-edge_weight;
    }
    
    unsigned int Ci;
    double Mpeak;
    double Mavail;
};

class Ctree{
  protected:
    vector<Cnode*> * nodes;
    unsigned int root_index;
    unsigned int root_count;
    unsigned int offset_id;
    unsigned int tree_id;
    
  public:

    Ctree(){
      root_index=0;
      root_count=0;
      offset_id = 0;
        tree_id = 1;
      nodes = new vector<Cnode*>();
    }

    Ctree(int N, int *prnts, double *nwghts, double *ewghts,double *mswghts){
      root_index=1;
      root_count=0;
      offset_id = 0;
      tree_id = 1;
      nodes = new vector<Cnode*>();

      this->AllocateNodes(N);

      for(int i = 1; i < N+1 ; i++){
        Cnode * cur_node = this->GetNode(i);
        cur_node->GetChildren()->clear();
        cur_node->SetEW(ewghts[i]);
        cur_node->SetNW(nwghts[i]);
        cur_node->SetMSW(mswghts[i]);
        cur_node->SetId(i);
        cur_node->SetLabel(i);
      }

      for(int i = 1; i <N+1 ; i++){
        Cnode * cur_node = this->GetNode(i);

        if(prnts[i]<N+1){//prnts[i] !=0
          cur_node->SetParentId(prnts[i]);
          cur_node->SetParent(this->GetNode(prnts[i]));
          this->GetNode(prnts[i])->AddChild(cur_node);
        }
        else{
          cur_node->SetParentId(0);
          this->SetRootId(i);
          this->SetTreeId(i);
        }
      }
        
    }


    ~Ctree(){
      if(root_index!=0 && nodes->size()>0){
        delete GetRoot();
      }
      delete nodes;
    }

    //TODO recoder ca en recursif
    void Print(ostream & out) const{
      //		if(root_index!=0 && nodes->size()>0){
      //			out<<nodes->size()<<endl;
      //			
      //			GetRoot()->Print(out);
      //		}
      //		

      out<<nodes->size()<<endl;

      for(vector<Cnode*>::iterator iter = nodes->begin();iter!=nodes->end();iter++){
        out<<max((unsigned int)0,(*iter)->GetParentId()/*+1-offset_id*/)<<" "<<(*iter)->GetNW()<<" "<<(*iter)->GetEW()<<endl;
      }

    }


    void AllocateNodes(int new_node_count){
      if(root_count>0 && nodes->size()>0){
        delete GetRoot();
      }

      nodes->resize(new_node_count);

      unsigned int i = 0;
      for(vector<Cnode*>::iterator iter = nodes->begin();iter!=nodes->end();iter++){
        *iter = new Cnode();
        (*iter)->SetId(i++);
      }	


      offset_id = nodes->front()->GetId();
    }

    void AddNode(Cnode * newNode){
      nodes->push_back(newNode);
    }
    
    void RemoveNode(Cnode *node){
        vector<Cnode*>::const_iterator iter= find(nodes->begin(), nodes->end(), node);
        if (iter==nodes->end()) {
            cout<<"Does not find the node on this tree!"<<endl;
        }
        nodes->erase(iter);
    }
    
    void RemoveNodes(vector<Cnode*>::const_iterator first, vector<Cnode*>::const_iterator last){
        nodes->erase(first, last);
    }
    
    void Clear(){
        nodes->clear();
    }
    
    void AddRoot(Cnode * newNode){
      root_count++;
      assert(root_count == 1);
      nodes->push_back(newNode);
      root_index = nodes->size()-1;
    }

    Cnode * GetRoot() const{
      return nodes->at(root_index-1);
    }

    unsigned int GetRootId() const{
      return root_index;
    }

    void SetRootId(unsigned int root_id){
      root_index = root_id;
    }

    Cnode * GetNode(unsigned int node_id) const{
      return nodes->at(node_id-1);
    }

    Cnode * GetNodeByPos(unsigned int node_idx) const{
      return nodes->at(node_idx);
    }

    const vector<Cnode*> * GetNodes() const{
      return nodes;
    }
    
    unsigned int GetTreeId(){
        return tree_id;
    }
    
    void SetTreeId(unsigned int _id){
        tree_id = _id;
    }
    
    Ctree * Cut(Cnode* node){
        Cnode *parent = node->GetParent();
        parent->RemoveChild(node);
        this->RemoveNode(node);
        node->SetParentId(0);
        
        double reduce=node->GetMSCost();// reduce the makespan
        while (parent!=this->GetRoot()) {
            parent->ResetMS_MSCcost(reduce);
            parent=parent->GetParent();
        }
        this->GetRoot()->ResetMS_MSCcost(reduce);

        Ctree *subtree = new Ctree();
        
        subtree->SetRootId(1);
        subtree->SetTreeId(node->GetId());
        subtree->AddNode(node);
        
        vector<Cnode*> visit_next;
        vector<Cnode*>::iterator first_node;
        Cnode* end_node;
        if (node->IsLeaf()) {
            //subtree->AddNode(NULL);
            return subtree;
        } else {
            visit_next=*(node->GetChildren());
            while (!visit_next.empty()) {
                subtree->AddNode(visit_next.back());this->RemoveNode(visit_next.back());
                end_node=visit_next.back();
                visit_next.pop_back();
                if (!end_node->IsLeaf()) {
                    visit_next.insert(visit_next.end(), end_node->GetChildren()->begin(), end_node->GetChildren()->end());
                }
            }
            //subtree->AddNode(NULL);
            return subtree;
        }
    }
    
    void Cut(Cnode* node, bool no_return){
        Cnode *parent = node->GetParent();
        parent->RemoveChild(node);
        this->RemoveNode(node);
        
        double reduce=node->GetMSCost();// reduce the makespan
        while (parent!=this->GetRoot()) {
            parent->ResetMS_MSCcost(reduce);
            parent=parent->GetParent();
        }
        this->GetRoot()->ResetMS_MSCcost(reduce);
        
        vector<Cnode*> visit_next;
        vector<Cnode*>::iterator first_node;
        Cnode* end_node;
        if (!node->IsLeaf()) {
            visit_next=*(node->GetChildren());
            while (!visit_next.empty()) {
                this->RemoveNode(visit_next.back());
                end_node=visit_next.back();
                visit_next.pop_back();
                if (!end_node->IsLeaf()) {
                    visit_next.insert(visit_next.end(), end_node->GetChildren()->begin(), end_node->GetChildren()->end());
                }
            }
        }
    }
    
    void Merge(Ctree* child_tree, Cnode* restored_edge){
        nodes->insert(nodes->end(),child_tree->GetNodes()->begin(), child_tree->GetNodes()->end());
        
        restored_edge->AddChild(child_tree->GetRoot());
        restored_edge->ResetMS_MSCcost(-child_tree->GetRoot()->GetMSCost());
        child_tree->GetRoot()->SetParentId(restored_edge->GetId());
    }
};

Ctree* SubtreeRooted(Cnode* node);

class Qnode{
private:
    unsigned int id;
    unsigned int corres_id;
    vector<Qnode*> * children;
    double node_weight;
    double edge_weight;
    double makespan;
    Qnode* parent;
    Cnode* corres_node;
    unsigned int parent_id;
    bool need_change_parent=false;
    bool is_solid=false;
    
public:
    Qnode(){
        id=0;
        corres_id=0;
        node_weight=0;
        edge_weight=0;
        makespan=0;
        children = new vector<Qnode*>();
    }
    
    Qnode(unsigned int id_number, double nw, double ew){
        id = id_number;
        node_weight=nw;
        edge_weight=ew;
        makespan=0;
        children = new vector<Qnode*>();
    }
    
    ~Qnode(){
        for(vector<Qnode*>::iterator iter = children->begin();iter!=children->end();iter++){
            delete *iter;
        }
        delete children;
    }

    void SetId(unsigned int id_value){
        id=id_value;
    }
    
    unsigned int GetId() const{
        return id;
    }
    
    void SetCorresId(unsigned int value){
        corres_id=value;
    }
    
    void SetCorresNode(Cnode* cnode){
        corres_node=cnode;
    }
    
    Cnode* GetCorresNode() const{
        return corres_node;
    }
    
    unsigned int GetCorresId() const{
        return corres_id;
    }
    
    double GetNW() const{
        return node_weight;
    }
    
    double GetEW() const{
        return edge_weight;
    }
    
    void SetEW(double ew){
        edge_weight=ew;
    }
    
    void SetNW(double nw){
        node_weight=nw;
    }
    
    void SetMS(double ms){
        makespan = ms;
    }
    
    double GetMS(){
        //if (makespan==0) {
            double largest=0;
            for (vector<Qnode*>::iterator ite=children->begin(); ite!=children->end(); ++ite) {
                if((*ite)->GetMS()>largest){largest=(*ite)->GetMS();}
            }
            makespan=node_weight+edge_weight+largest;
        //}
        return makespan;
    }
    
    void GetMSasap(double & sequential, double & parallel){
        sequential=node_weight;
        parallel=0;
        double child_seq=0,child_para=0;
        for (vector<Qnode*>::iterator ite=children->begin(); ite!=children->end(); ++ite) {
            if (!(*ite)->IsSolid()) {
                (*ite)->GetMSasap(child_seq,child_para);
                sequential=sequential+child_seq;
                if (child_para>parallel) {
                    parallel=child_para;
                }
            }else{
                (*ite)->GetMSasap(child_seq,child_para);
                child_para=child_seq+(*ite)->GetEW()+child_para;//take communication into account
                if(child_para>parallel){parallel=child_para;}
            }
        }
    }
    
    void AddChild(Qnode* pnode){
        this->children->push_back(pnode);
    }
    
    void RemoveChild(Qnode *pchild){
        vector<Qnode*>:: const_iterator ite= find(children->begin(), children->end(), pchild);
        children->erase(ite);
    }
    
    void SetParent(Qnode* pparent){
        this->parent=pparent;
    }
    
    void SetParent_id(unsigned int pid){
        this->parent_id=pid;
    }
    
    Qnode* GetParent() const{
        return this->parent;
    }
    
    unsigned int GetParent_id() const{
        return this->parent_id;
    }
    
    vector<Qnode*>* GetChildren() const{
        return children;
    }
    
    bool IsLeaf() const{
        if (this->children->empty()) {
            return true;
        }
        return false;
    }
    
    bool NeedNewParent() const{
        return need_change_parent;
    }
    
    void SetChangeParent(){
        need_change_parent = !need_change_parent;
    }
    
    bool IsSolid(){
        return is_solid;
    }
    
    void SetSolid(){
        is_solid=!is_solid;
    }
};

class Qtree{
private:
    vector<Qnode*> * nodes;
    vector<Qnode*> * solid_nodes;

public:
    Qtree(Cnode* root,unsigned long number_nodes,double beta){
        vector<Cnode*> add_weight_queue;
        list<Qnode*> Qnode_queue;
        Cnode* cnode;
        Qnode* cur_node;
        Qnode* child_node;
        double nw=0;
        
        nodes = new vector<Qnode*>();
        this->AllocateNodes(number_nodes);
        
        cur_node=this->GetNode(1);
        cur_node->SetId(1);
        cur_node->SetCorresId(1);
        cur_node->SetEW(0);
        cur_node->SetParent_id(0);
        cur_node->SetCorresNode(root);
        
        Qnode_queue.push_back(cur_node);
        unsigned int id_node=1;
        while (!Qnode_queue.empty()) {
            cnode= Qnode_queue.front()->GetCorresNode();
            nw=cnode->GetMSW();
            
            for (vector<Cnode*>::iterator iter=cnode->GetChildren()->begin(); iter!=cnode->GetChildren()->end(); ++iter) {
                add_weight_queue.push_back(*(iter));
            }
            
            while(!add_weight_queue.empty()){
                cnode = add_weight_queue.back();
                add_weight_queue.pop_back();
                if (cnode->GetCut()==true) {
                    id_node++;
                    cnode->SetLabel(id_node);//label stores corresponding id of Qnode
                    child_node=this->GetNode(id_node);
                    child_node->SetId(id_node);
                    child_node->SetCorresId(cnode->GetId());
                    child_node->SetCorresNode(cnode);
                    Qnode_queue.front()->AddChild(child_node);
                    child_node->SetParent(Qnode_queue.front());
                    child_node->SetParent_id(Qnode_queue.front()->GetId());
                    child_node->SetEW(cnode->GetEW()/beta);
                    Qnode_queue.push_back(child_node);
                    }else{
                    nw = nw + cnode->GetMSW();
                    for (vector<Cnode*>::iterator child=cnode->GetChildren()->begin();child!= cnode->GetChildren()->end(); ++child) {
                        add_weight_queue.push_back(*child);}
                    }
            }
            Qnode_queue.front()->SetNW(nw);
            Qnode_queue.pop_front();
        }
        solid_nodes = new vector<Qnode*>();
    }
    
    Qtree(Cnode* root,unsigned long number_nodes,double beta, vector<Ctree*>& subtrees){
        vector<Cnode*> add_weight_queue;
        list<Qnode*> Qnode_queue;
        Cnode* node;
        Qnode* cur_node;
        Qnode* child_node;
        double nw=0;
        
        nodes = new vector<Qnode*>();
        this->AllocateNodes(number_nodes);
        
        cur_node=this->GetNode(1);//root of Qtree
        cur_node->SetId(1);
        cur_node->SetCorresId(1);
        cur_node->SetEW(0);
        cur_node->SetParent_id(0);
        cur_node->SetCorresNode(root);
        
        //produce the corresponding subtree
        Ctree* subtree = new Ctree();
        subtree->SetRootId(1);
        subtree->SetTreeId(1);
        subtree->AddNode(root);
        subtrees.push_back(subtree);
        
        Qnode_queue.push_back(cur_node);
        unsigned int i=1,subtree_index=0;
        
        while (!Qnode_queue.empty()) {
            node= Qnode_queue.front()->GetCorresNode();
            nw=node->GetMSW();
            
            for (vector<Cnode*>::iterator iter=node->GetChildren()->begin(); iter!=node->GetChildren()->end(); ++iter) {
                add_weight_queue.push_back(*(iter));
            }
            
            while(!add_weight_queue.empty()){
                node = add_weight_queue.back();
                //cout<<"---pop up node "<<node->GetId()<<" from add_weight_queue"<<endl;
                add_weight_queue.pop_back();
                if (node->GetCut()==true){
                    i++;
                    //cout<<"i: "<<i<<" "<<"node "<<node->GetId()<<endl;
                    child_node=this->GetNode(i);
                    child_node->SetId(i);
                    child_node->SetCorresId(node->GetId());
                    child_node->SetCorresNode(node);
                    Qnode_queue.front()->AddChild(child_node);
                    child_node->SetParent(Qnode_queue.front());
                    child_node->SetParent_id(Qnode_queue.front()->GetId());
                    child_node->SetEW(node->GetEW()/beta);
                    Qnode_queue.push_back(child_node);
                    
                    //node->GetParent()->RemoveChild(node); //do not change the tree
                    node->SetParentId(0);
                    
                    //produce the corresponding subtree
                    subtree = new Ctree();
                    subtree->SetRootId(1);
                    subtree->SetTreeId(node->GetId());
                    subtree->AddNode(node);
                    subtrees.push_back(subtree);
                    
                }else{
                    nw = nw + node->GetMSW();
                    subtrees[subtree_index]->AddNode(node);
                    for (vector<Cnode*>::iterator child=node->GetChildren()->begin();child!= node->GetChildren()->end(); ++child) {
                        add_weight_queue.push_back(*child);
                        //cout<<"add node "<<(*child)->GetId()<<" in to add_weight_queue"<<endl;
                    }
                }
            }
            Qnode_queue.front()->SetNW(nw);
            Qnode_queue.pop_front();
            subtree_index++;
        }
        solid_nodes = new vector<Qnode*>();
    }
    
    void AddNodes(Cnode* root, unsigned long number_nodes, double beta){
        vector<Cnode*> add_weight_queue;
        list<Qnode*> Qnode_queue;
        Cnode* cnode;
        Qnode* cur_node;
        Qnode* child_node;
        double nw=0;
        
        Qnode* new_qnode;
        for (unsigned int i=0; i<number_nodes; ++i) {
            new_qnode = new Qnode();
            this->nodes->push_back(new_qnode);
        }
        //this->AllocateNodes(number_nodes);
        
        cur_node=this->GetNode(1);
        cur_node->SetId(1);
        cur_node->SetCorresId(1);
        cur_node->SetEW(0);
        cur_node->SetParent_id(0);
        cur_node->SetCorresNode(root);
        
        Qnode_queue.push_back(cur_node);
        unsigned int id_node=1;
        while (!Qnode_queue.empty()) {
            cnode= Qnode_queue.front()->GetCorresNode();
            nw=cnode->GetMSW();
            
            for (vector<Cnode*>::iterator iter=cnode->GetChildren()->begin(); iter!=cnode->GetChildren()->end(); ++iter) {
                add_weight_queue.push_back(*(iter));
            }
            
            while(!add_weight_queue.empty()){
                cnode = add_weight_queue.back();
                add_weight_queue.pop_back();
                if (cnode->GetCut()==true) {
                    id_node++;
                    cnode->SetLabel(id_node);//label stores corresponding id of Qnode
                    child_node=this->GetNode(id_node);
                    child_node->SetId(id_node);
                    child_node->SetCorresId(cnode->GetId());
                    child_node->SetCorresNode(cnode);
                    Qnode_queue.front()->AddChild(child_node);
                    child_node->SetParent(Qnode_queue.front());
                    child_node->SetParent_id(Qnode_queue.front()->GetId());
                    child_node->SetEW(cnode->GetEW()/beta);
                    Qnode_queue.push_back(child_node);
                }else{
                    nw = nw + cnode->GetMSW();
                    for (vector<Cnode*>::iterator child=cnode->GetChildren()->begin();child!= cnode->GetChildren()->end(); ++child) {
                        add_weight_queue.push_back(*child);}
                }
            }
            Qnode_queue.front()->SetNW(nw);
            Qnode_queue.pop_front();
        }
        
        if (root->GetCut()==false) {
            root->SetCut();//set it cut
        }
        
        if (!solid_nodes->empty()) {
            for (vector<Qnode*>::iterator iter=solid_nodes->begin(); iter!=solid_nodes->end(); ++iter) {
                if ((*iter)->NeedNewParent()) {
                    cnode=(*iter)->GetCorresNode()->GetParent();
                    while (!cnode->GetCut()) {
                        cnode=cnode->GetParent();
                    }
                    (*iter)->SetParent_id(cnode->GetLabel());
                    (*iter)->SetParent(this->GetNode(cnode->GetLabel()));
                    this->GetNode(cnode->GetLabel())->AddChild((*iter));
                }
            }
        }
    }
    
    ~Qtree(){
        if (!nodes->empty()) {
            delete this->GetNode(1);//delete root
        }
        delete nodes;

        delete solid_nodes;
    }
    
    void AllocateNodes(unsigned long new_node_count){
        nodes->resize(new_node_count);
        
        for(vector<Qnode*>::iterator iter = nodes->begin();iter!=nodes->end();iter++){
            *iter = new Qnode();
        }
    }
    
    Qnode * GetNode(unsigned int node_id) const{
        return nodes->at(node_id-1);
    }
    
    vector<Qnode*> * GetNodes() {
        return nodes;
    }
    
    Qnode * GetNodebyCorID(unsigned int corres_id) const{
        vector<Qnode*>::iterator iter=nodes->begin();
        while (iter!=nodes->end()) {
            if ((*iter)->GetCorresId()==corres_id) {
                return *iter;
            }
            ++iter;
        }
        return *iter;
    }
    
    double GetMakeSpan() const{
        return this->GetNode(1)->GetMS();
    }
    
    double GetMakeSpan(unsigned int num_p){
        if (this->nodes->size()>num_p) {
            double parallel_part=0;
            double sequential_part= this->GetNode(1)->GetNW();
            vector<Qnode*>* children=this->GetNode(1)->GetChildren();
            sort(children->begin(), children->end(), cmp_more);
            vector<unsigned int> setnw;
            vector<double> nwbefore;
            for (unsigned int i=0; i<this->nodes->size()-num_p; ++i) {
                if (!children->at(i)->IsSolid()) {
                    sequential_part=sequential_part+children->at(i)->GetNW();
                }
                if (!children->at(i)->IsLeaf()) {
                    nwbefore.push_back(children->at(i)->GetNW());
                    children->at(i)->SetNW(0);
                    setnw.push_back(i);
                }
            }
            
            for (vector<Qnode*>::iterator iter=children->begin(); iter!=children->end(); ++iter) {
                if ((*iter)->GetMS()>parallel_part) {
                    parallel_part=(*iter)->GetMS();
                }
            }
            
            if (!setnw.empty()) {
                children->at(setnw.back())->SetNW(nwbefore.back());
                setnw.pop_back();
                nwbefore.pop_back();
            }
            
            return (sequential_part+parallel_part);
        }else{
            return this->GetNode(1)->GetMS();
        }
    }
    
    void MergeChain(Qnode* parent, Qnode* child){
        double new_node_weight = parent->GetNW()+child->GetNW();
        parent->SetNW(new_node_weight);
        parent->GetChildren()->pop_back();
        child->GetCorresNode()->SetCut();//set it as non-cut
        for (vector<Qnode*>::iterator ite=child->GetChildren()->begin(); ite!=child->GetChildren()->end(); ++ite) {
            (*ite)->SetParent_id(parent->GetId());
            (*ite)->SetParent(parent);
            parent->AddChild(*ite);
        }
        this->RemoveNode(child);
    }
    
    void RemoveNode(Qnode* pnode){
        vector<Qnode*>::iterator it = find(nodes->begin(), nodes->end(), pnode);
        nodes->erase(it);
    }
    
    void AddNode(Qnode* pnode){
        nodes->push_back(pnode);
    }
    
    unsigned long GetSize() const{
            return (nodes->size()+solid_nodes->size());
    }
    
    void MoveNode(Cnode* corCnode){
        for (vector<Qnode*>::iterator iter=nodes->begin(); iter!=nodes->end(); ++iter) {
            if ((*iter)->GetCorresNode()==corCnode) {
                (*iter)->SetChangeParent();//make it need find a parent
                (*iter)->SetSolid();//set it as solid node
                for (vector<Qnode*>::iterator child=(*iter)->GetChildren()->begin(); child!=(*iter)->GetChildren()->end(); ++child) {
                    (*child)->SetChangeParent();//make it's children don't need to change parent
                }
                solid_nodes->push_back(*iter);
                nodes->erase(iter);
                break;
            }
        }
    }
    
    void Clear(){
        for (vector<Qnode*>::iterator iter=nodes->begin(); iter!=nodes->end(); ++iter) {
            (*iter)->GetChildren()->clear();
            delete *iter;
        }
        nodes->clear();
    }
    
    void deleteCorCnode(){
        for (vector<Qnode*>::iterator iter=solid_nodes->begin(); iter!=solid_nodes->end(); ++iter) {
            delete (*iter)->GetCorresNode();
        }
    }
    
    void Split(Qnode* prev,vector<Cnode*>* cut_set,double beta){
        double decrement=0;
        unsigned long size_before = this->solid_nodes->size();
        solid_nodes->resize(size_before+cut_set->size());
        for (vector<Cnode*>::iterator ite=cut_set->begin(); ite!=cut_set->end(); ++ite,++size_before) {
            solid_nodes->at(size_before) = new Qnode(size_before+1,(*ite)->GetMSCost(),(*ite)->GetEW()/beta);
            solid_nodes->at(size_before)->SetParent(prev);
            solid_nodes->at(size_before)->SetParent_id(prev->GetId());
            solid_nodes->at(size_before)->SetCorresId((*ite)->GetId());
            solid_nodes->at(size_before)->SetCorresNode((*ite));
            prev->AddChild(solid_nodes->at(size_before));
            solid_nodes->at(size_before)->SetSolid();//set it as a solid node
            decrement=decrement+(*ite)->GetMSCost();
        }
        prev->SetNW(prev->GetNW()-decrement);
    }
    
    void Split(Qnode* parent,Cnode* child,double beta){
        double decrement=0;
        unsigned long size_before = this->nodes->size();
        nodes->resize(size_before+1);
        
        nodes->at(size_before) = new Qnode(size_before+1,child->GetMSCost(),child->GetEW()/beta);
        nodes->at(size_before)->SetParent(parent);
        nodes->at(size_before)->SetParent_id(parent->GetId());
        nodes->at(size_before)->SetCorresId(child->GetId());
        nodes->at(size_before)->SetCorresNode(child);
        parent->AddChild(nodes->at(size_before));
        nodes->at(size_before)->SetSolid();//set it as a solid node
        decrement=decrement+child->GetMSCost();

        parent->SetNW(parent->GetNW()-decrement);
    }
    
    static bool cmp_more(Qnode *a, Qnode *b){
        return a->GetNW()<b->GetNW();//non-communication and non-decreasing sort
    }
    
    void popupPQ(Cnode* largestNode,unsigned int num_p,double beta){
        Qnode* root = this->GetNode(1);
        root->SetMS(0);//it will re-calculate
        root->SetNW(root->GetNW()+largestNode->GetMSW());
        
        Qnode* headPQ = this->GetNodebyCorID(largestNode->GetId());
        vector<Qnode*>* Qchildren = headPQ->GetChildren();
        root->RemoveChild(headPQ);
        this->RemoveNode(headPQ);
        for (vector<Qnode*>::iterator ite=Qchildren->begin(); ite!=Qchildren->end(); ++ite) {
            (*ite)->SetParent(root);
            root->AddChild(*ite);
        }
        Qchildren->clear();
        delete headPQ;
        
        largestNode->SetCut();//set it as non-cut
        vector<Cnode*>* Cchildren = largestNode->GetChildren();
        unsigned long size_before=nodes->size();
        nodes->resize(size_before+Cchildren->size());
        for (vector<Cnode*>::iterator ite=Cchildren->begin(); ite!=Cchildren->end(); ++ite,++size_before) {
            nodes->at(size_before) = new Qnode(size_before+1,(*ite)->GetMSCost(),(*ite)->GetEW()/beta);
            nodes->at(size_before)->SetParent(root);
            nodes->at(size_before)->SetParent_id(1);
            nodes->at(size_before)->SetCorresId((*ite)->GetId());
            nodes->at(size_before)->SetCorresNode((*ite));
            root->AddChild(nodes->at(size_before));
        }
    }
    
};

struct s_node_t {
  int parent;
  vector<int> children;
  double edge_weight;
  double node_weight;
  double Mpeak;
  double Mavail;
};

struct s_io_t {
  unsigned int node;
  double unloaded_data;
};



typedef list<int> schedule_t;
typedef list<Cnode*> cut_t;

typedef pair<unsigned int, double> io_t;
typedef map<unsigned int, double> io_map;
typedef pair<unsigned int, unsigned int> node_sche;
typedef pair<unsigned int, double> node_ew;

struct OrdoLiu_t ;
struct val_seg_t {
  schedule_t::iterator begin;
  unsigned int begin_index;
  schedule_t::iterator end;
  unsigned int end_index;
  OrdoLiu_t * orig_ordo;
  double value;
};

struct OrdoLiu_t {
  double max_pebble_cost;
  double fi;
  list<val_seg_t> val_seg;
  schedule_t schedule;
}; 


void ConvertToLiu(const Ctree * tree_us, Ctree * tree_liu) ;
void ConvertToLiu(const int * oldprnts,const double * oldnwghts,const double * oldewghts, int N,const int* chstart,const int * children, int ** pprnts, double ** pnwghts, double ** pewghts);
void parse_tree(const char *filename,Ctree * tree);
void parse_tree(const char *filename,int * N ,int **prnts,double **nwghts,double **ewghts,double **mswghts);



extern "C" {
#endif
  void po_construct(const int N, const int * prnts, int **chstart,int **chend,int **children, int * root);
  void poaux(const int * chstart, const int * children, int N, int r, int * por, int * label);
#ifdef __cplusplus
} /* closing brace for extern "C" */



bool check_schedule(int * prnts,int * sched,int N);

double MaxOutDegree(Ctree * tree,int quiet);
double MaxOutDegree(int N, double * nwghts, double * ewghts, int * chstart,int * children);

void NextValley(Cnode * node, double available_memory,  double & cut_value, list<Cnode*> & min_sub_cut, list<unsigned int> & sub_schedule, double & Inc, int quiet, int depth,int & count);
double IOCounter(Ctree & tree, schedule_t & sub_schedule, double available_memory,bool divisible,int quiet);
double IOCounter(Ctree* tree,int N, double * nwghts, double * ewghts, int * chstart,int * children, int * schedule, double available_memory,bool divisible,int quiet, unsigned int & com_freq, io_method_t method=FURTHEST_NODE);



#endif
#endif
