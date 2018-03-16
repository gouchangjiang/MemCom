#!/usr/bin/python

import random
import math
import os

def rand_nb_of_children():
    nb_of_chilren_int = [1,2,21,22] # change this to have a different distribution of degree
    return random.choice(nb_of_chilren_int)

def rand_edge_size():
    return int(round(min(max(random.expovariate(1.0)*100,10),2000))) # edge size is between 10 to 2,000

def generate_tree(start_index, stop_index, parent_index, output_buffer, dot_buffer=None,rec=0):
    # Generates a tree whose node indices go from start_index to stop_index-1, and whose root node as parent_index as parent      
    outgoing_edge_size = rand_edge_size()
    processing_size = outgoing_edge_size * 3
    processing_time = rand_edge_size() 
    print >> output_buffer, start_index, parent_index, processing_size, processing_time, outgoing_edge_size
    if dot_buffer is not None:
        print >> dot_buffer, start_index, "[shape=circle, label=",start_index,", width=", math.sqrt(processing_size), ", fixedsize=true, style=filled, fillcolor=\"#8EC13A\"]"
        if parent_index>0:
            print >> dot_buffer, parent_index, "->", start_index, "[dir=back]"

    if stop_index>start_index+1:
        nb_of_chilren = min(rand_nb_of_children(), stop_index - start_index-1)
        children_id_list = sorted(random.sample(xrange(start_index+2, stop_index), nb_of_chilren-1))
        # child 0: range [start_index+1, children_id_list[0]-1[
        # child 1: range [children_id_list[0], children_id_list[1]-1[
        # child i: range [children_id_list[i-1], children_id_list[i]-1[
        # child n-1: range [children_id_list[n-2], stop_index[

        children_id_list.append(stop_index)
        for i in xrange(nb_of_chilren):
            if i==0:
                start = start_index+1
            else:
                start = children_id_list[i-1]
            stop = children_id_list[i]
            generate_tree(start, stop, start_index, output_buffer, dot_buffer,rec+1)
            


output_dir="./trees_random_10/" # change this
for nb_of_nodes in [1000, 2000 ,3000 ,4000 ,5000, 6000]:
    for iteration in xrange(500):
        tree_buffer = open(output_dir+"tree_%07d_%02d.txt" % (nb_of_nodes, iteration),"w")
        generate_tree(1, nb_of_nodes+1, 0, tree_buffer)
        tree_buffer.close()

        
        
