# encoding: utf-8

import math
import sys
import itertools
import time

#######################################################################
## This class represents an unoriented cycle-graph.

## For instance, assume a permutation (5, 2, ...). The configuration
## in the graph would be.

## node(0).ap  = -5      ## node(0).ab  = Null
## node(-5).ab = 5       ## node(-5).ap = 0  
## node(5).ap  = -2      ## node(5).ab  = -5 
## node(-2).ab = 2       ## node(-2).ap = 5
########################################################################

class cycle_configuration_graph() :
    def __init__(self, cycles) :

        ## The number of cycles was never calculated
        self.__num_cycles           = False
        self.__num_odd_cycles       = False
        self.__first_indice_shift   = 0 # Tell us the index of the
                                        # edge that is pointed by
                                        # begin_node. If it is
                                        # negative, than we know that
                                        # a mirror have occurred.

        ## self.n is the number of black edges
        self.n = 0
        for cycle in cycles :
            self.n = self.n + len(cycle)
        n = self.n
          
        # Creating nodes
        node_list = []
        node_list = [cycle_graph_node(i) for i in range(2*n)]
        self.begin_node = node_list[0 ]
        self.end_node   = node_list[-1]

        # Creating ap
        for i in range(0,2*n,2) :
            node_list[i  ].ap  = node_list[i+1]
            node_list[i+1].ap  = node_list[i  ]

        # Creating ab
        for i in range(1,2*n-1,2) :
            node_list[i  ].ab  = node_list[i+1]
            node_list[i+1].ab  = node_list[i  ]

        # Creating ac 
        for cycle in cycles :
            for i in range(len(cycle)) :               
                front, back  = -1,-1 # from left of the black edge.
                j = (i+1)%len(cycle)

                if cycle[i] > 0 :
                    front = 2*( cycle[i]  -1 )
                else :
                    front = 2*(-cycle[i]) -1

                if cycle[j] > 0 :
                    back = 2*( cycle[j]) -1
                else :
                    back = 2*(-cycle[j]  -1 )

                node_list[front].ac = node_list[back]
                node_list[back].ac  = node_list[front]

        self.set_values()
        for i in range(0,2*n,1) :
            positive =  min(abs(node_list[i].value),
                            abs(node_list[i].ac.value))
            
############################################################                
################ Rearrangement Operation  #################
############################################################                

    def k4cut(self, i, j, k, l, sigma) :

        original_blocks = [None, [None, None], [None, None], [None, None], None]
        use_block = [0, i!=j, j!=k, k!=l]

        count     = 0

        node = self.begin_node
        # Find the black edges
        while node :
            count = count + 1

            if count == i :
                original_blocks[1][0] = node.ap
                original_blocks[0] = node
            if count == j :
                original_blocks[2][0] = node.ap
                original_blocks[1][1] = node
            if count == k :
                original_blocks[3][0] = node.ap
                original_blocks[2][1] = node
            if count == l :
                original_blocks[4] = node.ap
                original_blocks[3][1] = node
            node = node.ap.ab

        final_blocks = [original_blocks[0],]

        for aux in sigma :
            if use_block[abs(aux)] :
                if aux < 0 :
                    final_blocks += original_blocks[-aux][::-1]
                else :
                    final_blocks += original_blocks[aux]

        final_blocks += [original_blocks[-1]]

        lista_blocos = [andre.index for andre in final_blocks] 

        # Change the edges
        for aux in range(1,len(final_blocks),2) :
            final_blocks[aux-1].ap = final_blocks[aux]
            final_blocks[aux].ap = final_blocks[aux-1]

        self.__num_cycles           = False
        self.__num_odd_cycles       = False
        self.reset_indices()
        return #unp_i, unp_j, unp_k, unp_l


############################################################                
######################## Bounds  ###########################
############################################################                
    def lower_bound(self, num_cycles, hasShort) :
        if hasShort :
            return (float(self.n - num_cycles)+1)/3
        else:
            return (float(self.n - num_cycles))/3

############################################################                
################### Auxiliary Methods  #####################
############################################################                 

    def get_cycles(self, want_vertices = False) :
        self.clean_visit()        

        node = self.end_node        
        cycles    = []
        vertices  = []

        while node :
            if not node.visit :
                cycle_node  = node
                cycle       = []
                cycle_nodes = []

                while not cycle_node.visit :
                    if cycle_node.index % 2 == 0 :
                        cycle.append( -(cycle_node.index+2)/2 )
                    else :
                        cycle.append( +(cycle_node.index+1)/2 )


                    cycle_node.visit = True
                    cycle_nodes.append(cycle_node)
                    cycle_node = cycle_node.ap
                    cycle_node.visit  = True
                    cycle_nodes.append(cycle_node)
                    cycle_node = cycle_node.ac
                    
                cycles.append(tuple(cycle))
                vertices.append(cycle_nodes)

            node = node.ap.ab
        if want_vertices :
            return tuple(cycles), vertices
        return tuple(cycles)

    ## Here we try to give values to the nodes to check if they are a
    ## permutation.
    def set_values(self) :
        node = self.begin_node        

        node.value = 0
        node = node.ac
        node.value = -1

        while node.ab :
            node = node.ab
            node.value = -(node.ab.value)
            node = node.ac

            if node.ac.value > 0 :
                node.value = -(node.ac.value + 1)
            else :
                node.value = -(node.ac.value - 1)
    
    def reset_indices(self) :
        node  = self.begin_node 
        count = 0

        while node :
            node.index = count
            node       = node.ap
            
            count      = count + 1
            
            node.index = count
            node       = node.ab
            
            count      = count + 1

    def num_cycles(self) :
        if type(self.__num_cycles) == bool :
            self.calculate_cycles()
        return self.__num_cycles

    def num_odd_cycles(self) :
        if type(self.__num_odd_cycles) == bool :
            self.calculate_cycles()
        return self.__num_odd_cycles


    def calculate_cycles(self) :
        cycles, vertices = self.get_cycles(want_vertices = True)
        num_cycles = len(cycles)
        num_odd = 0

        for cycle in cycles :
            if len(cycle) % 2 == 1 :
                num_odd = num_odd + 1

        self.__num_cycles           = num_cycles
        self.__num_odd_cycles       = num_odd

        for i in range(len(cycles)) :            
            cycle = cycles[i]
            size = len(cycle)

            direction = 1
            for el in cycle :
                if (el < 0) :
                    direction = 2
                    break

            vertice_set = vertices[i]
            for vertex in vertice_set :
                vertex.size      = size
                vertex.direction = direction                
                vertex.cycle = i

    def clean_visit(self) :
        node = self.begin_node
        
        while node :
            node.visit = False
            node       = node.ap
            node.visit = False
            node       = node.ab



#####################################################################
################# REPRESENTS THE NODE OF A GRAPH  ###################
#####################################################################

class cycle_graph_node :
    #index  : stores the black edge i, 0 <= i <= n+1
    #value  : stores pi_i    
    #ap     : points to the record that stores pi_{i-1}, 1 <= i <= n+1
    #ab     : points to the record that stores pi_{i+1}, 0 <= i <= n
    #ac     : points to the record that stores i + 1,    0 <= i <= n
    #visit  : indicates if the record has been visited
    def __init__(self, index) :
        self.index, self.value    = index, 0
        self.cycle    = 0

        ## 0 = unset 
        ## i = num of black edges
        self.size        = 0 

        ## 0 = unset
        ## 1 = convergent
        ## 2 = divergent
        self.direction = 0

        self.ap, self.ab, self.ac = 0,0,0
        self.visit  = False


#####################################################################
## This class was developed to sort cycle-graph configurations using
# 4-cuts. We apply several rules in order to make the configuration 
# closer to the identity. If no rule is found, so the program returns
# with no error message. 
#####################################################################

class Cut_4 :
    def __init__(self, cycles) :
        self.input_cycles  = str(cycles).replace(" ", "")
        self.graph      = cycle_configuration_graph(cycles)
    
    def has_short_cycle(self, vertices) :
        for vertice in vertices :
            if vertice[0].size == 2 or vertice[0].size == 3 :
                return True
        return False
    
    def sort(self, arqout, inicio) :
        sequence = []
        graph = self.graph

        graph.calculate_cycles()
        _, vertices = graph.get_cycles(want_vertices = True)
        num_cycles = len(graph.get_cycles())
        short = self.has_short_cycle(vertices)
        
        lowerb = graph.lower_bound(num_cycles, short)
        max_approx = 1.5

        while True :
            graph.calculate_cycles()
            if len(sequence) == 0 :
                operations = self.search_k4cut(graph, True)
            else :
                operations = self.search_k4cut(graph, False)
            if not operations :
                break
            
            for op in operations :
                if callable(op[0]) :
                    operations = operations + op[0](*tuple(op[1:]))
                else :
                    op = self.format_node_to_operation(op)
                    graph.k4cut(*op)
                    sequence.append(tuple(op))

        graph.calculate_cycles()
        all_cycles, vertices = graph.get_cycles(want_vertices = True)
        num_cycles = len(all_cycles)

        if num_cycles == graph.n :
            if len(sequence) == 0 :
                lowerb = 1
            if (((float(len(sequence)))/lowerb) <= max_approx) :
                print('%d %s %f %f' % (#self.input_cycles,
                           len(sequence),
                           str(sequence).replace(" ", ""),
                           len(sequence)/lowerb,
                           time.time()-inicio), file = arqout)
            else :
                print('%s %s %s %d %s %f ERROR-LOWER-BOUND-HIGHER' % (self.input_cycles,
                               self.input_wgray,
                               self.input_wblack,
                               len(sequence),
                               str(sequence).replace(" ", ""),
                               len(sequence)/lowerb))
        else :
            print('%s %s %s %d %s %f ERROR-NOT-SORTED' % (self.input_cycles,
                               self.input_wgray,
                               self.input_wblack,
                               len(sequence),
                               str(sequence).replace(" ", ""),
                               len(sequence)/lowerb))

    def format_node_to_operation(self, operation) :
        if operation :
            i = int(math.ceil(float(operation[0].index + 1) / 2))
            j = int(math.ceil(float(operation[1].index + 1) / 2))
            k = int(math.ceil(float(operation[2].index + 1) / 2))
            l = int(math.ceil(float(operation[3].index + 1) / 2))
            return(i,j,k,l,operation[4])     

    def search_k4cut(self, graph, BF) :
        graph.calculate_cycles()
        operations = None
        
        ## 1 - Search a 4-cut on a long cycle that creates three new cycles
        if (BF) and (not operations) :
            operations = self.search_4cut_long_cycle(graph)
            if operations :
                return [operations]
        ## 2 - Search a convergent black edge crossing with another black edge
        if not operations :
            operations = self.search_4cut_convergent_black_edge(graph)
            if operations :
                return [operations]
        ## 3 - Search a long fully divergent cycle
        if not operations :
            operations = self.search_4cut_long_fully_divergent_cycle(graph)
            if operations :
                return [operations]
        ## 4 - Search two short fully divergent cycles that do not cross
        if not operations :
            operations = self.search_4cut_two_short_divergent_cycles(graph)
            if operations :
                return [operations]

##################################################
############### Several Steps ####################
##################################################

    def __list_complement_cycles(self, node) :
        list_of_vertices = [node, node.ap]
        aux = node.ap.ac
        while aux != node :
            list_of_vertices.append(aux)
            aux = aux.ap
            list_of_vertices.append(aux)
            aux = aux.ac
        
        new_ab = dict()
        aux = node.ap
        for _ in range(node.size-1) :
            aux2 = aux.ab
            while aux2.cycle != node.cycle :
                aux2 = aux2.ap.ab
            new_ab[aux.index] = aux2
            new_ab[aux2.index] = aux
            aux = aux2.ap
        new_ab[node.index] = aux
        new_ab[aux.index] = node
        
        cycles = []
        while len(list_of_vertices) > 0 :
            aux = list_of_vertices.pop()
            curr_cycle = [aux]
            if new_ab[aux.index].ap != aux :
                curr_cycle.append(new_ab[aux.index])
            list_of_vertices.remove(new_ab[aux.index])
            aux = new_ab[aux.index]
            while aux.ac != curr_cycle[0] :
                aux = aux.ac
                if aux.ap not in curr_cycle :
                    curr_cycle.append(aux)
                list_of_vertices.remove(aux)
                aux = new_ab[aux.index]
                if aux.ap not in curr_cycle :
                    curr_cycle.append(aux)
                list_of_vertices.remove(aux)
            cycles.append(curr_cycle)
        return cycles
        
    def __order_by_index(self, i,j,k,l) :
        if i.index > j.index:
            if j.index > k.index:
                if l.index > j.index:
                    if l.index > i.index:
                        return [k,j,i,l]
                    else:
                        return [k,j,l,i]
                else:
                    if l.index > k.index:
                        return [k,l,j,i]
                    else:
                        return [l,k,j,i]
            else:
                if i.index > k.index:
                    if l.index > k.index:
                        if l.index > i.index:
                            return [j,k,i,l]
                        else:
                            return [j,k,l,i]
                    else:
                        if l.index > j.index:
                            return [j,l,k,i]
                        else:
                            return [l,j,k,i]
                else:
                    if l.index > i.index:
                        if l.index > k.index:
                            return [j,i,k,l]
                        else:
                            return [j,i,l,k]
                    else:
                        if l.index > j.index:
                            return [j,l,i,k]
                        else:
                            return [l,j,i,k]
        else:
            if i.index > k.index:
                if l.index > i.index:
                    if l.index > j.index:
                        return [k,i,j,l]
                    else:
                        return [k,i,l,j]
                else:
                    if l.index > k.index:
                        return [k,l,i,j]
                    else:
                        return [l,k,i,j]
            else:
                if j.index > k.index:
                    if l.index > k.index:
                        if l.index > j.index:
                            return [i,k,j,l]
                        else:
                            return [i,k,l,j]
                    else:
                        if l.index > i.index:
                            return [i,l,k,j]
                        else:
                            return [l,i,k,j]
                else:
                    if l.index > j.index:
                        if l.index > k.index:
                            return [i,j,k,l]
                        else:
                            return [i,j,l,k]
                    else:
                        if l.index > i.index:
                            return [i,l,j,k]
                        else:
                            return [l,i,j,k]
        
    def __get_indices_sigma_from_complement_cycle(self, cycle, graph) :
        found = False

        sigmas = [  [[-1,-2,-3],[-1,-2,-3]] , [[-1,3,2],[-1,3,2]],
                    [[-1,-3,2],[-1,3,-2]]   , [[-1,3,-2],[-1,-3,2]],
                    [[2,1,-3],[2,1,-3]]     , [[2,-1,-3],[-2,1,-3]],
                    [[-2,1,-3],[2,-1,-3]]   , [[2,-3,1],[3,1,-2]],
                    [[-2,3,1],[3,-1,2]]     , [[2,-3,-1],[-3,1,-2]],
                    [[-2,3,-1],[-3,-1,2]]   , [[-2,-3,-1],[-3,-1,-2]],
                    [[3,1,-2],[2,-3,1]]     , [[3,-1,2],[-2,3,1]],
                    [[-3,-1,2],[-2,3,-1]]   , [[-3,1,-2],[2,-3,-1]],
                    [[-3,-1,-2],[-2,-3,-1]] , [[3,2,-1],[-3,2,1]],
                    [[3,-2,1],[3,-2,1]]     , [[-3,2,1],[3,2,-1]],
                    # I don't know if the following ones make sense, 
                    # but just to be sure that I included all possible cases.
                    [[-1,2,-3],[-1,2,-3]]   , [[-2,-3,1],[3,-1,-2]],
                    [[3,-1,-2],[-2,-3,1]]   , [[3,2,1],[3,2,1]],
                    [[-3,2,-1],[-3,2,-1]]   ]


        for aux in itertools.combinations(cycle, 4) :
            i, j, k, l = self.__order_by_index(*aux)
            
            if ( (i.ap == j or i.ap == k or i.ap == l) or
                 (j.ap == i or j.ap == k or j.ap == l) or
                 (k.ap == i or k.ap == j or k.ap == l) or
                 (l.ap == i or l.ap == j or l.ap == k)  ) :
                break
            #format_node_to_operation
            for sigma in sigmas :
                i2 = int(math.ceil(float(i.index + 1) / 2))
                j2 = int(math.ceil(float(j.index + 1) / 2))
                k2 = int(math.ceil(float(k.index + 1) / 2))
                l2 = int(math.ceil(float(l.index + 1) / 2))
                sigma2 = sigma[0]
                block_len = [j2-i2,k2-j2,l2-k2]
                i3 = i2
                j3 = i3+block_len[abs(sigma2[0])-1]
                k3 = j3+block_len[abs(sigma2[1])-1]
                l3 = l2
                sigma3 = sigma[1]
                curr_num_cycles = len(graph.get_cycles())
                graph.k4cut(i2,j2,k2,l2,sigma2)
                graph.calculate_cycles()
                if len(graph.get_cycles()) == curr_num_cycles + 3 :
                    found = True
                graph.k4cut(i3,j3,k3,l3,sigma3)
                if found :
                    return(i,j,k,l,sigma2)
        return None
        

    def search_4cut_long_cycle(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            if (cycle[0].size > 3) : ## Check if it is long
                isCycleValid = True
                if isCycleValid :
                    aux = cycle[0]
                    black_edges = [cycle[0]]
                    black_edges_indices = [cycle[0].index]
                    aux = cycle[0].ap.ac
                    while aux != cycle[0] :
                        black_edges.append(aux)
                        black_edges_indices.append(aux.index)
                        aux = aux.ap.ac
                    resultado = self.__get_indices_sigma_from_complement_cycle(black_edges, graph)
                    if resultado :
                        return resultado
        return None
        
    '''
    Types of gray edges:
    0 trivial
    1 convergent /_```_\
    2 divergent /_```\_
    3 convergent _/```\_
    4 divergent _/``_\
    '''
    def __check_gray_edge_type(self, graph, edge) :
        if edge.ac == edge.ap :
            if (  (edge.ab.index == graph.n) or 
                  (edge.index == edge.ab.index + 1)  ) :
                return (0, edge, edge.ap)
            return (0, edge.ap, edge)
        edgeL, edgeR = edge, edge.ac
        if edgeR.index < edgeL.index :
            edgeL, edgeR = edgeR, edgeL
        if (  (edgeL.ab == 0) or 
              (edgeL.index == edgeL.ab.index + 1)  ) : 
            if (edgeR.ab == 0) or (edgeR.index < edgeR.ab.index) :
                return (1, edgeL, edgeR)
            return (2, edgeL, edgeR)
        else:
            if ( (edgeR.ab == 0) or
                (edgeR.index + 1 == edgeR.ab.index)  ) :
                return (4, edgeL, edgeR)
            return (3, edgeL, edgeR)

    def __find_gray_edge_crossing_with_convergent(self, graph, etype, vLeft, vRight) :
        if etype == 1 :
            aux = vLeft.ap
            while aux != vRight :
                if aux.ac.index < vLeft.index or aux.ac.index > vRight.index :
                    otype, otherL, otherR = self.__check_gray_edge_type(graph, aux)
                    break
                aux = aux.ab
                if aux.ac.index < vLeft.index or aux.ac.index > vRight.index :
                    otype, otherL, otherR = self.__check_gray_edge_type(graph, aux)
                    break
                aux = aux.ap
        else :
            aux = vLeft.ab
            while aux != vRight :
                if aux.ac.index < vLeft.index or aux.ac.index > vRight.index :
                    otype, otherL, otherR = self.__check_gray_edge_type(graph, aux)
                    break
                aux = aux.ap
                if aux.ac.index < vLeft.index or aux.ac.index > vRight.index :
                    otype, otherL, otherR = self.__check_gray_edge_type(graph, aux)
                    break
                aux = aux.ab
         
        if otype == 1 or otype == 3 :
            sigma = [3,2,1]
            # first let us check if there are only three black edges
            if vLeft.ap == otherL :
                i = vLeft
                j = otherL
                k = vRight
                l = otherR
                if k.index > l.index :
                    k, l = l, k
            elif vLeft.ap == otherR :
                i = otherL
                j = vLeft
                k = otherR
                l = vRight
            elif vRight.ap == otherL :
                i = vLeft
                j = otherL
                k = vRight
                l = otherR
            elif vRight.ap == otherR :
                i = otherL
                j = vLeft
                k = vRight
                l = otherR
                if i.index > j.index :
                    i, j = j, i
            else : # four black edges
                i, j = vLeft, otherL
                if i.index > j.index :
                    i, j = j, i
                k, l = vRight, otherR
                if k.index > l.index :
                    k, l = l, k
        else : # otype = 2 or otype = 4
            if vLeft.ap == otherL :
                if otype == 2 :
                    i = otherL
                    j = vLeft
                    k = otherR
                    l = vRight
                    sigma = [-2,-3,1]
                else :
                    i = vLeft
                    j = otherL
                    k = vRight
                    l = otherR
                    sigma = [3,-2,1]
            elif vLeft.ap == otherR :
                i = otherL
                j = vLeft
                k = otherR
                l = vRight
                sigma = [-3,1,2]
            elif vRight.ap == otherL :
                i = vLeft
                j = vRight
                k = otherL
                l = otherR
                sigma = [3,-1,2]
            elif vRight.ap == otherR :
                if otype == 2 :
                    i = otherL
                    j = vLeft
                    k = otherR
                    l = vRight
                    sigma = [-2,1,3]
                else :
                    i = vLeft
                    j = otherL
                    k = vRight
                    l = otherR
                    sigma = [-1,-2,3]
            else : # four black edges
                i, j = vLeft, otherL
                if i.index > j.index :
                    i, j = j, i
                k, l = vRight, otherR
                if k.index > l.index :
                    k, l = l, k
                if i == otherL :
                    sigma = [-2,-3,1]
                else :
                    sigma = [3,-1,-2]
        return(i,j,k,l,sigma)


    def search_4cut_convergent_black_edge(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            if (cycle[0].size >= 2) : ## Check if it is not trivial
                aux = cycle[0].ap
                etype, vLeft, vRight = self.__check_gray_edge_type(graph, aux)
                if etype % 2 == 1 :
                    return self.__find_gray_edge_crossing_with_convergent(graph, etype, vLeft, vRight)
                aux = aux.ap.ac
                while aux != cycle[0].ap :
                    # find a convergent black edge
                    etype, vLeft, vRight = self.__check_gray_edge_type(graph, aux)
                    if etype % 2 == 1 :
                        return self.__find_gray_edge_crossing_with_convergent(graph, etype, vLeft, vRight)
                    aux = aux.ap.ac
        return None

    def search_4cut_long_fully_divergent_cycle(self, graph) :
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            if (cycle[0].size > 3) : ## Check if it is long
                aux = cycle[0]
                i = aux
                i2 = aux.index
                aux = aux.ap.ac
                j = aux
                j2 = aux.index
                aux = aux.ap.ac
                k = aux
                k2 = aux.index
                aux = aux.ap.ac
                l = aux
                l2 = aux.index
                if i2 > j2 > k2 > l2 :
                    return(l,k,j,i,[-1,2,-3])
                elif i2 > j2 > l2 > k2 :
                    return(k,l,j,i,[-1,3,-2])
                elif i2 > k2 > j2 > l2 :
                    return(l,j,k,i,[-2,3,-1])
                elif i2 > k2 > l2 > j2 :
                    return(j,l,k,i,[-3,1,-2])
                elif i2 > l2 > j2 > k2 :
                    return(k,j,l,i,[-2,1,-3])
                else : # i > l > k > j
                    return(j,k,l,i,[-3,2,-1])
        return None


    def search_4cut_two_short_divergent_cycles(self, graph) :
        cycle1 = None
        _, vertices = graph.get_cycles(want_vertices = True)
        for cycle in vertices :
            if (cycle[0].size > 1) :
                if (cycle1 == None) : ## Check if it is not trivial
                    cycle1 = cycle[0]
                    i = cycle1.ac
                    i2 = cycle1.ac.index
                    j = cycle1
                    j2 = cycle1.index
                else :
                    k = cycle[0].ac
                    k2 = cycle[0].ac.index
                    l = cycle[0]
                    l2 = cycle[0].index
                    if i2 < k2 < l2 < j2 : #between c1
                        return(i,k,l,j,[-3,2,-1])
                    elif k2 < i2 < j2 < l2 : #c1 is between
                        return(k,i,j,l,[-3,2,-1])
                    elif i2 < j2 < k2 < l2 : #c1 is between
                        return(i,j,k,l,[-1,2,-3])
                    elif k2 < l2 < i2 < j2 : #c1 is between
                        return(k,l,i,j,[-1,2,-3])
        # there is only one non-trivial cycle...
        if cycle1 :
            return(i,i,j,j,[1,-2,3])
        return None


#############################################################################
###   BASIC FUNCTIONS, TO TRANSFORM A PERMUTATION INTO A LIST OF CYCLES   ###
#############################################################################

def get_position(permutation) :
    n = len(permutation)-2
    position    = [-1 for i in range(0, n+2)]
    for i in range(0, n+2) :
        position[abs(permutation[i])] = i
    return position

def get_rightmost_element(cycle, position) :
    max_position = 0
    for i in range(len(cycle)) :
        if position[cycle[i]] > position[cycle[max_position]] :
            max_position = i
    return max_position

## The unordered cycle starts with a gray edge, we order them by
## making it start with the rightmost black edge.
def order_cycle(cycle, position) :
    index = get_rightmost_element(cycle, position)
    new   = []
    new.append(cycle[index])

    if index % 2 == 0 :
        iter_el  = (index-1) % len(cycle)
        while iter_el != index :
            new.append(cycle[iter_el])
            iter_el = (iter_el-1) % len(cycle)
    else :
        iter_el  = (index+1) % len(cycle)
        while iter_el != index :
            new.append(cycle[iter_el])
            iter_el = (iter_el+1) % len(cycle)
    return new

def canonical_representation(cycle, position) :
    cycle     = order_cycle(cycle, position)
    canonical = []

    for i in range(0,len(cycle),2) :
        if position[cycle[i]] < position[cycle[i+1]] :
            black = -position[cycle[i+1]]
            canonical.append(black )
        else :
            black = position[cycle[i]]
            canonical.append(black)
    return canonical


def construct_str_cycle(permutation) :
    n = len(permutation)

    permutation = [0] + permutation + [n+1]
    position    = [-1 for i in range(0, n+2)]
    #sign        = [-1 for i in range(0, n+2)]


    for i in range(1, n+2) :
        position[abs(permutation[i])] = i
    #    sign    [abs(permutation[i])] = permutation[i] / abs(permutation[i])

    ## 1 if the gray edge i,-(i+1) was never used.
    gray_available     = [1 for i in range(0, n+1)]
    #black_available    = [1 for i in range(0, n+1)]

    cycles = []

    for i in range(0, n+1) :

        ## 
      if gray_available[i] :

        start     = i
        cycle = [start]

        end   = start
        positive  = True
        
        while True :
            
            ## Will be used later, it says if after walking through
            ## the black edge we are in the right or in the left
            is_vertex_left = None

            if positive :
                ## Gray edge: we are looking for the edge ( end,-(end+1) )
                gray_available[end] = gray_available[end] - 1
                end = end + 1
                cycle.append(end)

                ## Black edge: we are at the vertex -end.
                if permutation[position[end]] > 0 :
                    # If the sign in that position is positive, than
                    # -end is in the left (cycle graph)                    
                    end = abs(permutation[position[end]-1])
                    is_vertex_left = False
                
                else :
                    # If the sign in that position is negative, than
                    # -end is in the right (cycle graph)
                    end = abs(permutation[position[end]+1])
                    is_vertex_left = True
            else :
                ## Gray edge: we are looking for the edge ( -end, end-1  )
                end = end - 1                                 ##  Note we swapped
                gray_available[end] = gray_available[end] - 1 ##  these lines
                cycle.append(end)


                ## Black edge: we are at the vertex +end.
                if permutation[position[end]] > 0 :
                    # If the sign in that position is positive, than
                    # +end is in the right (cycle graph)
                    end = abs(permutation[position[end]+1])
                    is_vertex_left = True
                else : 
                    # If the sign in that position is negative, than
                    # +end is in the left (cycle graph)
                    end = abs(permutation[position[end]-1])
                    is_vertex_left = False
                    
            if end == start :
                break
            else :
                cycle.append(end)
                
                if is_vertex_left :
                    if permutation[position[end]] < 0 :
                        positive = True
                    else :
                        positive = False
                else :                    
                    if permutation[position[end]] < 0 :
                        positive = False
                    else :
                        positive = True
        cycles.append(cycle)

    int_position = get_position(permutation)
    canonicals = []

    for cycle in cycles :
        canonicals.append(canonical_representation(cycle, int_position))

    str_canonicals = str(canonicals)
    str_canonicals = str_canonicals.replace(" ", "")
    str_canonicals = str_canonicals.replace("[", "(")
    str_canonicals = str_canonicals.replace("]", ",)")
    return str_canonicals

        
if __name__ == '__main__':

    inputPath = sys.argv[1]
    outputPath = sys.argv[2]

    arqperms = open(inputPath, 'r')
    arqout = open(outputPath, 'a')
        
    while True:
        line = arqperms.readline()
        if not line :
            break
        perm, ops = line.split(' ')
        permutation = eval("[%s]" % perm)
        str_cycles = construct_str_cycle(permutation)

        config = eval(str_cycles)
        sort = Cut_4(config)
        sort.sort(arqout, time.time())
        
    arqperms.close()
    arqout.close()

