import math
import sys
import time

class create_and_sort_configurations :

    def smallest_hcycle(self, hcycles) :
        minimum = 0
        for count in range(len(hcycles)) :
            if len(hcycles[count]) < len(hcycles[minimum]) :
                minimum = count
        return hcycles[minimum]


    def renumber_canonical_configuration(self, configuration) :
        index = []
        for cycle in configuration :
            for el in cycle :
                index.append(abs(el))
    
        index.sort()
        old_to_new = {}
        for i in range(len(index)) :
            old_to_new[index[i]] = i+1
    
        new = []
        for cycle in configuration :
            ncycle = []
            for el in cycle :
                if el < 0 :
                    ncycle.append(-old_to_new[-el])
                else :
                    ncycle.append(old_to_new[el])
            new.append(ncycle)
        return new, index


    def find_extension(self, renumbered, renumbered_index, extensions, original_configuration_size, k_value) :

        renumbered_graph = cycle_configuration_graph(renumbered)
        hcycles          = renumbered_graph.hamiltonian_cycles()
        hcycle           = self.smallest_hcycle(hcycles)
    
        ## This provides an index from the hamiltonian cycle for the
        ## subset to the full configuration
        hcycle_index = {}
        for i in range(1, len(renumbered_index)+1) :
            hcycle_index[2*(i-1)] = 2*(renumbered_index[i-1]-1)
            hcycle_index[2*(i)-1] = 2*(renumbered_index[i-1])-1

        in_set = []
        for i in range(0, len(hcycle), 2) :
            el = [hcycle[i], hcycle[i+1]]
            if hcycle[i] > hcycle[i+1] :
                el = [hcycle[i+1], hcycle[i]]
            
            if el[0] == 0 :
                ## Here, we break the region into two. One for the
                ## beggining of the cycle graph and the other for the
                ## end.
                el1 = [0,0]
                el2 = [0,0]
    
                el1[0] = 0
                el1[1] = hcycle_index[el[0]]
    
                el2[0] = hcycle_index[el[1]]
                el2[1] = 2*original_configuration_size
                in_set.append(el1)
                in_set.append(el2)
            else :
                el[0] = hcycle_index[el[0]]
                el[1] = hcycle_index[el[1]]
                in_set.append(el)

        if len(hcycle) == 2*renumbered_graph.n :
            return self.full_extension(in_set, extensions)
        else :
            return self.partial_extension(in_set, extensions)
 
    def full_extension(self, in_set, extensions) :
        cycle_to_extend = None
        count = 0
        
        while count < len(extensions) :
            cycle    = extensions[count]
            if len(cycle) > 1 :
                cycle_to_extend = cycle
                break
            count = count + 1
        return cycle_to_extend, count

    def partial_extension(self, in_set, extensions) :    
        ## Ja verifiquei os intervalos que estao dentro e for do
        ## ciclo hamiltoniano que eu quero tratar. Agora eu
        ## extendo com qualquer ciclo que tem aresta preta dentro
        ## e fora do intervalo.

        cycle_to_extend = None
        count = 0
        while count < len(extensions) :

            cycle = extensions[count]
            in_cycle  = False
            out_cycle = False
            for el in cycle : 
                aux_in  = False
                for in_range in in_set :
                    if in_range[0] < 2*abs(el)-1 < in_range[1] :
                        aux_in = True
                if aux_in :
                    in_cycle  = True
                else :
                    out_cycle = True

            if in_cycle and out_cycle :
                cycle_to_extend = cycle
                break
            count = count + 1
        return cycle_to_extend, count
    

    def extend_and_sort(self, graph, k) :
        cycles = graph.get_cycles()

        ## 1 - Achar um ciclo pelo qual comecar o processo de extensao e consulta
        configuration = None
        
        i = 0
        while i < len(cycles) and len(cycles[i]) < 2 :
            i += 1
        if i == len(cycles) :
            return None

        j = 0
        while j < len(cycles) and len(cycles[j]) < 3 :
            j += 1
        if j < len(cycles) :
            i = j


        configuration = [cycles[i]]
        curr_k = len(configuration[0])

        extensions = list(cycles) 
        extensions.__delitem__(i)
        sorting    = []
        

        ## Dentro do proprio loop eu verifico se tem sorting e se dah
        ## para extender.
        while curr_k < k - 1 :
            ## Transformar o nosso conjunto de ciclos em uma
            ##  configuracao valida.
            renumbered, renumbered_index = self.renumber_canonical_configuration(configuration)
            
            # Estendemos com mais um ciclo...
            # A ordem não sei se importa. Começo com convergentes
            cycle_to_extend, count = self.find_extension(renumbered, renumbered_index, extensions, graph.n, k)
            if not cycle_to_extend :
                break        
            
            ## Ainda eh possivel estender, entao vamos nos preparar
            ## para a proxima iteracao do laco.
            curr_k = curr_k + len(cycle_to_extend)
            configuration = configuration + [cycle_to_extend]
            extensions.__delitem__(count)

        renumbered, renumbered_index = self.renumber_canonical_configuration(configuration)
        return self.config_to_kcut(renumbered, renumbered_index, k)

    
    def config_to_kcut(self, renumbered, renumbered_index, k) :
        new_perm = []
        ignorado = None
        hciclos = sorted(cycle_configuration_graph(renumbered).hamiltonian_cycles())
        for hciclo in hciclos[::-1] :
            for i in range(1,len(hciclo),2) :
                if hciclo[i] > 0 :
                    if hciclo[i] % 2 == 0 :
                        new_perm.append(int(hciclo[i]/2))
                    else :
                        new_perm.append(-int((hciclo[i]+1)/2))
        for i in range(len(new_perm)) :
            if new_perm[i] > 0 :
                new_perm[i] += 1
            else :
                new_perm[i] -= 1
        new_perm = [1] + new_perm + [len(new_perm)+2]
        # caso temos uma permutação na configuração
        if len(hciclos) == 1 :
            # caso seja possível cortar todas as arestas pretas, vamos
            # apenas organizar as coisas...
            if len(renumbered_index) <= k :
                while len(new_perm) <= k :
                    new_perm.append(len(new_perm)+1)
                indices = [int(i) for i in renumbered_index]
                while len(indices) < k :
                    indices.append(indices[-1])
            # caso contrário, vamos ignorar a última aresta preta
            # deixando o último elemento em sua posição atual
            else :
                ignorado = len(new_perm) - 1
                new_perm_ig = [el for el in new_perm if abs(el) != ignorado]
                new_perm_ig[-1] = new_perm_ig[-1]-1
                indices = [int(i) for i in renumbered_index[:-1]]
                new_perm = new_perm_ig
        # agora não temos uma permutação. primeiro vamos tratar
        # o caso onde um 2-ciclo irá permanecer
        else :
            par1 = None
            par2 = None
            par3 = None
            if len(hciclos) == 2 :
                tamanho = len(hciclos[0]) + len(hciclos[1])
                for i in range(1,tamanho,2) :
                    if ( (i in hciclos[0] and (i-1) in hciclos[1]) or
                         (i in hciclos[1] and (i-1) in hciclos[0]) ) :
                        if i in hciclos[1] :
                            hciclos = [hciclos[1], hciclos[0]]
                        par1 = [i-1,i]
                        posicao_i = hciclos[0].index(i)
                        posicao_j = hciclos[1].index(i-1)
                        if posicao_i == len(hciclos[0])-1 :
                            par2 = [hciclos[0][0],-1]
                        elif posicao_i == 0 :
                            par2 = [hciclos[0][len(hciclos[0])-1],-1]
                        elif posicao_i % 2 == 0 :  #par, -1
                            par2 = [hciclos[0][posicao_i-1],-1]
                        else :  #impar, +1
                            par2 = [hciclos[0][posicao_i+1],-1]
                        if posicao_j == len(hciclos[1])-1 :
                            par2[1] = hciclos[1][0]
                        elif posicao_j == 0 :
                            par2[1] = hciclos[1][len(hciclos[1])-1]
                        elif posicao_j % 2 == 0 :  #par, -1
                            par2[1] = hciclos[1][posicao_j-1]
                        else :  #impar, +1
                            par2[1] = hciclos[1][posicao_j+1]
                        break
            else :
                tamanho = len(hciclos[0]) + len(hciclos[1]) + len(hciclos[2])
                if 0 in hciclos[1] :
                    hciclos = [hciclos[1], hciclos[0], hciclos[2]]
                elif 0 in hciclos[2] :
                    hciclos = [hciclos[2], hciclos[0], hciclos[1]]
                for i in range(0,tamanho-1,2) :
                    if ( (i in hciclos[0] and (i+1) in hciclos[1]+hciclos[2]) or
                         ((i+1) in hciclos[0] and i in hciclos[1]+hciclos[2]) ) :
                        if i in hciclos[1] :
                            hciclos = [hciclos[1], hciclos[0], hciclos[2]]
                        elif i in hciclos[2] :
                            hciclos = [hciclos[2], hciclos[0], hciclos[1]]
                        if (i+1) in hciclos[2] :
                            hciclos = [hciclos[0], hciclos[2], hciclos[1]]
                        par1 = [i,i+1]
                        posicao_i = hciclos[0].index(i)
                        posicao_j = hciclos[1].index(i+1)
                        if posicao_i == len(hciclos[0])-1 :
                            par2 = [hciclos[0][0],-1]
                        elif posicao_i == 0 :
                            par2 = [hciclos[0][len(hciclos[0])-1],-1]
                        elif posicao_i % 2 == 0 :  #par, -1
                            par2 = [hciclos[0][posicao_i-1],-1]
                        else :  #impar, +1
                            par2 = [hciclos[0][posicao_i+1],-1]
                        if posicao_j == len(hciclos[1])-1 :
                            par2[1] = hciclos[1][0]
                        elif posicao_j == 0 :
                            par2[1] = hciclos[1][len(hciclos[1])-1]
                        elif posicao_j % 2 == 0 :  #par, -1
                            par2[1] = hciclos[1][posicao_j-1]
                        else :  #impar, +1
                            par2[1] = hciclos[1][posicao_j+1]
                        break
                par2, par3 = [par2[0], hciclos[2][0]], [par2[1], hciclos[2][-1]]
            brancas = [[0,0],]
            for j in range(len(hciclos)) :
                for i in range(0,len(hciclos[j])-1,2) :
                    if hciclos[j][i] != 0 and hciclos[j][i] != (tamanho-1) :
                        brancas.append([hciclos[j][i],hciclos[j][i+1]])
            brancas.append([tamanho-1,tamanho-1])
            cinzas = []
            for j in range(len(hciclos)) :
                cinzas.append([hciclos[j][0],hciclos[j][-1]])
                for i in range(1,len(hciclos[j])-2,2) :
                    cinzas.append([hciclos[j][i],hciclos[j][i+1]])
            pretas = [[0,0]]
            for i in range(1,tamanho-1,2) :
                pretas.append([i,i+1])
            pretas.append([tamanho-1,tamanho-1])
            nova_ordem = [[0,0]]
            brancas.remove(nova_ordem[0])
            while len(brancas) > 0 :
                if ( (par1 and nova_ordem[-1][1] in par1) or
                     (par2 and nova_ordem[-1][1] in par2) or
                     (par3 and nova_ordem[-1][1] in par3) ) :
                    if (par2 and nova_ordem[-1][1] in par2) :
                        par2, par1 = par1, par2
                    if (par3 and nova_ordem[-1][1] in par3) :
                        par3, par1 = par1, par3
                    if nova_ordem[-1][1] == par1[1] :
                        prox_adj = par1[0]
                    else :
                        prox_adj = par1[1]
                    i = 0
                    while prox_adj not in brancas[i] :
                        i += 1
                    if prox_adj == brancas[i][0] :
                        nova_ordem.append(brancas[i])
                        brancas.remove(brancas[i])
                    else :
                        nova_ordem.append(brancas[i][::-1])
                        brancas.remove(brancas[i])
                    par1 = None
                else :
                    i = 0
                    while nova_ordem[-1][1] not in cinzas[i] :
                        i += 1
                    if nova_ordem[-1][1] == cinzas[i][0] :
                        prox_adj = cinzas[i][1]
                    else :
                        prox_adj = cinzas[i][0]
                    i = 0
                    while prox_adj not in brancas[i] :
                        i += 1
                    if prox_adj == brancas[i][0] :
                        nova_ordem.append(brancas[i])
                        brancas.remove(brancas[i])
                    else :
                        nova_ordem.append(brancas[i][::-1])
                        brancas.remove(brancas[i])
            perm_bloco = [1]
            for el in nova_ordem[1:-1] :
                if el in pretas :
                    perm_bloco.append(pretas.index(el)+1)
                else :
                    perm_bloco.append(-pretas.index(el[::-1])-1)

            perm_bloco.append(len(nova_ordem))
            indices = [int(i) for i in renumbered_index]
            new_perm = perm_bloco[:]
            
            if len(indices) <= k :
                while len(new_perm) <= k :
                    new_perm.append(len(new_perm)+1)
                while len(indices) < k :
                    indices.append(indices[-1])

            else :
                encontrou = False
                indice = 1
                while not(encontrou) :
                    if indice >= len(new_perm) :
                        sys.exit(0)
                    if new_perm[indice] != new_perm[indice-1] + 1 :
                        indice += 1
                    else :
                        encontrou = True
                        if new_perm[indice] < 0 :
                            ignorado = -new_perm[indice]
                        else :
                            ignorado = new_perm[indice-1]
                        new_perm_ig = []
                        for el in new_perm :
                            if abs(el) != abs(ignorado) :
                                if abs(el) < abs(ignorado) :
                                    new_perm_ig.append(el)
                                else :
                                    if el > 0 :
                                        new_perm_ig.append(el-1)
                                    else :
                                        new_perm_ig.append(el+1)
                if not ignorado :
                    ignorado = 99999999
                indices = [int(i) for i in renumbered_index]
                indices = indices[:ignorado-1] + indices[ignorado:]
                new_perm = new_perm_ig

        return(indices, new_perm)
 

############################################################################
############### Do not need to sort a real permutation #####################
############################################################################

## This class was developed to sort cycle-graph configurations. We
## apply several rules in order to make the configuration closer to
## the identity. If no rule is found, so the program returns with no
## error message. 

## If the input cycle configuration is a full component, then we
## guarantee that the final permutation is the identity.

class k_cut_generic :
    def __init__(self, cycles, k_value) :
        self.input = str(cycles).replace(" ", "")
        if type(cycles) == str :
            cycles           = eval(cycles)
        self.k               = k_value
        self.graph           = cycle_configuration_graph(cycles)
        self.graph1          = cycle_configuration_graph(cycles)
        self.config           = create_and_sort_configurations()

    def sort(self, perm, inicio, saida) :
        sequence = []
        graph = self.graph
        graph.calculate_cycles()     
        k_value = self.k
        initial_cycles = len(graph.get_cycles())
        min1c = float(graph.n - graph.num_cycles())/float(k_value-1)
        min2c = float(math.ceil(min1c))
        if k_value % 2 == 0:
            min1o = float(graph.n - graph.num_odd_cycles())/float(k_value)
        else :
            min1o = float(graph.n - graph.num_odd_cycles())/float(k_value-1)
        min2o = float(math.ceil(min1o))

        while True :

            operations = None
            if not operations :
                graph.transform_into_simple_permutation()
                graph.calculate_cycles()
                cycles      = graph.get_cycles()
                operations  =  self.config.extend_and_sort(graph, k_value) 
                if not operations :
                    break
                igs, pos = graph.kcut(operations[0],operations[1],k_value)
                if igs :
                    sequence.append(tuple([igs[:],pos[:]]))
                    self.graph1.kcut(igs, pos, k_value)
                cycles = self.graph1.get_cycles()    
                graph  = cycle_configuration_graph(cycles)

            if not operations :
                break

        if graph.n != len(graph.get_cycles()) :
            print("%d %s %.3f %.3f %.3f %.3f %f ERROR-NOT-SORTED" % (len(sequence), str(sequence).replace(" ", ""), len(sequence)/min1c, len(sequence)/min2c, len(sequence)/min1o, len(sequence)/min2o, time.time()-inicio), file = saida)
        else :
            print("%d %s %.3f %.3f %.3f %.3f %f" % (len(sequence), str(sequence).replace(" ", ""), len(sequence)/min1c, len(sequence)/min2c, len(sequence)/min1o, len(sequence)/min2o, time.time()-inicio), file = saida)


############################################################################
############ Do not need to represent a real permutation ###################
############################################################################

## This class represents an unoriented cycle-graph. Keep in mind that
## this class is not the same as used for the sorting by transposition
## problem. Here, the white edges are not the reverse of the black
## edges as before. That is because the black edges can go from righ
## to left as well as to left to right.

## For instance, assume a permutation (5, 2, ...). The configuration
## in the graph would be.

## node(0).ap  = -5      ## node(0).ab  = Null
## node(-5).ab = 5       ## node(-5).ap = 0  
## node(5).ap  = -2      ## node(5).ab  = -5 
## node(-2).ab = 2       ## node(-2).ap = 5

class cycle_configuration_graph() :
    def __init__(self, cycles) : ## Ignore tells us which black edges
                                 ## we should ignore.

        ## The number of cycles was never calculated
        self.__num_cycles           = False
        self.__num_odd_cycles       = False

        ## self.n is the number of black edges. Remember this graph
        ## might not be a permutation
        self.n = 0
        for cycle in cycles :
            self.n = self.n + len(cycle)
        n = self.n
          
        # Creating nodes
        node_list = []
        node_list = [cycle_graph_node(i, False) for i in range(2*n)]
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
                    front = int(2*( cycle[i]  -1 ))
                else :
                    front = int(2*(-cycle[i]) -1)

                if cycle[j] > 0 :
                    back = int(2*( cycle[j]) -1)
                else :
                    back = int(2*(-cycle[j]  -1 ))

                node_list[front].ac = node_list[back]
                node_list[back].ac  = node_list[front]




############################################################                
################ Rearrangement Operations  #################
############################################################                

    def kcut(self, indices, novas_posicoes, k) :
        nodes  = [[None,None] for i in range(k)]
        unps   = [     0      for i in range(k)]
        pos_unps = novas_posicoes[:]
        count, unp_count = 0, 0

        node = self.begin_node

        while node :
            count = count + 1
            if not node.padded :
                unp_count = unp_count + 1

            if count in indices :
                pos = indices.index(count)
                nodes[pos] = [node, node.ap]
                unps[pos]  = unp_count
            node = node.ap.ab
        
        blocks = [[0,nodes[0][0]]] + [[nodes[i-1][1],nodes[i][0]] for i in range(1,k)] + [[nodes[-1][1],0]]

        pos_orig = novas_posicoes[:]
        for i in range(k+1) :
            if novas_posicoes[i] < 0 :
                novas_posicoes[i] = -novas_posicoes[i]
                blocks[novas_posicoes[i]-1] = blocks[novas_posicoes[i]-1][::-1]

        ordered_blocks = []
        for el in novas_posicoes :
            ordered_blocks.append(blocks[el-1])

        for i in range(len(ordered_blocks)-1) :
            if ordered_blocks[i][1] and ordered_blocks[i+1][0] :
                ordered_blocks[i][1].ap = ordered_blocks[i+1][0]
                ordered_blocks[i+1][0].ap = ordered_blocks[i][1]

        while len(unps) > 0 and unps[0] == 0 :
            unps = unps[1:]
            if pos_unps.count(1) > 0 :
                pos_unps.remove(1)
            else :
                pos_unps.remove(-1)
            pos_unps = [el-1 if el > 0 else el+1 for el in pos_unps]
        while len(unps) > 0 and unps[-1] == 0 :
            unps = unps[:-1]
            el = len(pos_unps)
            if pos_unps.count(el) > 0 :
                pos_unps.remove(el)
            else :
                pos_unps.remove(-el)
        if len(unps) > 1 and unps[0] != unps[-1] :
            tamanho = len(unps)
            unps2 = unps[:]
            removedp = 0
            for i in range(1,len(unps2)) :
                if unps2[i] == unps2[i-1] :
                    if pos_unps.count(i+1-removedp) > 0 :
                        pos_unps.remove(i+1-removedp)
                    else :
                        pos_unps.remove(-(i+1-removedp))
                    pos_unps = [el if abs(el) < (i+1-removedp) else (el-1 if el > 0 else el+1) for el in pos_unps]
                    unps = unps[:i-removedp] + unps[i+1-removedp:]
                    removedp += 1
            while len(unps) < k :
                unps.append(unps[-1])
            while len(pos_unps) <= k :
                pos_unps.append(len(pos_unps)+1)

        else :
            unps, pos_unps = None, None

        self.__num_cycles           = False
        self.__num_odd_cycles       = False
        self.reset_indices()
        self.get_cycles()

        return unps, pos_unps

############################################################                
############## Cycle Graph Transformations  ################
############################################################                
    def transform_into_simple_permutation(self) :
        node = self.end_node                
        while node :
            b3 = node
            b1 = b3.ap.ac
            b2 = b1.ap.ac
            vg  = b2.ap.ac

            if b3 != b1 and b3 != b2 and b3 != vg :
                ## The cycle has more than 3 black edges
                v = cycle_graph_node(None, True)
                w = cycle_graph_node(None, True)

                ## Setting edges for v and w
                w.ab = v
                v.ab = w
                w.ap = b3.ap
                v.ap = b3
                w.ac = b2.ap
                v.ac = vg 

                #Setting other pointers
                b3.ap.ap = w
                b3.ap    = v
                b2.ap.ac = w
                vg.ac    = v
            else :
                node = node.ap.ab

        ## Fazer esse aqui
        self.__num_cycles           = False
        self.__num_odd_cycles       = False
        self.reset_indices()

        self.n = 0
        node = self.begin_node
        while node :
            self.n += 1
            node = node.ap.ab
        


############################################################                
################### Auxiliary Methods  #####################
############################################################                

    def hamiltonian_cycles(self) :
        self.clean_visit()

        outer_node = self.begin_node
        hamiltonian_cycles = []

        while outer_node :
            if not outer_node.visit :
                cycle = []

                node = outer_node
                node.visit = True
                cycle.append(node.index)

                node = node.ac
                node.visit = True
                cycle.append(node.index)

                while (node.ab  and not node.ab.visit) :
                    node = node.ab
                    node.visit = True
                    cycle.append(node.index)
                    node = node.ac
                    node.visit = True
                    cycle.append(node.index)

                cycle = cycle[1:] + [cycle[0]]
                hamiltonian_cycles.append(cycle)
            outer_node = outer_node.ap.ab
        return hamiltonian_cycles


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
            vertice_set = vertices[i]
            for vertex in vertice_set :
                vertex.cycle = i

    def clean_visit(self) :
        node = self.begin_node
        
        while node :
            node.visit = False
            node       = node.ap
            node.visit = False
            node       = node.ab


#####################################################################
################## REPRESENTS A NODE OF A GRAPH #####################
#####################################################################

class cycle_graph_node :
    #index  : stores the black edge i, 0 <= i <= n+1
    #value  : stores pi_i    
    #ap     : points to the record that stores pi_{i-1}, 1 <= i <= n+1
    #ab     : points to the record that stores pi_{i+1}, 0 <= i <= n
    #ac     : points to the record that stores i + 1,    0 <= i <= n
    #visit  : indicates if the record has been visited
    def __init__(self, index, padded) :
        self.index, self.value    = index, 0        
        self.padded   = padded
        self.cycle    = 0 
        self.ap, self.ab, self.ac = 0,0,0
        self.visit  = False

        ## 0 = unset 
        ## i = num of black edges
        self.size        = 0

#####################################################################

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
    sign        = [-1 for i in range(0, n+2)]


    for i in range(1, n+2) :
        position[abs(permutation[i])] = i
        sign    [abs(permutation[i])] = permutation[i] / abs(permutation[i])

    ## 1 if the gray edge i,-(i+1) was never used.
    gray_available     = [1 for i in range(0, n+1)]
    black_available    = [1 for i in range(0, n+1)]

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
    k_value = int(sys.argv[3])

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

        sort = k_cut_generic(config, k_value)
        sort.sort(permutation, time.time(), arqout)
    arqperms.close()
    arqout.close()


