import random
import sys


# the custom seed is mandatory
custom_seed = eval(sys.argv[1])
random.seed(custom_seed)
random.randint(0, 100)

# desired size of each permutation
n = eval(sys.argv[2])
# desired number of permutations
m = eval(sys.argv[3])
# desired number of 4-cuts
d = eval(sys.argv[4])
# number of cuts
k_size = eval(sys.argv[5])

for b in range(m):
    iota = list(range(0,n+2))
    ops = []
    for kcut in range(d):
        cuts = None
        while not cuts :
            cuts = sorted(random.choices(list(range(1,n+2)),k=k_size))
            i = 0
            while len(cuts) > 1 and i < len(cuts)-1 :
                if cuts[i] != cuts[i+1] :
                    i += 1
                else :
                    cuts[i:i+2] = [cuts[i]]
            if len(cuts) < 2 :
                cuts = None
        perm_cut = list(range(2,len(cuts)+1))
        iota_perm_cut = perm_cut[:]
        signs = random.choices([1,-1], weights = [1, 1], k=len(perm_cut))
        signs[0] = signs[-1] = 1
        perm_cut = [perm_cut[i]*signs[i] for i in range(len(perm_cut))]
        random.shuffle(perm_cut)
        perm_cut = [1] + perm_cut + [len(cuts)+1]
        # just to avoid an useless k-cut,when sigma is the identity permutation,
        # we set all signs to negative in that case (-1,-2,...)
        if perm_cut[1:-1] == iota_perm_cut :
            signs = random.choices([-1,-1],k=len(perm_cut))
            perm_cut = [perm_cut[i]*signs[i] for i in range(len(perm_cut))]
            perm_cut = [1] + perm_cut[1:-1] + [len(cuts)+1]
        blocks = [iota[:cuts[0]]] + [iota[cuts[i-1]:cuts[i]] for i in range(1,len(cuts))] + [iota[cuts[-1]:]]
        iota = []
        for el in perm_cut :
            if el > 0 :
                iota += blocks[el-1]
            else :
                iota += [-x for x in blocks[abs(el+1)][::-1]]
        while len(cuts) < k_size :
            cuts.append(cuts[-1])
        while len(perm_cut) <= k_size :
            perm_cut.append(len(perm_cut)+1)
        ops.append(list([cuts,perm_cut]))
    print(str(iota[1:-1]).replace(', ',',').replace('[','').replace(']',''), str(ops).replace(', ',','))


