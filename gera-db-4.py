import random
import sys

# all possible sigmas for a 4-cut. we do not use 1,2,3 since it always
# results in the same permutation
allS3 = [
    [1,2,-3],[1,-2,3],[-1,2,3],[1,-2,-3],[-1,-2,3],[-1,2,-3],[-1,-2,-3],#[1,2,3]
    [1,3,2],[1,3,-2],[1,-3,2],[-1,3,2],[1,-3,-2],[-1,-3,2],[-1,3,-2],[-1,-3,-2],
    [2,1,3],[2,1,-3],[2,-1,3],[-2,1,3],[2,-1,-3],[-2,-1,3],[-2,1,-3],[-2,-1,-3],
    [2,3,1],[2,3,-1],[2,-3,1],[-2,3,1],[2,-3,-1],[-2,-3,1],[-2,3,-1],[-2,-3,-1],
    [3,1,2],[3,1,-2],[3,-1,2],[-3,1,2],[3,-1,-2],[-3,-1,2],[-3,1,-2],[-3,-1,-2],
    [3,2,1],[3,2,-1],[3,-2,1],[-3,2,1],[3,-2,-1],[-3,-2,1],[-3,2,-1],[-3,-2,-1],
    ]

# the custom seed is mandatory
custom_seed = eval(sys.argv[1])
random.seed(custom_seed)

# desired size of each permutation
n = eval(sys.argv[2])
# desired number of permutations
m = eval(sys.argv[3])
# desired number of 4-cuts
d = eval(sys.argv[4])

# The following line was added by mistake. However,
# to generate the exact data set used in the article
# it cannot be removed, so just ignore it
iota = [random.randint(0, 100) for a in range(n+2)]


for b in range(m):
    iota = list(range(0,n+2))
    ops = []
    for kcut in range(d):
        ijkl = random.choices(range(1, n+1), k=4)
        ijkl.sort()
        i,j,k,l = ijkl
        indiceS3 = allS3[random.randint(0,46)]
        blocos = [iota[:i], iota[i:j], iota[j:k], iota[k:l], iota[l:]]
        iota = blocos[0]
        for el in indiceS3 :
            if el < 0 :
                iota += [-revel for revel in blocos[-el][::-1]]
            else :
                iota += blocos[el]
        iota += blocos[-1]
        ops.append([i,j,k,l,indiceS3])
    # We print the permutation followed by the list of 4-cuts applied
    print(str(iota[1:-1]).replace(', ',',').replace('[','').replace(']',''), str(ops).replace(', ',','))


