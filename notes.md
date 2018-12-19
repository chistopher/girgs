

# Improvement Ideas
- [ ] use level local index for cells
- [x] in WeightLayer: use array_view concept to provide access to all nodes in ViA
    - e.g. only one method that returns pointer to first element and size of collection
    - removes need for multiple calls to k-th point
- [x] only save one directed edge instead of both directions
- [x] use omp to process recursive calls in parallel (somehow aggregate on edge vectors?
    - circumvented aggregation by grouping calls for parallel execution
    - possible because of not saving both directions for edges
- [x] produce deterministic results with same seed even when parallel
    - needs static scheduling and a random generators for each thread
- [x] stop building empty weight_layers
- [x] set levels to (L-2)/D+1 ? because we need no coord_helper if we stop recursion after (L-2)/D
    - does not work because we also need coord helper for node insertion
- [x] make graph member to reduce stack size
    - removed graph from recursion
- [x] move all index helper functions to coordinate helper
- [ ] pre-compute expressions used in for loop conditions
- [x] checkEdgeExplicit make alpha == inf check return
    - prevents a call to random generator
- [x] remove type 2 samples in threshold model
- [x] remove bruteforce over all ij for type1
- [x] remove bruteforce over all ij for type2
    - it is correct to sample all pairs in m_layer_pairs[k] for k in range(level, m_levels)
    - doing so gives no perf. improvement over current implementation but is less code
- [x] remove computations for assertions in release builds


# TODO's
- [ ] look at average size of the V_i^A (I guess they are very small)
- [ ] try use c again for transition into Erdos-Renyi
- [ ] profile time spend in type 1 / type 2
- [ ] clean up array view concept in WeightLayer
- [x] adapt exporter to directed saving of edges


# kd-tree

- must provide two operations:
    1. how many points of layer i are in a cell of level l<=i (in constant time with prefix sums)
    2. k-th point of layer i and cell of level l<=i (also constant with cool lookup array and prefix sums)

- cells have global index
    - ordering is by concatenating geometric ordering of all levels
        - because for a cell on level i all descendants on level j>i must be consecutive
        - examples ordering see below
    - index of first cell on level l is: ``((2^(d*l))-1)/((2^d)-1)``
    - children of cell i are: (2^d)*i+[1..2^d]
    - parent of cell i is: (i-1)/(2^d)
    - on level l are $(2^d)^l = 2^{d \cdot l}$ cells
    - each cell has 2^d children



# enumeration of cell-pairs and sampling of point-pairs (In Paper)

- for all layer pairs (i<=j)
    - **goal:** we want to sample all edges between points in layer i and j
    - search target level t determined by (w_i*w_j)/W
        - to enable all queries (level t <= layer i,j) must hold
    - for all cell pairs A,B that partition the space
        - pairs are of two types
            - type I: A,B are on level t and are touching
            - type II: A,B are not necessarily on level t but their parents touch and A,B do not touch
        - sample edges (according to type) between points in
            - cell A and layer i
            - cell B and layer j
        - if i=j only add edge (u,v) if u < v



# enumeration of cell-pairs and sampling of point-pairs (recursive idea)

notes:
- type 1 pairs include (A,A)
- in the original paper type 1 and type 2 pairs include A,B and B,A
- we use "visit" for cell-pairs and "sample" for point-pairs/layer-pairs


constraints:
1. all type 1 (touching and identical) cell pairs (of each level) must be visited. (show by induction over level)
    - observe: for a touching pair the parents must be identical or touching
    - for a pair (A,B) we make recursive call for
        - all children pairs in the cell if A=B
        - all children pairs (a,b) with a in A and b in B if A and B touch
2. all type 2 cell pairs must be visited.
    - should follow from first constraint and the fact that type 2 must have touching parents (since we sample all child pairs of touching parents)
3. no cell pair should be visited twice (again induction)
    - pairs with same parent
        - let c0,c1....ck be the children of the current cell
        - only make pairs ci,cj with 0<=i<=j<=k
    - pairs with different parent
        - for pair A,B we only make pairs (a,b) for a in A and b in B
        - we never have the pairs A,B and B,A like in the paper (except for A,A)
    - identical pairs (A,A)
        - we apply the edge-case to " remove all edges with u > v sampled in this iteration" if A=B and the layer are the same (i=j)
            - ATTENTION: for a cell pair we must sample all layer pairs and not only i<=j pairs
4. all cell pairs of the orig algo are visited once (omitting the double visit of the original algo)
    - follows from 1-3
5. for a cell pair A,B (of type 1 or 2), the recursive algo must sample the same pairs of weight-layers as the original algo
    - (in fact we must sample all layer-pairs that were sampled by cell pairs A,B and B,A combined)
    - type 1
        - original algo sampled only on one target level (for a pair of weigh-layers)
        - in recursive impl. the target level is determined like before and we try all possible layer-pairs that could have this level as target level
    - type 2
        - we sample all layer pairs that include this cell pair in their  partitioning data structure described in the paper (hopefully)



# global cell index examples

### dimension 1
- children=2*i+[1..2]
- parent=(i-1)/2
```
level: 4
cells: 15
l       cells   indexRange
0       1       [0..0]
1       2       [1..2]
2       4       [3..6]
3       8       [7..14]
```

### dimension 2
```
level: 4
cells: 85
l       cells   indexRange
0       1       [0..0]
1       4       [1..4]
2       16      [5..20]
3       64      [21..84]
```

### dimension 3
- children=8*i+[1..8]
- parent=(i-1)/8
```
level: 4
cells: 585
l       cells   indexRange
0       1       [0..0]
1       8       [1..8]
2       64      [9..72]
3       512     [73..584]
```

### general d
- children=(2^d)*i+[1..2^d]
- parent=(i-1)/(2^d)
```
2^0d    [0]
2^1d    [1..2^d]
2^2d    [2^d+1..2^d+2^2d]
2^3d    [2^d+2^2d+1..2^d+2^2d+2^3d]
```
