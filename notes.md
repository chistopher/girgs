# Naming

layer:
    - partitioning of points in weight-layers

level:
    - depth in kd-tree


# improvement ideas
- use level local index for cells
- in WeightLayer: use array_view concept to provide access to all nodes in ViA
    - e.g. only one method that returns pointer to first element and size of collection
    - removes need for multiple calls to k-th point
- use omp to process recursive calls in parallel (somehow aggregate on edge vectors?)
- stop building empty weight_layers
- set levels to (L-2)/D+1 ? because we need no coord_helper if we stop recursion after (L-2)/D
- make graph member to reduce stack size ? needs benchmarking
- move index helper functions to common.h
- pre-compute expressions used in for loop condition



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
6. all point-pairs should be sampled once
    - follows from constraints 4,5 and the fact that original algo samples all point-pairs
    - explanation for type 2 is in comment after "sample edges between V_i^A and V_j^B"


Alternative proof for 4-6
we sample all point-pairs once
    let a,b be two points from layer i and j respectively

    we now show that the pair a,b is sampled exactly once.

    let k be the deepest (with highest index) level in which their cells touch
    -> their cells touch in level k, thus they also touch in all levels <=k

    observation 1: there are type 1 cell-pairs containing our points on all levels <=k
    observation 2: there is only one type 2 cell-pair containing our points and it is on level k+1
        (unless the points are in touching cells on the deepest level; for bounded number of layers we change i+j to min(maxLayer, i+j)) TODO verify that remark 1 still holds

    1) (i+j<=k): k is deeper than (or equal to) i+j
        1.1) we sample them in exactly one type 1 cell-pair:
            since i+j<=k, observation 1 states that there is a type 1 cell-pair on level i+j
            this cell-pair will sample our point-pair
            the point-pair is not sampled by other type 1 cell-pairs, because type 1 cell-pairs on level l sample only layers where i+j=l
            -> we sample the point-pair only on level i+j
        1.2) we do not sample them in any type 2 cell-pair:
            (observation 2) the only type 2 cell-pair containing our points is on level k+1
            for this type 2 cell-pair we sample all layer-pairs i',j' >= k+1 > k => i+j
            inequality yields: i',j'>i,j
            -> the layers i,j are not sampled by that (or any other) type 2 cell-pair
    2) (k<i+j): i+j is deeper than k
        2.1) we do not sample them in any type 1 cell-pair:
            a layer-pair can only be sampled by a type 1 cell-pair when the cell-pair is on level i+j
            all type 1 cell-pairs from observation 1 are on levels <=k<i+j (since the cells containing a,b do not touch in level i+j)
            -> they are not sampled in any type 1 cell-pair
        2.2) we sample them in exactly one type 2 cell-pair:
            (observation 2) the only type 2 cell-pair containing our points is on level k+1
            for this type 2 cell-pair we sample all layer-pairs i',j' >= k+1 > k
            since i,j > k, our two points are one of those pairs


    Remark: can we change the level l where layer-pairs i,j are sampled by type 1 cell-pairs (currently l=i+j)?
    - proof holds if l satisfies:
        - unique for i,j    (needed for 1.1)
        - l>i and l>j       (needed for for inequality in 1.2)
    - max(i,j) is ok
    Remark 2: can be verified empirically by changing the edge probability to 1.0 and expecting the complete graph




```
sample(cell A, cell B)
    let l be level of cells

    if(touching A,B or A=B)
        // sample all type 1 occurrences with this cell pair
        sample edges for weightlayers i,j with (wi*wj/W) somehow determines l as target level // is i+j = l right ??? -> see "we sample all point-pairs" proof
    else // not touching
        // sample all type 2 occurrences with this cell pair
        // the type 2 occurrences can come from any pair of lower levels i,j>=l
        // assert(parent must be touching) to be sure we are type II pair; otherwise we would not have done this recursive call
        for all weightlayer pairs i,j>=l
            sample edges between V_i^A and V_j^B
            // for two points x in V_i^A and y in V_j^B this can happen only one time, because there is just one pair of cells (A', B') in the hierarchy with x in A' and y in B' where:
            // 1. A' and B' are not touching
            // 2: parent(A') touches parent(B')
            // It is easy to see, that no children of these cells will sample those points again, because the recursion ends here

    if(A=B)
        recursive call for all children pairs of cell A=B

    if(touching A,B and A!=B)
        recursive call for all children pairs (a,b) where a in A and b in B
        // these will be type 1 if a and b touch or type 2 if they dont
```



# global cell index examples

### dimension 1
- children=2*i+[1..2]
- parent=(i-1)/2
```
level: 4
cells: 85
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




# kromer

entry:
    HyperbolicLinear::linearSampling
        - samples positions
        - calls sample edges

    HyperbolicLinear::sampleEdges
        - constructs weight layers (vector<set<Point>>) to be held by HyperbolicLinear
        - constructs GIRGS data structure (does not use its member GeometricDS)
            - for layer 0 to L-1 construct a layer of cells with PointsFromCells constructor

DEFPRECATED:
TreePoints constructs from set of points
    has vector<vector<set<Point>>>: layer cell points_set

GeometricDS<TreePoints>
    init:
        computes weight layers
    has vector<TreePoints>
