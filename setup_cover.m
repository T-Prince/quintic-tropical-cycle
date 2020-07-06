Z2 :=  FiniteField(2);
Z := Integers();
S7 := Sym(7);

//////////////////////////////////////////////////////////////////////
// the simplicial complex
//////////////////////////////////////////////////////////////////////

// There are 5 vertices, numbered 1..5
// There are 10 edges, numbered 1..10
// There are 10 triangles, numbered 1..10
// There are 5 tetrahedra, numbered 1..5

// faces_of_ts records the triangles that are the faces of each tetrahedron
faces_of_ts := [[2,4,8,9],[1,3,8,10],[3,4,5,6],[1,5,7,9],[2,6,7,10]];
// edges_of_ts records the edges of each tetrahedron
edges_of_ts := [[3,5,6,8,9,10],[1,2,4,5,9,10],[2,3,4,5,6,7],[1,2,3,7,8,9],[1,4,6,7,8,10]];
// vertices_of_ts records the vertices of each tetrahedron
vertices_of_ts := [[2,3,4,5],[1,3,4,5],[1,2,4,5],[1,2,3,5],[1,2,3,4]];

// edges_of_fs records the edges of each triangle
edges_of_fs := [[1,9,2],[8,10,6],[4,5,2],[6,5,3],[7,3,2],[7,6,4],[7,8,1],[10,5,9],[8,9,3],[1,10,4]];
// vertices_of_fs records the vertices of each triangle
vertices_of_fs := [[1,3,5],[2,3,4],[1,4,5],[2,4,5],[1,2,5],[1,2,4],[1,2,3],[3,4,5],[2,3,5],[1,3,4]];

// vertices_of_es records the vertices of each edge
vertices_of_es := [[1,3],[1,5],[2,5],[1,4],[4,5],[2,4],[1,2],[2,3],[3,5],[3,4]];

// ts_of_fs records the tetrahedra incident to each triangle
ts_of_fs := [[i : i in [1..5] | j in faces_of_ts[i]] : j in [1..10]];

// ts_of_fs records the tetrahedra incident to each edge
ts_of_es := [[i : i in [1..5] | j in edges_of_ts[i]] : j in [1..10]];

// es_of_vs records the edges incident to each vertex
es_of_vs := [[i : i in [1..10] | j in vertices_of_es[i]] : j in [1..5]];

// edges_tk records the edges in each face of the k-th tetrahedron
// each triple of edges is ordered so that the first edge
// has vertices v1 and v2, the second v2 and v3, and the third v1 and v3.
// The vertices v1, v2 and v3 are ordered by their index in [1..5].
edges_t1 := [[4,6,3],[3,2,1],[6,2,5],[4,5,1]];
edges_t2 := [[1,5,2],[3,4,2],[6,4,5],[1,6,3]];
edges_t3 := [[3,4,1],[5,4,2],[6,2,1],[6,5,3]];
edges_t4 := [[1,6,2],[4,3,2],[4,5,1],[5,6,3]];
edges_t5 := [[5,6,3],[4,3,2],[4,5,1],[1,6,2]];

// edge_lists is a list of the edges_tk
edge_lists := [edges_t1 ,edges_t2, edges_t3, edges_t4, edges_t5];

// v1_tk records the vertices in each face of the k-th tetrahedron.
v1_t1 := [[1,2,3],[1,3,4],[2,3,4],[1,2,4]];
v1_t2 := [[1,2,4],[1,3,4],[2,3,4],[1,2,3]];
v1_t3 := [[1,3,4],[2,3,4],[1,2,4],[1,2,3]];
v1_t4 := [[1,3,4],[1,2,4],[1,2,3],[2,3,4]];
v1_t5 := [[2,3,4],[1,2,4],[1,2,3],[1,3,4]];

// v_lists_1 is a list of the v1_tk
v_lists_1 := [v1_t1,v1_t2,v1_t3,v1_t4,v1_t5];

// v2_tk records the vertices in each edge of the k-th tetrahedron.
v2_t1 := [[1,4],[3,4],[1,3],[1,2],[2,4],[2,3]];
v2_t2 := [[1,2],[1,4],[1,3],[3,4],[2,4],[2,3]];
v2_t3 := [[1,4],[2,4],[1,3],[3,4],[2,3],[1,2]];
v2_t4 := [[1,3],[1,4],[2,4],[1,2],[2,3],[3,4]];
v2_t5 := [[1,3],[1,4],[2,4],[1,2],[2,3],[3,4]];

// v_lists_2 is a list of the v2_tk
v_lists_2 := [v2_t1,v2_t2,v2_t3,v2_t4,v2_t5];

// sanity check
&and[edge_lists[i] eq [[Index(edges_of_ts[i],e) : e in edges_of_fs[f]] : f in faces_of_ts[i]] : i in [1..5]];


// We store tuples of indices of faces which intersect.
// Note that the label 'be' denotes 'bad edge' and stores an edge which
// does not appear in an intersection. For example Meets_333 stores a vertex v
// together with a single edge meeting v excluded from the intersection.

Meets_1233 := [[t,f,be] : t in [1..5], f in [1..10], be in [1..10] | be in edges_of_fs[f] and f in faces_of_ts[t]];
Meets_1333 := [[t,v] : t in [1..5], v in [1..5] | v in vertices_of_ts[t]];
Meets_3333 := [1..5];

Meets_233 := [[f,be] : f in [1..10], be in [1..10] | be in edges_of_fs[f]];
Meets_123 := [[t,f,e] :  t in [1..5], f in [1..10], e in [1..10] | e in edges_of_fs[f] and f in faces_of_ts[t]];
Meets_133 := [[t,f,be] :  t in [1..5], f in [1..10], be in [1..10] | be in edges_of_fs[f] and f in faces_of_ts[t]];
Meets_333 := [[v,be] :  v in [1..5], be in [1..10] | v in vertices_of_es[be]];

Meets_12 := [[t,f] :  t in [1..5], f in [1..10] | f in faces_of_ts[t]];
Meets_13 := [[t,e] :  t in [1..5], e in [1..10] | e in edges_of_ts[t]];
Meets_23 := [[f,e] :  f in [1..10], e in [1..10] | e in edges_of_fs[f]];
Meets_33 := [[f,be] : f in [1..10], be in [1..10] | be in edges_of_fs[f]];

//////////////////////////////////////////////////////////////////////
// the open cover
//////////////////////////////////////////////////////////////////////

// On each face there are 25 negative vertices.  We subdivide each face into hexagons, called blobs, and pentagons.
// blob_i records the negative vertices, in counterclockwise order, met by each blob.
blob_1 := [2,3,4,12,11,10];
blob_2 := [4,5,6,14,13,12];
blob_3 := [6,7,8,16,15,14];
blob_4 := [11,12,13,19,18,17];
blob_5 := [13,14,15,21,20,19];
blob_6 := [18,19,20,24,23,22];

// bn_matrix is a list of the blob_i
bn_matrix := [blob_1,blob_2,blob_3,blob_4,blob_5,blob_6];

// There is one open set in the cover for each tetrahedron, one for each blob, one for each negative vertex, one
// for each positive vertex, and one for each vertex of the simplicial complex.

// Over an open set of type t ("tetrahedron") the sheaf has rank 7.
// Over an open set of type b ("blob") it has rank 7.
// Over an open set of type n it has rank 4.
// Over an open set of type p it has rank 4.
// Over an open set of type v ("vertex") it has rank 7.

// There are no intersections of type tt.
// Over a tb intersection the sheaf has rank 7.
// Over a tn intersection it has rank 7.
// Over a tp intersection it has rank 7.
// Over a tv intersection it has rank 7.

// There are no intersections of type bb.
// Over a bn intersection the sheaf has rank 7.
// There are no intersections of type bp.
// There are no intersections of type bv.

// Over an nn intersection the sheaf has rank 5.
// Over an np intersection it has rank 5 or 7, depending on the intersection; see singular_23 below.
// Over an nv intersection it has rank 7.

// Over a pp intersection the sheaf has rank 7.
// Over a pv intersection the sheaf has rank 7.

// There are no vv intersections.

// The sheaf has rank 7 over all triple intersections.  The non-empty triple intersections are of type tbn, tnn, tnp, tnv, tpp, tpv, bnn, nnp, npp, npv.
// It has rank 7 over all quadruple intersections. The non-empty quadruple intersections are of type tbnn, tnnp, tnpp, tnpv.
// There are no quintuple intersections.

// We refer to the t type open sets as class 1, the n and b as class 2, the p type open sets as class 3, and v type open sets as class 4.
// We extend this numbering to intersections, thus a tnn intersection is of class 122, and a tnpv intersection of class 1234.

// Place a small disc about each negative vertex.  The discs for adjacent negative vertices intersect, also in discs.
// We enumerate these intersections. The k-th entry in nn_matrix records the indices of the negative-vertex-discs
// that meet in the k-th intersection.  
nn_matrix :=
[
[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[2,10],[4,12],[6,14],[8,16],[10,11],[11,12],[12,13],
[13,14],[14,15],[15,16],[11,17],[13,19],[15,21],[17,18],[18,19],[19,20],[20,21],[18,22],[20,24],[22,23],[23,24],[23,25]
];

// Each blob determines six nn-intersections, in cyclic order -- those given by the edges of the hexagon that defines
// the blob.  bnn_k gives the indices in nn_matrix of the edges of the k-th blob
bnn_1 := [9,2,3,10,14,13];
bnn_2 := [10,4,5,11,16,15];
bnn_3 := [11,6,7,12,18,17];
bnn_4 := [19,14,15,20,23,22];
bnn_5 := [20,16,17,21,25,24];
bnn_6 := [26,23,24,27,29,28];

// bnn_matrix is a list of the bnn_k
bnn_matrix := [bnn_1,bnn_2,bnn_3,bnn_4,bnn_5,bnn_6];

// sanity check
for i in [1..6] do
    blob := bn_matrix[i];
    edges := [[blob[6],blob[1]]] cat [[blob[i],blob[i+1]] : i in [1..5]];
    bnn_matrix[i] eq [Index(nn_matrix, Sort(e)) : e in edges];
end for;


// Each triangle contains 15 positive vertices, 5 on each edge.  There is an open ball in the cover for each positive vertex.
// There are 24 negative-negative-positive (nnp) triple intersections, eight along each edge.  Each nn intersection along
// an edge meets precisely one positive-vertex ball; hence the 24 nnp intersections are uniquely determined by the
// nn-intersections along the edges of the triangle.  In counter-clockwise order, these are given by nnp_lists.
nnp_lists := [[1,2,3,4,5,6,7,8],[8,12,18,21,25,27,29,30],[1,9,13,19,22,26,28,30]];

// There are three possible directions that the line segment determining (and determined by) an nn-intersection can point:
// (1) N-S; (2) SW-NE; (3) SE-NW.  nn_direction records these directions.
nn_direction := [2,3,2,3,2,3,2,3,1,1,1,1,2,3,2,3,2,3,1,1,1,2,3,2,3,1,1,2,3,1];

// We now fix trivializations of the sheaf over the 5 open sets of type v and the 5 open sets of type t.
// Paste[[i,j]] records the identification between the trivialization over the i-th open set of type t
// and the j-th open set of type v.

Paste := AssociativeArray();


// Without loss of generality we can insist that everything glues to the 5-th vertex via the identity map
for t in [1..4] do
    Paste[[t,5]] := (S7 ! 1);
end for;

// Without loss of generality we may assume that tetrahedron 1 glues to all vertices via the identity map,
// and that tetrahedron 3 glues to vertex 1 by the identity map, and that tetrahedron 5 glues to vertex 1
// via the identity map.  (Here we have chosen a maximal tree in the incidence graph and set the
// corresponding gluing maps to be the identity.)

// The other gluing maps then follow, by considering the standard model of the affine manifold structure
// on tetrahedron 4 and then doing a lot of monodromy calculations.
Paste[[2,1]] := S7 ! (1,2)(4,5);
Paste[[3,1]] := S7 ! 1;
Paste[[4,1]] := S7 ! (1,2)(6,7);

Paste[[1,2]] := S7 ! 1;
Paste[[3,2]] := S7 ! (1,6)(3,4);
Paste[[4,2]] := S7 ! (3,4)(2,7);

Paste[[1,3]] := S7 ! 1;
Paste[[2,3]] := S7 ! (2,3)(5,6);
Paste[[4,3]] := S7 ! (2,3)(4,7);

Paste[[1,4]] := S7 ! 1;
Paste[[2,4]] := S7 ! (2,5)(3,6);
Paste[[3,4]] := S7 ! (1,4)(3,6);


Paste[[5,1]] := S7 ! 1;
Paste[[5,2]] := Paste[[4,2]]*Paste[[4,1]]^-1*(S7 ! (2,6)(3,5));
Paste[[5,3]] := Paste[[4,3]]*Paste[[4,1]]^-1*(S7 ! (2,3)(5,6));
Paste[[5,4]] := Paste[[2,4]]*Paste[[2,1]]^-1*(S7 ! (3,7)(2,4));

// tfe_monos[[i,j,k]] records the monodromy transformation given by transport, with respect to
// the trivialization in the i-th open set of type t, around the five 'legs' specified by edge
// k in triangle j.
tfe_monos := AssociativeArray();
for i in [1..5] do
    for j in faces_of_ts[i] do
	for k in edges_of_fs[j] do
	    other_t := [l : l in [1..5] | j in faces_of_ts[l] and l ne i][1];
	    v1, v2 := Explode(vertices_of_es[k]);
	    tfe_monos[[i,j,k]] := Paste[[i,v2]]^-1*
				  Paste[[other_t, v2]]*
				  Paste[[other_t,v1]]^-1*
				  Paste[[i,v1]];
	end for;
    end for;
end for;

MPaste_Inv := AssociativeArray();
MPaste := AssociativeArray();
for A in Keys(Paste) do
	MPaste_Inv[A] := SparseMatrix(PermutationMatrix(Z2,Paste[A]^-1));
	MPaste[A] := SparseMatrix(PermutationMatrix(Z2,Paste[A]));
end for;

fe_monos := AssociativeArray();
for i in [1..10] do
	for j in edges_of_fs[i] do
		tetra := ts_of_fs[i];
		fe_monos[[i,j]] := Paste[[tetra[2],vertices_of_es[j][1]]]*
			Paste[[tetra[2],vertices_of_es[j][2]]]^-1*
				Paste[[tetra[1],vertices_of_es[j][2]]]*
					Paste[[tetra[1],vertices_of_es[j][1]]]^-1;
	end for;
end for;

blob_paste := AssociativeArray();
for i in [1..5] do
    for j in faces_of_ts[i] do
	// For each face, glue each blob on that face to precisely one of the 
	// tetrahedra incident to that face using the identity map.
	if i eq ts_of_fs[j][1] then
	    for k in [1..6] do
		blob_paste[[i,j,k]] := (S7 ! 1);
	    end for;
	else
	    // for each of the remaining tetrahedra the paste map is determined by the monodromy around a branch of the discriminant locus.
		// In fact, this matrix is such that, after pasting to a neighbouring vertex (with 'Paste') and then back to the first tetrahedron
		// (using the inverse'Paste') we recover precise a monodomry matrix.
	    blob_paste[[i,j,1]] := Paste[[ts_of_fs[j][1],vertices_of_fs[j][3]]]^-1*Paste[[i,vertices_of_fs[j][3]]]*tfe_monos[[i,j,edges_of_fs[j][1]]];
	    blob_paste[[i,j,2]] := Paste[[ts_of_fs[j][1],vertices_of_fs[j][1]]]^-1*Paste[[i,vertices_of_fs[j][1]]]*tfe_monos[[i,j,edges_of_fs[j][3]]];
	    blob_paste[[i,j,3]] := Paste[[ts_of_fs[j][1],vertices_of_fs[j][2]]]^-1*Paste[[i,vertices_of_fs[j][2]]]*tfe_monos[[i,j,edges_of_fs[j][3]]];
	    blob_paste[[i,j,4]] := Paste[[ts_of_fs[j][1],vertices_of_fs[j][1]]]^-1*Paste[[i,vertices_of_fs[j][1]]]*tfe_monos[[i,j,edges_of_fs[j][1]]];
	    blob_paste[[i,j,5]] := Paste[[ts_of_fs[j][1],vertices_of_fs[j][2]]]^-1*Paste[[i,vertices_of_fs[j][2]]]*tfe_monos[[i,j,edges_of_fs[j][1]]];
	    blob_paste[[i,j,6]] := Paste[[ts_of_fs[j][1],vertices_of_fs[j][3]]]^-1*Paste[[i,vertices_of_fs[j][3]]]*tfe_monos[[i,j,edges_of_fs[j][1]]];
	end if;
    end for;
end for;

// We fix a trivialization over each of the segments of the singular locus and record the paste maps to neighbouring tetrahedra
// in the array 'sing_paste'.
sing_paste := AssociativeArray();

// We form a projection matrix by looking at the orbits of the monodromy matrix.
// There is parity issue: we cannot simulatneously fix a single trivialisation for 
// every segment in a two dimensional face with the same direction, the ordering 
// of the final two columns has to be recorded seperately.

for j in [1..10] do
	i := ts_of_fs[j][1];
	for k in [1..3] do
			X := Orbits(sub<S7 | tfe_monos[[i,j,edges_of_fs[j][k]]]> );
			entries := [<a,b,1> : a in [1..7], b in [1..5] | a in X[b]];
			sing_paste[[i,j,k,1]] := SparseMatrix(Z2,7,5,entries);
			
			// Monodromy around loops based at a point in the discriminant locus can interchange 
			// columns recording which pairs of sheets combine over the given segement of the discriminant locus.
			// Therefore we record each possibility using a fourth index in the associate array.
			non_trivial_pair := [b : b in [1..5] | #[c : c in [1..7] | <c,b,1> in entries] eq 2];
			sing_paste[[i,j,k,2]] := sing_paste[[i,j,k,1]]*SparseMatrix(PermutationMatrix(Z2,Sym(5)!(non_trivial_pair[1],non_trivial_pair[2])));
	end for;
end for;

// Restriction maps to the remaining tetrahedra can be related to those already defined using paste maps between open sets
// corresponding to tetrahedra and vertices.
for j in [1..10] do
	i := ts_of_fs[j][2];
	iprime := ts_of_fs[j][1];
		sing_paste[[i,j,1,1]] := MPaste_Inv[[i,vertices_of_fs[j][1]]]*MPaste[[iprime,vertices_of_fs[j][1]]]*sing_paste[[iprime,j,1,1]];
		sing_paste[[i,j,2,1]] := MPaste_Inv[[i,vertices_of_fs[j][2]]]*MPaste[[iprime,vertices_of_fs[j][2]]]*sing_paste[[iprime,j,2,1]];
		sing_paste[[i,j,3,1]] := MPaste_Inv[[i,vertices_of_fs[j][3]]]*MPaste[[iprime,vertices_of_fs[j][3]]]*sing_paste[[iprime,j,3,1]];
		sing_paste[[i,j,1,2]] := MPaste_Inv[[i,vertices_of_fs[j][3]]]*MPaste[[iprime,vertices_of_fs[j][3]]]*sing_paste[[iprime,j,1,2]];
		sing_paste[[i,j,2,2]] := MPaste_Inv[[i,vertices_of_fs[j][1]]]*MPaste[[iprime,vertices_of_fs[j][1]]]*sing_paste[[iprime,j,2,2]];
		sing_paste[[i,j,3,2]] := MPaste_Inv[[i,vertices_of_fs[j][2]]]*MPaste[[iprime,vertices_of_fs[j][2]]]*sing_paste[[iprime,j,3,2]];
end for;

// Indexing the segments of the singular locus which appear in a two dimensional face
// we record which ordering of the final two columns in the corrsponding sing_paste
// matrix are required.
nn_parity := [1,2,2,1,1,2,2,1,2,2,2,2,2,2,1,1,2,2,1,1,1,1,2,2,1,2,2,2,2,1];

// There are 13 matrices corresponding to intersections between class 2 and class 3 open sets.
// The array singular_23 records which of these meet the discriminant locus.
singular_23 := [1,4,7,10,13];

// Consider the triangular faces of class 23 (np) in the acyclic cover which do not meet the 
// singular locus. Trivialisations of \pi_\star\ZZ_2 can be identified with the trivialisation
// over one of the vertices of the edge associated with the class 3 open set in the intersection.
// If this is the vertex has the lower index of the two vertices of this edge, we record a zero, if
// not we record a 1.
offset_23 := [0,1,1,0,0,0,0,1,1,0,0,0,0];
offset_22 := [1,1,0,0,1,1,0,0];

// The vertices v1, v2 and v3 of the given face are numbered anti-clockwise,
// in ascending order of their index in {1,...,5}.
// Edges are ordered as follows: edge 1 contains v1 and v2;
// edge 2 contains v2 and v3;
// edge 3 contains v1 and v3.
// Given vertex i of edge j of a face v_lookup[(j-1)*2+i] is the corresponding index of the vertex of the face.
// Note that we always orient our edges from lowest index vertex to highest and order.
v_lookup := [1,2,2,3,1,3];

n_sing := SparseMatrix(Z2,5,4,[<1,1,1>,<2,2,1>,<3,3,1>,<4,4,1>,<5,4,1>]);

// We describe the restriction map from an open set containing a positive vertex
// to an open subset contained inside the interior of a tetrahedron.

p_sing := AssociativeArray();

for i in [1..10] do
	t1, t2, t3 := Explode(ts_of_es[i]);
	v1, v2 := Explode(vertices_of_es[i]);
	faces := [f : f in faces_of_ts[t1] | i in edges_of_fs[f]];
	X := Orbits(sub<S7 | tfe_monos[[t1,faces[1],i]],tfe_monos[[t1,faces[2],i]] > );
	entries := [<a,b,1> : a in [1..7], b in [1..4] | a in X[b]];
	p_sing[[t1,i]] := SparseMatrix(Z2,7,4,entries);
	p_sing[[t2,i]] := MPaste_Inv[[t2,v1]]*MPaste[[t1,v1]]*p_sing[[t1,i]];
	p_sing[[t3,i]] := MPaste_Inv[[t3,v1]]*MPaste[[t1,v1]]*p_sing[[t1,i]];
end for;

// The 'beef_up' function will replace the entries '1' of a matrix with seven by seven identity blocks.
// This will often be used to pass from 'incidence' matrices, recording the intersections of the open
// sets to Cech differentials where the transition function from one open set to the next is the identity.
beef_up := function(A)
	beef_A := SparseMatrix(Z2,0,NumberOfColumns(A)*7,[]);
	beef_row := SparseMatrix(Z2,7,0,[]);
	for i in [1..NumberOfRows(A)] do
		for j in [1..NumberOfColumns(A)] do
			if A[i,j] eq 1 then
				beef_row := HorizontalJoin(beef_row,IdentitySparseMatrix(Z2,7));
			else
				beef_row := HorizontalJoin(beef_row,SparseMatrix(Z2,7,7,[]));			
			end if;
		end for;
		beef_A := VerticalJoin(beef_A,beef_row);
		beef_row := SparseMatrix(Z2,7,0,[]);
	end for;
	
	return beef_A;
end function;

// oneify replaces every non-zero entry in a matrix
// with the value 1, and converts it to a matrix over Z2.
oneify := function(A)
	for i in [1..NumberOfRows(A)] do
		for j in [1..NumberOfColumns(A)] do
			if A[i,j] ne 0 then
				A[i,j] := 1;
			end if;
		end for;
	end for;
	return ChangeRing(A,Z2);
end function;