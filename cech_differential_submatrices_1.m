////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// Here we construct the (non-zero) submatrices of the Cech differential between
// double and triple intersections. These rely on the arrays and matrices defined
// in setup_cover.m. We construct two sets of matrices, incidence_ab_cde are 
// submatrices of the differential for the constant sheaf Z_2 on S^3 
// (and can then be used as a sanity check). The matrices full_ab_cde are 
// the submatrices of the cech differential of \pi_\star Z_2. Submatrices are 
// constructed row by row. We recall that sets in the given open cover of S^3,
// acyclic with respect to \pi_\star Z_2 are divided into 4 classes:
//
// Class 1: Interiors of tetrahedra.
// Class 2: Products of a small open interval with a open set contained in a 2-dimensional face.
// These form two subclasses, 'n' and 'b'. Open sets of class 'n' restrict to neightbourhoods
// of a negative vertex, while those of class 'b' correspond to regions enclosed by \Delta.
// Class 3: These are neighbourhoods of positive vertices (contained in edges).
// Class 4: These are nieghbourhoods of one of the 5 vertices of the boundary.
//
// Submatrices are labelled ab_cde for a,b,c,d, and e in {1,2,3,4}. These are 
// zero unless {a,b} is contained in {c,d,e}.
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
// Incidence_12_122
////////////////////////////////////////////////////////////////////////////////////

// The incidence matrix 12_122 for a fixed maximal cell is divided into four cases
// depending on the subclass of the sets of class 2 involved. Note that no intersection
// with a 'b' class set will produce a triple intersection involving a set of class 'nn'.
// The three remaining cases we label A (b sets to bn sets) , B: (n to bn) and C: (n sets to nn sets).

mtx_A := SparseMatrix(Z2,0,0,[]);
mtx_B := SparseMatrix(Z2,0,0,[]);
mtx_C := SparseMatrix(Z2,0,0,[]);

// Build the block diagonal matrix for A.
temp_A := SparseMatrix(Z2,6,1,[<1,1,1>,<2,1,1>,<3,1,1>,<4,1,1>,<5,1,1>,<6,1,1>]);
for i in [1..6] do
mtx_A := DiagonalJoin(mtx_A,temp_A);
end for;

// Matrix B: look up the negative nodes meeting the desired bn type intersection.
mtx_B := SparseMatrix(Z2,36,25,[]);
for i in [1..6] do
	for j in [1..6] do
		mtx_B[6*(i-1)+j,bn_matrix[i][j]] := 1;
	end for;
end for;

// Matrix C: look up the pair of negative nodes meeting the desired 'nn' type intersection.
mtx_C := SparseMatrix(Z2,30,25,[]);
for i in [1..30] do
	mtx_C[i,nn_matrix[i][1]] := 1;
	mtx_C[i,nn_matrix[i][2]] := 1;
end for;

// We join A, B and C and copy 20 times (once for each of the 5 tetraedra and 10 four faces of each tetrahedra).
temp := VerticalJoin(HorizontalJoin(mtx_A,mtx_B),HorizontalJoin(SparseMatrix(Z2,30,6,[]),mtx_C));
incidence_12_122 := SparseMatrix(Z2,0,0,[]);
for i in [1..20] do
	incidence_12_122 := DiagonalJoin(incidence_12_122,temp);
end for;

// Since trivializations over sets of class 12 and 122 agree (induced by the open set of class 1), we can
// simply apply the function 'beef_up'. Note this will always be the case for submatrices 'full_1a_1bc'.
full_12_122 := beef_up(incidence_12_122);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_22_122
////////////////////////////////////////////////////////////////////////////////////

// Recall that sets in class 22 are of the form 'nn' ('n' type meets 'n' type),
// or 'bn' ('b' type meets 'n' type). Note the case 'bb' does not appear.

incidence_22_122 := SparseMatrix(Z2,0,660,[]);
row := SparseMatrix(Z2,66,0,[]);
full_22_122 := SparseMatrix(Z2,0,4020,[]);
rf := SparseMatrix(Z2,66*7,0,[]);

for i in [1..5] do
	for j in [1..4] do
		for k in [1..10] do
			if k eq faces_of_ts[i][j] then

				// Divide the cases into the 'bn' and 'nn' type 22 intersections.
				block_bn := SparseMatrix(Z2,0,0,[]);
				block_nn := SparseMatrix(Z2,0,0,[]);

				// Iterate over the 36 'bn' intersections (first fixing the 'b' set and iterating over the 6 neighbouring 'n' sets)
				for l in [1..36] do

					// Sections over 'bn' intersections are trivialized with the trivialization from the 'b' open set.
					// Sections over a 123 intersection are trivialized with the trivialization from the set of class 1.
					// The transition function is recorded in 'blob_paste'. Note (l-1) div 6 + 1 recovers the region
					// (between 1 and 6) determined by l.

					block_bn := DiagonalJoin(block_bn, SparseMatrix(PermutationMatrix(Z2,blob_paste[[i,k,(l-1) div 6 + 1]]^-1)));

				end for;
				
				for l in [1..30] do
				
					// Sections over 'nn' intersections carry their own trivialization, and sing_paste
					// records the transition map to the 123 trivialization. Note nn_parity is needed here
					// as the transition map depends not only on the direction of the segement of \Delta
					// contained in the 'nn' intersection, but the location of that segment in the face.

					block_nn := DiagonalJoin(block_nn, sing_paste[[i,k,nn_direction[l],nn_parity[l]]]);
					
				end for;
				
				// Join the 'bn' and 'nn' type blocks.
				block := DiagonalJoin(block_bn,block_nn);

				// 'row' forms part of a row of the matrix 'incidence' (which is insensitive to the various trivializations).
				// 'rf' (row_full) forms part of the matrix incidnece_full.				
				row := HorizontalJoin(row,IdentitySparseMatrix(Z2,66));
				rf := HorizontalJoin(rf,block);
			else
				row := HorizontalJoin(row,SparseMatrix(Z2,66,66,[]));
				rf := HorizontalJoin(rf,SparseMatrix(Z2,66*7,30*5+36*7,[]));
			end if;
		end for;

		incidence_22_122 := VerticalJoin(incidence_22_122,row);
		full_22_122 := VerticalJoin(full_22_122,rf);
		row := SparseMatrix(Z2,66,0,[]);
		rf := SparseMatrix(Z2,66*7,0,[]);
	end for;
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_12_123
////////////////////////////////////////////////////////////////////////////////////

// The following three matrices 'block_i' record which of the 25 'n' type open sets are involved in each of the 13 intersections of type '23', running
// along the each edge of a two dimensional face.
block_1 := SparseMatrix(Z2,13,25,[<1,1,1>,<2,2,1>,<3,2,1>,<4,3,1>,<5,4,1>,<6,4,1>,<7,5,1>,<8,6,1>,<9,6,1>,<10,7,1>,<11,8,1>,<12,8,1>,<13,9,1>]);
block_2 := SparseMatrix(Z2,13,25,[<1,9,1>,<2,8,1>,<3,8,1>,<4,16,1>,<5,15,1>,<6,15,1>,<7,21,1>,<8,20,1>,<9,20,1>,<10,24,1>,<11,23,1>,<12,23,1>,<13,25,1>]);
block_3 := SparseMatrix(Z2,13,25,[<1,1,1>,<2,2,1>,<3,2,1>,<4,10,1>,<5,11,1>,<6,11,1>,<7,17,1>,<8,18,1>,<9,18,1>,<10,22,1>,<11,23,1>,<12,23,1>,<13,25,1>]);
block := VerticalJoin(VerticalJoin(block_1,block_2),block_3);
block := HorizontalJoin(SparseMatrix(Z2,39,6,[]),block);

incidence_12_123 := SparseMatrix(Z2,0,0,[]);

for i in [1..20] do
	incidence_12_123 := DiagonalJoin(incidence_12_123,block);
end for;
full_123 := beef_up(incidence_12_123);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_23_123
////////////////////////////////////////////////////////////////////////////////////

incidence_23_123 := SparseMatrix(Z2,0,390,[]);
full_23_123 := SparseMatrix(Z2,0,2430,[]);
row := SparseMatrix(Z2,39,0,[]);
rf := SparseMatrix(Z2,39*7,0,[]);

// Recall that there are 13 intersections of the form '23' along each edge,
// of these 13 a total of 5 contain a segment of the discriminant locus,
// while the remaining 8 do not.

for i in [1..5] do
	for j in [1..4] do
		for k in [1..10] do
			if k eq faces_of_ts[i][j] then
				row := HorizontalJoin(row,IdentitySparseMatrix(Z2,39));
				
				block := SparseMatrix(Z2,0,0,[]);
				// Running along an edge of a triangular face we encounter 13 faces of class '23',
				// the index of which is stored in 'pos'.
				for l in [1..39] do
						pos := (l-1) mod 13 + 1;
						edge_index := (l-1) div 13 +1;
						
						// If the current face contains a piece of discrimiant locus,
						// we insert a submatrix from 'sing_paste'.  
						if pos in singular_23 then
							block := DiagonalJoin(block, sing_paste[[i,k,edge_index,1]]);
							
						// Otherwise we compute the map relating the between the trivialisation of \pi_\star\ZZ_2
						// over the given '23' face with the tetrahedron appearing in the class '123' triple intersection.
						// The trivialisation over a type '23' face is made by recalling that the sheets of \pi over this face 
						// can be identified with those over a vertex. The index of this vertex is computed using v_lookup.
						// The array offset_23 indicates which vertex of the given edge we use.
						else
							v_index := v_lookup[(edge_index-1)*2+offset_23[pos]+1];							
							block := DiagonalJoin(block, MPaste_Inv[[i,vertices_of_fs[k][v_index]]]);
						end if;
				end for;
				rf := HorizontalJoin(rf,block);
			else
				row := HorizontalJoin(row,SparseMatrix(Z2,39,39,[]));
				rf := HorizontalJoin(rf,SparseMatrix(Z2,39*7,243,[]));
			end if;
		end for;
		incidence_23_123 := VerticalJoin(incidence_23_123,row);
		row := SparseMatrix(Z2,39,0,[]);

		full_23_123 := VerticalJoin(full_23_123,rf);
		rf := SparseMatrix(Z2,39*7,0,[]);
	end for;
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_13_123
////////////////////////////////////////////////////////////////////////////////////

incidence_13_123 := SparseMatrix(Z2,0,0,[]);

// For each of the 13 '23' intersections we record which
// type '3' set appears in the intersection.
block := SparseMatrix(Z2,13,5,[<1,1,1>,<2,1,1>,<3,2,1>,<4,2,1>,<5,2,1>,<6,3,1>,<7,3,1>,<8,3,1>,<9,4,1>,<10,4,1>,<11,4,1>,<12,5,1>,<13,5,1>]);
big_block := SparseMatrix(Z2,0,30,[]);
row := SparseMatrix(Z2,13,0,[]);
for h in [1..5] do
	for i in [1..4] do
		for j in [1..3] do
			for k in [1..6] do
				if k eq edge_lists[h][i][j] then
					row := HorizontalJoin(row,block);
				else
					row := HorizontalJoin(row,SparseMatrix(Z2,13,5,[]));
				end if;
			end for;
			big_block := VerticalJoin(big_block,row);
			row := SparseMatrix(Z2,13,0,[]);
		end for;
	end for;
	incidence_13_123 := DiagonalJoin(incidence_13_123,big_block);
	big_block := SparseMatrix(Z2,0,30,[]);
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_12_124
////////////////////////////////////////////////////////////////////////////////////

incidence_12_124 := SparseMatrix(Z2,0,0,[]);

// This block matrix records which of the 25 open sets of type 'n' in class 2 meet
// an open set of class 4. Horizontal join an empty matrix corresponding to the
// type 'b' sets in class 2.

block := SparseMatrix(Z2,3,25,[<1,1,1>,<2,9,1>,<3,25,1>]);
block := HorizontalJoin(SparseMatrix(Z2,3,6,[]),block);

for i in [1..20] do
	incidence_12_124 := DiagonalJoin(incidence_12_124,block);
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_24_124
////////////////////////////////////////////////////////////////////////////////////

// Recall that faces of class 24 always use the trivialisation fixed for 
// the corresponding vertex, while class 124 uses the trivialisation fixed over
// the corresponding tetrahedron.

incidence_24_124 := SparseMatrix(Z2,0,30,[]);
full_24_124 := SparseMatrix(Z2,0,30*7,[]);

// 'row' records part of the 'incidence' matrix.
// 'rf' records part of the 'incidence_full' matrix.
row := SparseMatrix(Z2,3,0,[]);
rf := SparseMatrix(Z2,21,0,[]);

// Iterate over [i,j] pairs where i runs over tetrahedra.
// and j runs over faces of the ith tetrahedron. Then iterate k
// over the 2-faces, and check whether this is the 
// face fixed by i and j.
for i in [1..5] do
	for j in [1..4] do
		for k in [1..10] do

			// If k does index the face determined by [i,j],
			// append a 3x3 identiy block, representing the
			// three vertices of the kth face.
			if k eq faces_of_ts[i][j] then
				row := HorizontalJoin(row,IdentitySparseMatrix(Z2,3));

				// Since the class '24' open sets are trivialized using the trivialization
				// inherited from the class '4' (neighbourhood of a vertex) open set
				// we need to paste from this trivialization to that over the the class '1' open set.
				// Since Paste stores maps in the other direction (class '1' to '4') we use MPaste_Inv.
				v_pastes := [MPaste_Inv[[i,vertices_of_fs[k][l]]] : l in [1..3]];
				block := DiagonalJoin(DiagonalJoin(v_pastes[1],v_pastes[2]),v_pastes[3]);
				rf := HorizontalJoin(rf,block);
			else
				row := HorizontalJoin(row,SparseMatrix(Z2,3,3,[]));
				rf := HorizontalJoin(rf,SparseMatrix(Z2,21,21,[]));
			end if;
		end for;
		incidence_24_124 := VerticalJoin(incidence_24_124,row);
		full_24_124 := VerticalJoin(full_24_124,rf);
		row := SparseMatrix(Z2,3,0,[]);
		rf := SparseMatrix(Z2,21,0,[]);
	end for;
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_14_124
////////////////////////////////////////////////////////////////////////////////////

block := SparseMatrix(Z2,12,4,[]);
incidence_14_124 := SparseMatrix(Z2,0,0,[]);

// Since both class 14 and 124 open sets have trivializations inherited from the same
// (class 1) open set we use the 'beef_up' function to obtain the submatrix of incidence_full.
for h in [1..5] do
	for i in [1..4] do
		block[3*i-2,v_lists_1[h][i][1]] := 1;
		block[3*i-1,v_lists_1[h][i][2]] := 1;
		block[3*i,v_lists_1[h][i][3]] := 1;
	end for;
	incidence_14_124 := DiagonalJoin(incidence_14_124,block);
	block := SparseMatrix(Z2,12,4,[]);
end for;


////////////////////////////////////////////////////////////////////////////////////
// Incidence_33_133
////////////////////////////////////////////////////////////////////////////////////

incidence_33_133 := SparseMatrix(Z2,0,40,[]);
full_33_133 := SparseMatrix(Z2,0,40*7,[]);

row := SparseMatrix(Z2,4,0,[]);
rf := SparseMatrix(Z2,28,0,[]);

for i in [1..5] do
	for j in [1..6] do
		for k in [1..10] do
			if k eq edges_of_ts[i][j] then
				row := HorizontalJoin(row,IdentitySparseMatrix(Z2,4));
				
				// The trivilization of a '33' intersection is the same as that appearing
				// over one of the two end points of the corresponding edge.
				// The map MPaste_Inv records the map between these trivilizations.
				// There are 4 '33' type intersections along each edge, which we combine in 'block'.
				pastes := [MPaste_Inv[[i,vertices_of_es[k][l]]] : l in [1,2]];
				block := DiagonalJoin(DiagonalJoin(DiagonalJoin(pastes[2],pastes[1]),pastes[2]),pastes[1]);
				rf := HorizontalJoin(rf,block);
			else
				row := HorizontalJoin(row,SparseMatrix(Z2,4,4,[]));
				rf := HorizontalJoin(rf,SparseMatrix(Z2,28,28,[]));
			end if;
		end for;
		incidence_33_133 := VerticalJoin(incidence_33_133,row);
		full_33_133 := VerticalJoin(full_33_133,rf);
		row := SparseMatrix(Z2,4,0,[]);
		rf := SparseMatrix(Z2,28,0,[]);
	end for;
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_13_133
////////////////////////////////////////////////////////////////////////////////////

incidence_13_133 := SparseMatrix(Z2,0,0,[]);

// In this block we record the indices of the type 3 open sets which appear in each type '33' intersection.
block := SparseMatrix(Z2,4,5,[<1,1,1>,<1,2,1>,<2,2,1>,<2,3,1>,<3,3,1>,<3,4,1>,<4,4,1>,<4,5,1>]);

for i in [1..30] do
	incidence_13_133 := DiagonalJoin(incidence_13_133,block);
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_34_134
////////////////////////////////////////////////////////////////////////////////////

incidence_34_134 := SparseMatrix(Z2,0,20,[]);
full_34_134 := SparseMatrix(Z2,0,20*7,[]);
row := SparseMatrix(Z2,2,0,[]);
rf := SparseMatrix(Z2,14,0,[]);

for i in [1..5] do
	for j in [1..6] do
		for k in [1..10] do
			if k eq edges_of_ts[i][j] then
				row := HorizontalJoin(row,IdentitySparseMatrix(Z2,2));
				
				// Open sets of type '34' have trivialization induced the
				// one fixed over the larger class '4' open set. Those of
				// type 134 are trivialized using the open set of class '1',
				// these are related by MPaste_Inv.
				pastes := [MPaste_Inv[[i,vertices_of_es[k][l]]] : l in [1,2]];
				rf := HorizontalJoin(rf,DiagonalJoin(pastes[1],pastes[2]));				
			else
				row := HorizontalJoin(row,SparseMatrix(Z2,2,2,[]));
				rf := HorizontalJoin(rf,SparseMatrix(Z2,14,14,[]));
			end if;
		end for;
		incidence_34_134 := VerticalJoin(incidence_34_134,row);
		row := SparseMatrix(Z2,2,0,[]);
		
		full_34_134 := VerticalJoin(full_34_134,rf);
		rf := SparseMatrix(Z2,14,0,[]);
	end for;
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_13_134
////////////////////////////////////////////////////////////////////////////////////

incidence_13_134 := SparseMatrix(Z2,0,0,[]);

// This block records the indicies of the open sets of 
// class '3' which appear in each of the two '34' intersections
// along each edge.
block := SparseMatrix(Z2,2,5,[<1,1,1>,<2,5,1>]);

for i in [1..30] do
	incidence_13_134 := DiagonalJoin(incidence_13_134,block);
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_14_134
////////////////////////////////////////////////////////////////////////////////////

block := SparseMatrix(Z2,12,4,[]);
incidence_14_134 := SparseMatrix(Z2,0,0,[]);

// Trivilizations induced by class '1' open set for both '14' and '134' type open sets.
// For each tetrahedron h, and each of its 6 edges i add a one to the block
// corresponding to the two vertices of the edge represented by [h,i].
// The array v_lists_2 records, for each of the 5 tetrahedra, the vertices in {1,2,3,4} of its
// six edges.
for h in [1..5] do
	for i in [1..6] do
		block[2*i-1,v_lists_2[h][i][1]] := 1;
		block[2*i,v_lists_2[h][i][2]] := 1;
	end for;
	incidence_14_134 := DiagonalJoin(incidence_14_134,block);
	block := SparseMatrix(Z2,12,4,[]);
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_22_222
////////////////////////////////////////////////////////////////////////////////////

l_block_pre := SparseMatrix(Z2,6,6,[<1,1,1>,<1,6,1>,<2,1,1>,<2,2,1>,<3,2,1>,<3,3,1>,<4,3,1>,<4,4,1>,<5,4,1>,<5,5,1>,<6,5,1>,<6,6,1>]);
l_block := SparseMatrix(Z2,0,0,[]);
for i in [1..6] do
	l_block := DiagonalJoin(l_block, l_block_pre);
end for;
l_block_full := beef_up(l_block);

r_block := SparseMatrix(Z2,36,30,[]);

for i in [1..6] do
	for j in [1..6] do
		r_block[6*(i-1)+j,bnn_matrix[i][j]] := 1;
	end for;
end for;

block := HorizontalJoin(l_block,r_block);

incidence_22_222 := SparseMatrix(Z2,0,0,[]);
full_22_222 := SparseMatrix(Z2,0,0,[]);
r_block_row := SparseMatrix(Z2,7,0,[]);
r_block_full := SparseMatrix(Z2,0,150,[]);

for i in [1..10] do
	incidence_22_222 := DiagonalJoin(incidence_22_222,block);
	
	for j in [1..36] do
		for k in [1..30] do
			if r_block[j,k] eq 1 then

				// tetra records the index of a nieghbouring tetrahedron, and
				// blob_index, the index of the region enclosed by \Delta currently.
				tetra := ts_of_fs[i][1];
				blob_index := (j-1) div 6 + 1;

				// To pass from an nn intersection to bnn (222) intersection we must pass from the 
				// trivialization over the nn intersection to the trivilization over the type 'b'
				// open set appearing in the 222 intersection. To do this we first move to a neighbouring tetrahedron
				// using sing_paste, then back to the type 'b' open set with 'blob_paste'.
				r_block_row := HorizontalJoin(r_block_row,SparseMatrix(PermutationMatrix(Z2,
					blob_paste[[tetra,i,blob_index]]))
						*sing_paste[[tetra,i,nn_direction[k],nn_parity[k]]]);
			else
				r_block_row := HorizontalJoin(r_block_row,SparseMatrix(Z2,7,5,[]));
			end if;
		end for;
			r_block_full := VerticalJoin(r_block_full, r_block_row);
			r_block_row := SparseMatrix(Z2,7,0,[]);
	end for;

	full_22_222 := DiagonalJoin(full_22_222,HorizontalJoin(l_block_full,r_block_full));
	r_block_full := SparseMatrix(Z2,0,150,[]);
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_23_233
////////////////////////////////////////////////////////////////////////////////////

incidence_23_233 := SparseMatrix(Z2,0,0,[]);
// This block records, for each '233' intersection along each edge the indices of the '23' open sets
// incident to it.
block := SparseMatrix(Z2,4,13,[<1,2,1>,<1,3,1>,<2,5,1>,<2,6,1>,<3,8,1>,<3,9,1>,<4,11,1>,<4,12,1>]);

for i in [1..30] do
	incidence_23_233 := DiagonalJoin(incidence_23_233,block);
end for;
A := incidence_23_233;

// Although there is no non-trivial transition map some of these open sets
// contain singular locus, so we cannot simply employ 'beef_up'.
beef_A := SparseMatrix(Z2,0,NumberOfColumns(A)*7 - 5*30*2,[]);
beef_row := SparseMatrix(Z2,7,0,[]);
for i in [1..NumberOfRows(A)] do
	for j in [1..NumberOfColumns(A)] do
		if A[i,j] eq 1 then
			beef_row := HorizontalJoin(beef_row,IdentitySparseMatrix(Z2,7));
		else
			if ((j-1) mod 13 +1) in singular_23 then
				beef_row := HorizontalJoin(beef_row,SparseMatrix(Z2,7,5,[]));
			else
				beef_row := HorizontalJoin(beef_row,SparseMatrix(Z2,7,7,[]));	
			end if;
		end if;
	end for;
	beef_A := VerticalJoin(beef_A,beef_row);
	beef_row := SparseMatrix(Z2,7,0,[]);
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_33_233
////////////////////////////////////////////////////////////////////////////////////

incidence_33_233 := SparseMatrix(Z2,0,40,[]);
row := SparseMatrix(Z2,4,0,[]);
for i in [1..10] do
	for j in [1..3] do
		for k in [1..10] do
			if k eq edges_of_fs[i][j] then
				row := HorizontalJoin(row,IdentitySparseMatrix(Z2,4));
			else
				row := HorizontalJoin(row,SparseMatrix(Z2,4,4,[]));
			end if;
		end for;
		incidence_33_233 := VerticalJoin(incidence_33_233,row);
		row := SparseMatrix(Z2,4,0,[]);
	end for;
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_23_234
////////////////////////////////////////////////////////////////////////////////////

incidence_23_234 := SparseMatrix(Z2,0,0,[]);
full_23_234 := SparseMatrix(Z2,0,0,[]);

// Say the 1st and 13th set in 2-3 along an edge meets a type 4 (vertex) set.
block := SparseMatrix(Z2,2,13,[<1,1,1>,<2,13,1>]);

for i in [1..30] do
	incidence_23_234 := DiagonalJoin(incidence_23_234,block);
	face_index := (i-1) div 3 +1;
	e_of_f := (i-1) mod 3 +1;
	edge_index := edges_of_fs[face_index][e_of_f];	

	tetra := ts_of_fs[face_index][1];
	block_pair := [MPaste[[tetra,vertices_of_es[edge_index][l]]]*sing_paste[[tetra,face_index,e_of_f,1]] : l in [1,2]];
	block_full := DiagonalJoin(HorizontalJoin(block_pair[1],SparseMatrix(Z2,7,8*7+15,[])),block_pair[2]);
	full_23_234 := DiagonalJoin(full_23_234,block_full);
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_24_234
////////////////////////////////////////////////////////////////////////////////////

// This is a rather arbitrary labelling of the ordering of the vertices of a triangle.
// That is of possible pairs (n1,n2) n1 ne n2.
pairs := [[1,2],[2,3],[1,3]];
block := SparseMatrix(Z2,6,3,[]);
incidence_24_234 := SparseMatrix(Z2,0,0,[]);

for i in [1..3] do
	block[2*i-1,pairs[i][1]] := 1;
	block[2*i,pairs[i][2]] := 1;
end for;

for i in [1..10] do
	incidence_24_234 := DiagonalJoin(incidence_24_234,block);
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_34_234
////////////////////////////////////////////////////////////////////////////////////

// Both 34 and 234 open sets carry trivialization inherited from the 
// larger class 4 open set. i runs over the 10 faces, j the three edges
// of this face, k runs over the 10 edges.
incidence_34_234 := SparseMatrix(Z2,0,20,[]);
row := SparseMatrix(Z2,2,0,[]);
for i in [1..10] do
	for j in [1..3] do
		for k in [1..10] do
			if k eq edges_of_fs[i][j] then
				row := HorizontalJoin(row,IdentitySparseMatrix(Z2,2));
			else
				row := HorizontalJoin(row,SparseMatrix(Z2,2,2,[]));
			end if;
		end for;
		incidence_34_234 := VerticalJoin(incidence_34_234,row);
		row := SparseMatrix(Z2,2,0,[]);
	end for;
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_22_223
////////////////////////////////////////////////////////////////////////////////////

incidence_22_223 := SparseMatrix(Z2,0,0,[]);
full_22_223 := SparseMatrix(Z2,0,0,[]);

// The incidence matrix uses only the matrix 'block'.
// i runs over the three edges of a given face.
// j runs over the 8 223 intersections along a given edge of a face.
block := SparseMatrix(Z2,24,30,[]);
for i in [1..3] do
	for j in [1..8] do
		block[8*(i-1) + j, nnp_lists[i][j]] := 1;
	end for;
end for;
block := HorizontalJoin(SparseMatrix(Z2,24,36,[]),block);

block_full := SparseMatrix(Z2,0,402,[]);
block_row := SparseMatrix(Z2,7,0,[]);
for i in [1..10] do
	incidence_22_223 := DiagonalJoin(incidence_22_223,block);

	for j in [1..24] do
		block_row := HorizontalJoin(block_row,SparseMatrix(Z2,7,36*7,[]));
		for k in [1..30] do
			if block[j,k+36] eq 1 then
				tetra := ts_of_fs[i][1];
				v_index := v_lookup[((j-1) div 8)*2+offset_22[(j-1) mod 8 + 1]+1];
				vertex := vertices_of_fs[i][v_index];
				block_row := HorizontalJoin(block_row,MPaste[[tetra,vertex]]*sing_paste[[tetra,i,nn_direction[k],nn_parity[k]]]);
			else
				block_row := HorizontalJoin(block_row,SparseMatrix(Z2,7,5,[]));
			end if;
		end for;
			block_full := VerticalJoin(block_full, block_row);
			block_row := SparseMatrix(Z2,7,0,[]);
	end for;

	full_22_223 := DiagonalJoin(full_22_223,block_full);
	block_full := SparseMatrix(Z2,0,402,[]);
end for;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_23_223
////////////////////////////////////////////////////////////////////////////////////

row_ids := [1,3,4,6,7,9,10,12];
incidence_23_223 := SparseMatrix(Z2,0,0,[]);
full_23_223 := SparseMatrix(Z2,0,0,[]);

block := SparseMatrix(Z2,8,13,[]);
for i in [1..8] do
	block[i,row_ids[i]] := 1;
	block[i,row_ids[i]+1] := 1;
end for;

block_full := SparseMatrix(Z2,0,81,[]);
block_row := SparseMatrix(Z2,7,0,[]);

for i in [1..30] do
	incidence_23_223 := DiagonalJoin(incidence_23_223,block);

	for j in [1..8] do
		for k in [1..13] do
			if block[j,k] eq 1 then
				face_index := (i-1) div 3 +1;
				e_of_f := (i-1) mod 3 + 1;
				tetra := ts_of_fs[face_index][1];
				if k in singular_23 then
					// v_index records the value in {1,2,3} of a vertex of of the face associated
					// to the class 2 open set. This vertex is chosen so that the trivialsiation
					// of \pi_\star\ZZ_2 over the face currently considered may be identified with the 
					// fixed trivialsiation at this vertex.
					v_index := v_lookup[(e_of_f-1)*2+offset_22[j]+1];
					vertex := vertices_of_fs[face_index][v_index];
					block_row := HorizontalJoin(block_row,MPaste[[tetra,vertex]]*sing_paste[[tetra,face_index,e_of_f,1]]);
				else
					block_row := HorizontalJoin(block_row,IdentitySparseMatrix(Z2,7));
				end if;
			else
				if k in singular_23 then
					block_row := HorizontalJoin(block_row,SparseMatrix(Z2,7,5,[]));
				else
					block_row := HorizontalJoin(block_row,SparseMatrix(Z2,7,7,[]));
				end if;
			end if;
		end for;
			block_full := VerticalJoin(block_full, block_row);
			block_row := SparseMatrix(Z2,7,0,[]);
	end for;

	full_23_223 := DiagonalJoin(full_23_223,block_full);
	block_full := SparseMatrix(Z2,0,81,[]);	
end for;