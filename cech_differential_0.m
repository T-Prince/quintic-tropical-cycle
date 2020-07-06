////////////////////////////////////////////////////////////////////////////////////
// Incidence_1_12
////////////////////////////////////////////////////////////////////////////////////

col := SparseMatrix(Z2,124,1,[]);
for i in [1..124] do
	col[i,1] := 1;
end for;
incidence_1_12 := SparseMatrix(Z2,0,0,[]);
for i in [1..5] do
	incidence_1_12 := DiagonalJoin(incidence_1_12,col);
end for;
incidence_1 := incidence_1_12;
incidence_1 := VerticalJoin(incidence_1,SparseMatrix(Z2,660+390+30+40+20,5,[]));

// Obtain the first full column by beefing up the old one.
full_1 := beef_up(incidence_1_12);
full_1 := VerticalJoin(full_1,SparseMatrix(Z2,4020+2430+30*7+40*7+20*7,5*7,[]));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_1_13
////////////////////////////////////////////////////////////////////////////////////

col := SparseMatrix(Z2,30,1,[]);
for i in [1..30] do
	col[i,1] := 1;
end for;
incidence_1_13 := SparseMatrix(Z2,0,0,[]);
for i in [1..5] do
	incidence_1_13 := DiagonalJoin(incidence_1_13,col);
end for;

incidence_1 := VerticalJoin(incidence_1,incidence_1_13);
full_1 := VerticalJoin(full_1,beef_up(incidence_1_13));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_1_14
////////////////////////////////////////////////////////////////////////////////////

col := SparseMatrix(Z2,4,1,[]);
for i in [1..4] do
	col[i,1] := 1;
end for;
incidence_1_14 := SparseMatrix(Z2,0,0,[]);
for i in [1..5] do
	incidence_1_14 := DiagonalJoin(incidence_1_14,col);
end for;

incidence_1 := VerticalJoin(incidence_1,incidence_1_14);
full_1 := VerticalJoin(full_1,beef_up(incidence_1_14));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_2_12
////////////////////////////////////////////////////////////////////////////////////

incidence_2_12 := SparseMatrix(Z2,0,310,[]);
full_2_12 := SparseMatrix(Z2,0,1420,[]); // 31=6+25, 6*7+25*4 = 142

row := SparseMatrix(Z2,31,0,[]);
row_full := SparseMatrix(Z2,31*7,0,[]);
for i in [1..5] do
	for j in [1..4] do
		for k in [1..10] do
			if k eq faces_of_ts[i][j] then
				row := HorizontalJoin(row,IdentitySparseMatrix(Z2,31));

				// 31=6+25 (6*7+25*4), each block is
				l_block := SparseMatrix(Z2,0,0,[]);
				for l in [1..6] do
					l_block := DiagonalJoin(l_block,SparseMatrix(PermutationMatrix(Z2,blob_paste[[i,k,l]]^-1)));
				end for;
				r_block := SparseMatrix(Z2,0,0,[]);
				for l in [1..25] do
					r_block := DiagonalJoin(r_block,sing_paste[[i,k,1,1]]*n_sing);
				end for;
				block := DiagonalJoin(l_block,r_block);
				row_full := HorizontalJoin(row_full,block);
			else
				row := HorizontalJoin(row,SparseMatrix(Z2,31,31,[]));
				row_full := HorizontalJoin(row_full,SparseMatrix(Z2,31*7,142,[]));
			end if;
		end for;
		incidence_2_12 := VerticalJoin(incidence_2_12,row);
		full_2_12 := VerticalJoin(full_2_12,row_full);
		row := SparseMatrix(Z2,31,0,[]);
		row_full := SparseMatrix(Z2,217,0,[]);
	end for;
end for;
col_2 := incidence_2_12;
fcol_2 := full_2_12;


////////////////////////////////////////////////////////////////////////////////////
// Incidence_2_22
////////////////////////////////////////////////////////////////////////////////////

incidence_2_22 := Submatrix(incidence_12_122,1,1,660,310);
full_2_22 := SparseMatrix(Z2,0,0,[]);

block := SparseMatrix(Z2,0,7,[]);
block_b_bn := SparseMatrix(Z2,0,0,[]);
block_n_bn := SparseMatrix(Z2,0,0,[]);
block_n_nn := SparseMatrix(Z2,0,0,[]);

for i in [1..6] do
	block := VerticalJoin(block, IdentitySparseMatrix(Z2,7));
end for;
for i in [1..6] do
	block_b_bn := DiagonalJoin(block_b_bn, block);
end for;

block_row := SparseMatrix(Z2,7,0,[]);
for i in [1..10] do
	block_n_bn := SparseMatrix(Z2,0,100,[]);
	block_n_nn := SparseMatrix(Z2,0,100,[]);
	
	for j in [1..36] do
		block_row := SparseMatrix(Z2,7,0,[]);
		for k in [1..25] do
			if incidence_2_22[j,k+6] eq 1 then
				tetra := ts_of_fs[i][1];
				block_row := HorizontalJoin(block_row,SparseMatrix(
					PermutationMatrix(Z2,blob_paste[[tetra,i,(j-1) mod 6 + 1]]))
						*sing_paste[[tetra,i,1,1]]*n_sing);
			else
				block_row := HorizontalJoin(block_row,SparseMatrix(Z2,7,4,[]));
			end if;
		end for;
		block_n_bn := VerticalJoin(block_n_bn, block_row);
	end for;

	for j in [1..30] do
		block_row := SparseMatrix(Z2,5,0,[]);
		for k in [1..25] do
			if incidence_2_22[j+36,k+6] eq 1 then
				block_row := HorizontalJoin(block_row,n_sing);
			else
				block_row := HorizontalJoin(block_row,SparseMatrix(Z2,5,4,[]));
			end if;
		end for;
		block_n_nn := VerticalJoin(block_n_nn, block_row);
	end for;
	block := HorizontalJoin(VerticalJoin(block_b_bn,SparseMatrix(Z2,30*5,6*7,[])),VerticalJoin(block_n_bn,block_n_nn));
	full_2_22 := DiagonalJoin(full_2_22,block);
end for;

col_2 := VerticalJoin(col_2,incidence_2_22);
fcol_2 := VerticalJoin(fcol_2,full_2_22);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_2_23
////////////////////////////////////////////////////////////////////////////////////

incidence_2_23 := Submatrix(incidence_12_123,1,1,390,310);
full_2_23 := SparseMatrix(Z2,0,0,[]);

block := SparseMatrix(Z2,0,100,[]);
block_row_1 := SparseMatrix(Z2,7,0,[]);
block_row_2 := SparseMatrix(Z2,5,0,[]);

for i in [1..10] do
	for j in [1..39] do
		is_sing := (j-1) mod 13 + 1 in singular_23;
		for k in [7..31] do
			if incidence_2_23[j,k] eq 1 then
				tetra := ts_of_fs[i][1];
				if is_sing then
					block_row_2 := HorizontalJoin(block_row_2,n_sing);
				else
					vertex := vertices_of_fs[i][ v_lookup[((j-1) div 13)*2+ offset_23[(j-1) mod 13 +1]+1]];
					block_row_1 := HorizontalJoin(block_row_1,MPaste[[tetra,vertex]]*sing_paste[[tetra,i,1,1]]*n_sing);
				end if;
			else
				block_row_1 := HorizontalJoin(block_row_1,SparseMatrix(Z2,7,4,[]));
				block_row_2 := HorizontalJoin(block_row_2,SparseMatrix(Z2,5,4,[]));
			end if;
		end for;
		if is_sing then
			block := VerticalJoin(block, block_row_2);
		else 
			block := VerticalJoin(block, block_row_1);
		end if;
		block_row_1 := SparseMatrix(Z2,7,0,[]);
		block_row_2 := SparseMatrix(Z2,5,0,[]);
	end for;
	full_2_23 := DiagonalJoin(full_2_23,HorizontalJoin(SparseMatrix(Z2,39*7-2*3*5,42,[]),block));
	block := SparseMatrix(Z2,0,100,[]);
end for;

col_2 := VerticalJoin(col_2,incidence_2_23);
fcol_2 := VerticalJoin(fcol_2,full_2_23);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_2_24
////////////////////////////////////////////////////////////////////////////////////

incidence_2_24 := Submatrix(incidence_12_124,1,1,30,310);
full_2_24 := SparseMatrix(Z2,0,0,[]);

block := SparseMatrix(Z2,0,100,[]);
block_row := SparseMatrix(Z2,7,0,[]);
for i in [1..10] do
	for j in [1..3] do
		for k in [7..31] do
			if incidence_2_24[j,k] eq 1 then
				tetra := ts_of_fs[i][1];
				block_row := HorizontalJoin(block_row,MPaste[[tetra,vertices_of_fs[i][j]]]*sing_paste[[tetra,i,1,1]]*n_sing);
			else
				block_row := HorizontalJoin(block_row,SparseMatrix(Z2,7,4,[]));
			end if;
		end for;
		block := VerticalJoin(block, block_row);
		block_row := SparseMatrix(Z2,7,0,[]);
	end for;
	full_2_24 := DiagonalJoin(full_2_24,HorizontalJoin(SparseMatrix(Z2,21,42,[]),block));
	block := SparseMatrix(Z2,0,100,[]);
end for;

col_2 := VerticalJoin(col_2,incidence_2_24);
col_2 := VerticalJoin(col_2,SparseMatrix(Z2,230,310,[]));
fcol_2 := VerticalJoin(fcol_2,full_2_24);
fcol_2 := VerticalJoin(fcol_2,SparseMatrix(Z2,230*7,1420,[]));

incidence_1 := HorizontalJoin(incidence_1,col_2);
full_1 := HorizontalJoin(full_1,fcol_2);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_3_23
////////////////////////////////////////////////////////////////////////////////////

col_3 := SparseMatrix(Z2,620+660,50,[]);
fcol_3 := SparseMatrix(Z2,620*7+4020,200,[]);

block := SparseMatrix(Z2,13,5,[<1,1,1>,<2,1,1>,<3,2,1>,<4,2,1>,<5,2,1>,<6,3,1>,<7,3,1>,<8,3,1>,<9,4,1>,<10,4,1>,<11,4,1>,<12,5,1>,<13,5,1>]);
incidence_3_23 := SparseMatrix(Z2,0,50,[]);
row := SparseMatrix(Z2,13,0,[]);
for i in [1..10] do
	for j in [1..3] do
		for k in [1..10] do
			if k eq edges_of_fs[i][j] then
				row := HorizontalJoin(row,block);
			else
				row := HorizontalJoin(row,SparseMatrix(Z2,13,5,[]));
			end if;
		end for;
		incidence_3_23 := VerticalJoin(incidence_3_23,row);
		row := SparseMatrix(Z2,13,0,[]);
	end for;
end for;

block_row_1 := SparseMatrix(Z2,5,0,[]);
block_row_2 := SparseMatrix(Z2,7,0,[]);
full_3_23 := SparseMatrix(Z2,0,200,[]);
for i in [1..390] do
	face := (i-1) div 39 + 1;
	pos := (i-1) mod 13 +1;
	for j in [1..50] do
		if incidence_3_23[i,j] eq 1 then
			edge := (j-1) div 5 + 1;
			edge_index := [p : p in [1,2,3] | edges_of_fs[face][p] eq edge][1];
			tetra := ts_of_fs[face][1];
			if pos in singular_23 then
				block_row_1 := HorizontalJoin(block_row_1,
					oneify(Transpose(ChangeRing(sing_paste[[tetra,face,edge_index,1]],Z))*
							ChangeRing(p_sing[[tetra,edge]],Z)));
			else
				vertex := vertices_of_es[edge][offset_23[pos]+1];
				block_row_2 := HorizontalJoin(block_row_2,MPaste[[tetra,vertex]]*p_sing[[tetra,edge]]);
			end if;
		else
			block_row_1 := HorizontalJoin(block_row_1,SparseMatrix(Z2,5,4,[]));
			block_row_2 := HorizontalJoin(block_row_2,SparseMatrix(Z2,7,4,[]));
		end if;
	end for;
	if pos in singular_23 then
		full_3_23 := VerticalJoin(full_3_23, block_row_1);
	else 
		full_3_23 := VerticalJoin(full_3_23, block_row_2);
	end if;
	block_row_1 := SparseMatrix(Z2,5,0,[]);
	block_row_2 := SparseMatrix(Z2,7,0,[]);
end for;

col_3 := VerticalJoin(col_3,incidence_3_23);
col_3 := VerticalJoin(col_3,SparseMatrix(Z2,30,50,[]));

fcol_3 := VerticalJoin(fcol_3,full_3_23);
fcol_3 := VerticalJoin(fcol_3,SparseMatrix(Z2,30*7,50*4,[]));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_3_33
////////////////////////////////////////////////////////////////////////////////////

incidence_3_33  := Submatrix(incidence_13_133,1,1,40,50);
full_3_33  := SparseMatrix(Z2,0,0,[]);

block := SparseMatrix(Z2,0,20,[]);
block_row := SparseMatrix(Z2,7,0,[]);
for i in [1..10] do
	for j in [1..4] do
		for k in [1..5] do
			if incidence_3_33[j,k] eq 1 then
				tetra := ts_of_es[i][1];
				block_row := HorizontalJoin(block_row,MPaste[[tetra,vertices_of_es[i][j mod 2+1]]]*p_sing[[tetra,i]]);
			else
				block_row := HorizontalJoin(block_row,SparseMatrix(Z2,7,4,[]));
			end if;
		end for;
		block := VerticalJoin(block, block_row);
		block_row := SparseMatrix(Z2,7,0,[]);
	end for;
	full_3_33 := DiagonalJoin(full_3_33,block);
	block := SparseMatrix(Z2,0,20,[]);
end for;

col_3 := VerticalJoin(col_3,incidence_3_33);
fcol_3 := VerticalJoin(fcol_3,full_3_33);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_3_34
////////////////////////////////////////////////////////////////////////////////////

incidence_3_34  := Submatrix(incidence_13_134,1,1,20,50);
full_3_34 := SparseMatrix(Z2,0,0,[]);

col_3 := VerticalJoin(col_3,incidence_3_34);

block := SparseMatrix(Z2,0,20,[]);
block_row := SparseMatrix(Z2,7,0,[]);
for i in [1..10] do
	for j in [1..2] do
		for k in [1..5] do
			if incidence_3_34[j,k] eq 1 then
				tetra := ts_of_es[i][1];
				block_row := HorizontalJoin(block_row,MPaste[[tetra,vertices_of_es[i][j]]]*p_sing[[tetra,i]]);
			else
				block_row := HorizontalJoin(block_row,SparseMatrix(Z2,7,4,[]));
			end if;
		end for;
		block := VerticalJoin(block, block_row);
		block_row := SparseMatrix(Z2,7,0,[]);
	end for;
	full_3_34 := DiagonalJoin(full_3_34,block);
	block := SparseMatrix(Z2,0,20,[]);
end for;

fcol_3 := VerticalJoin(fcol_3,full_3_34);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_3_13
////////////////////////////////////////////////////////////////////////////////////

row := SparseMatrix(Z2,5,0,[]);
incidence_3_13 := SparseMatrix(Z2,0,50,[]);
rf := SparseMatrix(Z2,5*7,0,[]);
full_3_13 := SparseMatrix(Z2,0,200,[]);

for i in [1..5] do
	for j in [1..6] do
		for k in [1..10] do
			if k eq edges_of_ts[i][j] then
				row := HorizontalJoin(row,IdentitySparseMatrix(Z2,5));
				block := SparseMatrix(Z2,0,0,[]);
				face := [f : f in faces_of_ts[i] | k in edges_of_fs[f]][1];
				for l in [1..5] do
					block := DiagonalJoin(block,p_sing[[i,k]]);
				end for;
				rf := HorizontalJoin(rf,block);
			else
				row := HorizontalJoin(row,SparseMatrix(Z2,5,5,[]));
				rf := HorizontalJoin(rf,SparseMatrix(Z2,35,20,[]));
			end if;
		end for;
		incidence_3_13 := VerticalJoin(incidence_3_13,row);
		full_3_13 := VerticalJoin(full_3_13,rf);
		row := SparseMatrix(Z2,5,0,[]);
		rf := SparseMatrix(Z2,5*7,0,[]);
	end for;
end for;

col_3 := VerticalJoin(col_3,incidence_3_13);
col_3 := VerticalJoin(col_3,SparseMatrix(Z2,20,50,[]));
incidence_1 := HorizontalJoin(incidence_1,col_3);

fcol_3 := VerticalJoin(fcol_3,full_3_13);
fcol_3 := VerticalJoin(fcol_3,SparseMatrix(Z2,20*7,200,[]));
full_1 := HorizontalJoin(full_1,fcol_3);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_4_24
////////////////////////////////////////////////////////////////////////////////////

col_4 := SparseMatrix(Z2,660+620+390,5,[]);
fcol_4 := SparseMatrix(Z2,4020+620*7+2430,5*7,[]);

incidence_4_24 := SparseMatrix(Z2,30,5,[]);

for i in [1..10] do
	for j in [1..3] do
		for k in [1..5] do
			if k eq vertices_of_fs[i][j] then
				incidence_4_24[3*(i-1)+j,k] := 1;
			end if;
		end for;
	end for;
end for;

col_4 := VerticalJoin(col_4,incidence_4_24);
col_4 := VerticalJoin(col_4,SparseMatrix(Z2,40,5,[]));

fcol_4 := VerticalJoin(fcol_4,beef_up(incidence_4_24));
fcol_4 := VerticalJoin(fcol_4,SparseMatrix(Z2,40*7,5*7,[]));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_4_34
////////////////////////////////////////////////////////////////////////////////////

incidence_4_34 := SparseMatrix(Z2,20,5,[]);

for i in [1..10] do
	for j in [1..2] do
		for k in [1..5] do
			if k eq vertices_of_es[i][j] then
				incidence_4_34[2*(i-1)+j,k] := 1;
			end if;
		end for;
	end for;
end for;

col_4 := VerticalJoin(col_4,incidence_4_34);
col_4 := VerticalJoin(col_4,SparseMatrix(Z2,150,5,[]));

fcol_4 := VerticalJoin(fcol_4,beef_up(incidence_4_34));
fcol_4 := VerticalJoin(fcol_4,SparseMatrix(Z2,150*7,5*7,[]));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_4_14
////////////////////////////////////////////////////////////////////////////////////

incidence_4_14 := SparseMatrix(Z2,20,5,[]);
full_4_14 := SparseMatrix(Z2,20*7,5*7,[]);

for i in [1..5] do
	for j in [1..4] do
		for k in [1..5] do
			if k eq vertices_of_ts[i][j] then
				incidence_4_14[4*(i-1)+j,k] := 1;
				InsertBlock(~full_4_14,MPaste_Inv[[i,k]],(4*(i-1)+j)*7-6,k*7-6);
			end if;
		end for;
	end for;
end for;

col_4 := VerticalJoin(col_4,incidence_4_14);
fcol_4 := VerticalJoin(fcol_4,full_4_14);

incidence_1 := HorizontalJoin(incidence_1,col_4);
full_1 := HorizontalJoin(full_1,fcol_4);