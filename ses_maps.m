// This function is used to detemine the map G^\vee to G' stalk by stalk.
// It takes in a 7x5 block representing a matrix stored in sing_paste.
mono_index := function(A)
	pos := Matrix(Z2,[[0,1,1,1,0,0,1],[0,0,0,1,1,1,1],[1,1,0,0,0,1,1]]);
	
	indices := [i : i in [1..7] | A[i,4] ne 0];
	if #indices ne 2 then
		print "Error: number of identified sheets different from 2", Matrix(A);
	end if;
		big_index := [k : k in [1..7] | ColumnSubmatrix(pos,k,1) eq  ColumnSubmatrix(pos,indices[1],1)+ColumnSubmatrix(pos,indices[2],1)][1];
		index := [l : l in [1,2,3] | A[big_index,l] ne 0][1];
	return index;
end function;

// Compute C^1(\cU,G^\vee).
linear_maps_1 := SparseMatrix(Z2,0,0,[]);
block := SparseMatrix(Z2,7,3,[<1,1,1>,<2,1,1>,<6,1,1>,<7,1,1>,<2,2,1>,<3,2,1>,<4,2,1>,<7,2,1>,<4,3,1>,<5,3,1>,<6,3,1>,<7,3,1>]);

// Note that this assumes the columns with a pair of entries are always the final two.
block_sings := [SparseMatrix(Z2,5,2,[<2,1,1>,<3,1,1>,<4,1,1>,<4,2,1>,<5,2,1>]),SparseMatrix(Z2,5,2,[<1,1,1>,<3,1,1>,<4,1,1>,<4,2,1>,<5,2,1>]),SparseMatrix(Z2,5,2,[<1,1,1>,<2,1,1>,<4,1,1>,<4,2,1>,<5,2,1>])];

for i in [1..NumberOfColumns(incidence_12_122)] do
	linear_maps_1 := DiagonalJoin(linear_maps_1,block);
end for;

// Next determine whether we are in a nn or bn case etc.

linear_maps_1_22 := SparseMatrix(Z2,0,0,[]);
block_bn := SparseMatrix(Z2,0,0,[]);
block_nn := SparseMatrix(Z2,0,0,[]);
depth := 0;
height := 0;
for k in [1..10] do
	for l in [1..36] do
		block_bn := DiagonalJoin(block_bn, block);
	end for;
	linear_maps_1_22 := DiagonalJoin(linear_maps_1_22,block_bn);	
	for l in [1..30] do
		height := NumberOfRows(linear_maps_1_22)+5*(l-1)+1;
		depth := Min([d : d in [1..NumberOfRows(full_22_122)] | Submatrix(full_22_122,d,height,1,5) ne 0]);
		block_nn := DiagonalJoin(block_nn, block_sings[mono_index(Submatrix(full_22_122,depth,height,7,5))]);	
	end for;

	linear_maps_1_22 := DiagonalJoin(linear_maps_1_22,block_nn);
	block_bn := SparseMatrix(Z2,0,0,[]);
	block_nn := SparseMatrix(Z2,0,0,[]);
end for;
linear_maps_1 := DiagonalJoin(linear_maps_1,linear_maps_1_22);

// Now for the 23 case
linear_maps_1_23 := SparseMatrix(Z2,0,0,[]);
block_23 := SparseMatrix(Z2,0,0,[]);
local_height := 1;

for k in [1..10] do
	local_height := 1;
	for l in [1..39] do
		pos := (l-1) mod 13 + 1;

		// Detect whether the current np type face intersects the singular locus.
		if pos in singular_23 then

			// Store the current position in the matrix full_23_123
			height := NumberOfRows(linear_maps_1_23)+local_height;
			depth := Min([d : d in [1..NumberOfRows(full_23_123)] | Submatrix(full_23_123,d,height,1,5) ne 0]);
			
			// Add a 7x3 block corresponding to the given segment of discriminant locus.
			block_23 := DiagonalJoin(block_23,block_sings[mono_index(Submatrix(full_23_123,depth,height,7,5))]);
			local_height := local_height + 5;
		else
			block_23 := DiagonalJoin(block_23,block);
			local_height := local_height + 7;
		end if;
	end for;
	linear_maps_1_23 := DiagonalJoin(linear_maps_1_23,block_23);
	block_23 := SparseMatrix(Z2,0,0,[]);
end for;

linear_maps_1 := DiagonalJoin(linear_maps_1,linear_maps_1_23);

for i in [1..(NumberOfColumns(incidence) - NumberOfColumns(incidence_12_122) - NumberOfColumns(incidence_22_122) - NumberOfColumns(incidence_23_123))] do
	linear_maps_1 := DiagonalJoin(linear_maps_1,block);
end for;


// Now we compute C^2(\cU,G^\vee): this is straightforward.

linear_maps_2 := SparseMatrix(Z2,0,0,[]);
for i in [1..NumberOfRows(incidence)] do
	linear_maps_2 := DiagonalJoin(linear_maps_2,block);
end for;