////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// Assemble the matrices incidence and incidence_full from the submatrices 
// constructed in cech_differential_submatrices.m. We work one row at a time
// fixed triples (abc), and joining those submatrices incidence_de_abc and
// full_de_abc such that {d,e} is a subset of {a,b,c}.
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
// Incidence_12_122
////////////////////////////////////////////////////////////////////////////////////

row_122 := incidence_12_122;
full_122 := full_12_122;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_22_122
////////////////////////////////////////////////////////////////////////////////////

row_122 := HorizontalJoin(row_122,incidence_22_122);
row_122 := HorizontalJoin(row_122,SparseMatrix(Z2,1320,390+260,[]));

full_122 := HorizontalJoin(full_122,full_22_122);
full_122 := HorizontalJoin(full_122,SparseMatrix(Z2,1320*7,2430+260*7,[]));

incidence := row_122;
incidence_full := full_122;

////////////////////////////////////////////////////////////////////////////////////
// Incidence_12_123
////////////////////////////////////////////////////////////////////////////////////

row_123 := incidence_12_123;
row_123 := HorizontalJoin(row_123,SparseMatrix(Z2,780,660,[]));
full_123 := HorizontalJoin(full_123,SparseMatrix(Z2,780*7,360*7+300*5,[]));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_23_123
////////////////////////////////////////////////////////////////////////////////////

row_123 := HorizontalJoin(row_123,incidence_23_123);
row_123 := HorizontalJoin(row_123,SparseMatrix(Z2,780,90,[]));

full_123 := HorizontalJoin(full_123,full_23_123);
full_123 := HorizontalJoin(full_123,SparseMatrix(Z2,780*7,90*7,[]));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_13_123
////////////////////////////////////////////////////////////////////////////////////

row_123 := HorizontalJoin(row_123,incidence_13_123);
row_123 := HorizontalJoin(row_123,SparseMatrix(Z2,780,20,[]));

incidence := VerticalJoin(incidence,row_123);

full_123 := HorizontalJoin(full_123,beef_up(incidence_13_123));
full_123 := HorizontalJoin(full_123,SparseMatrix(Z2,780*7,20*7,[]));
incidence_full := VerticalJoin(incidence_full,full_123);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_12_124
////////////////////////////////////////////////////////////////////////////////////

row_124 := incidence_12_124;
row_124 := HorizontalJoin(row_124,SparseMatrix(Z2,60,660+390,[]));

full_124 := HorizontalJoin(beef_up(incidence_12_124),SparseMatrix(Z2,60*7,4020+2430,[]));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_24_124
////////////////////////////////////////////////////////////////////////////////////

row_124 := HorizontalJoin(row_124,incidence_24_124);
row_124 := HorizontalJoin(row_124,SparseMatrix(Z2,60,210,[]));

full_124 := HorizontalJoin(full_124,full_24_124);
full_124 := HorizontalJoin(full_124,SparseMatrix(Z2,60*7,210*7,[]));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_14_124
////////////////////////////////////////////////////////////////////////////////////

row_124 := HorizontalJoin(row_124,incidence_14_124);
incidence := VerticalJoin(incidence,row_124);
row_133 := SparseMatrix(Z2,120,620+660+390+30,[]);

full_124 := HorizontalJoin(full_124,beef_up(incidence_14_124));
incidence_full := VerticalJoin(incidence_full,full_124);
full_133 := SparseMatrix(Z2,120*7,620*7+4020+2430+210,[]);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_33_133
////////////////////////////////////////////////////////////////////////////////////

row_133 := HorizontalJoin(row_133,incidence_33_133);
row_133 := HorizontalJoin(row_133,SparseMatrix(Z2,120,20));

full_133 := HorizontalJoin(full_133,full_33_133);
full_133 := HorizontalJoin(full_133,SparseMatrix(Z2,120*7,20*7));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_13_133
////////////////////////////////////////////////////////////////////////////////////

row_133 := HorizontalJoin(row_133,incidence_13_133);
row_133 := HorizontalJoin(row_133,SparseMatrix(Z2,120,20,[]));
incidence := VerticalJoin(incidence, row_133);
row_134 := SparseMatrix(Z2,60,620+660+390+70,[]);

full_133 := HorizontalJoin(full_133,beef_up(incidence_13_133));
full_133 := HorizontalJoin(full_133,SparseMatrix(Z2,120*7,20*7,[]));
incidence_full := VerticalJoin(incidence_full, full_133);
full_134 := SparseMatrix(Z2,60*7,620*7+4020+2430+70*7,[]);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_34_134
////////////////////////////////////////////////////////////////////////////////////

row_134 := HorizontalJoin(row_134,incidence_34_134);
full_134 := HorizontalJoin(full_134,full_34_134);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_13_134
////////////////////////////////////////////////////////////////////////////////////

row_134 := HorizontalJoin(row_134,incidence_13_134);
full_134 := HorizontalJoin(full_134,beef_up(incidence_13_134));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_14_134
////////////////////////////////////////////////////////////////////////////////////

row_134 := HorizontalJoin(row_134,incidence_14_134);
incidence := VerticalJoin(incidence,row_134);
row_222 := SparseMatrix(Z2,360,620);

full_134 := HorizontalJoin(full_134, beef_up(incidence_14_134));
incidence_full := VerticalJoin(incidence_full,full_134);
full_222 := SparseMatrix(Z2,360*7,620*7);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_22_222
////////////////////////////////////////////////////////////////////////////////////

row_222 := HorizontalJoin(row_222,incidence_22_222);
row_222 := HorizontalJoin(row_222,SparseMatrix(Z2,360,390+260,[]));
incidence := VerticalJoin(incidence,row_222);

full_222 := HorizontalJoin(full_222,full_22_222);
full_222 := HorizontalJoin(full_222,SparseMatrix(Z2,360*7,2430+260*7,[]));
incidence_full := VerticalJoin(incidence_full,full_222);

row_233 := SparseMatrix(Z2,120,1280,[]);
full_233 := SparseMatrix(Z2,120*7,620*7+4020,[]);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_23_233
////////////////////////////////////////////////////////////////////////////////////

row_233 := HorizontalJoin(row_233,incidence_23_233);
row_233 := HorizontalJoin(row_233,SparseMatrix(Z2,120,30,[]));

full_233 := HorizontalJoin(full_233,beef_A);
full_233 := HorizontalJoin(full_233,SparseMatrix(Z2,120*7,30*7,[]));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_33_233
////////////////////////////////////////////////////////////////////////////////////

row_233 := HorizontalJoin(row_233,incidence_33_233);
row_233 := HorizontalJoin(row_233,SparseMatrix(Z2,120,190,[]));
incidence := VerticalJoin(incidence,row_233);

full_233 := HorizontalJoin(full_233,beef_up(incidence_33_233));
full_233 := HorizontalJoin(full_233,SparseMatrix(Z2,120*7,190*7,[]));
incidence_full := VerticalJoin(incidence_full,full_233);

row_234 := SparseMatrix(Z2,60,1280,[]);
full_234 := SparseMatrix(Z2,60*7,620*7+4020,[]);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_23_234
////////////////////////////////////////////////////////////////////////////////////

row_234 := HorizontalJoin(row_234,incidence_23_234);
full_234 := HorizontalJoin(full_234,full_23_234);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_24_234
////////////////////////////////////////////////////////////////////////////////////

row_234 := HorizontalJoin(row_234,incidence_24_234);
row_234 := HorizontalJoin(row_234,SparseMatrix(Z2,60,40,[]));

full_234 := HorizontalJoin(full_234,beef_up(incidence_24_234));
full_234 := HorizontalJoin(full_234,SparseMatrix(Z2,60*7,40*7,[]));

////////////////////////////////////////////////////////////////////////////////////
// Incidence_34_234
////////////////////////////////////////////////////////////////////////////////////

row_234 := HorizontalJoin(row_234,incidence_34_234);
row_234 := HorizontalJoin(row_234,SparseMatrix(Z2,60,170,[]));
incidence := VerticalJoin(incidence,row_234);
row_223 := SparseMatrix(Z2,240,620,[]);

full_234 := HorizontalJoin(full_234,beef_up(incidence_34_234));
full_234 := HorizontalJoin(full_234,SparseMatrix(Z2,60*7,170*7,[]));
incidence_full := VerticalJoin(incidence_full,full_234);
full_223 := SparseMatrix(Z2,240*7,620*7,[]);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_22_223
////////////////////////////////////////////////////////////////////////////////////

row_223 := HorizontalJoin(row_223,incidence_22_223);
full_223 := HorizontalJoin(full_223,full_22_223);

////////////////////////////////////////////////////////////////////////////////////
// Incidence_23_223
////////////////////////////////////////////////////////////////////////////////////

row_223 := HorizontalJoin(row_223,incidence_23_223);
row_223 := HorizontalJoin(row_223,SparseMatrix(Z2,240,260,[]));
incidence := VerticalJoin(incidence,row_223);

full_223 := HorizontalJoin(full_223,full_23_223);
full_223 := HorizontalJoin(full_223,SparseMatrix(Z2,240*7,260*7,[]));
incidence_full := VerticalJoin(incidence_full,full_223);