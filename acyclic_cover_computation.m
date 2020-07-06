////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// Here we carry out the computation of the mod 2 cohomology of the real Lagrangian.
// This uses the additional files: setup_cover.m, cech_differential_assembly_1.m,
// cech_differential_submatrices_1.m, cech_differential_0.m, ses_maps.m, and analysis.m.
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

// First we set up the arrays and matrices required to compute the Cech differential.
// In particular we record the monodromies around the singular locus.
print "We first setup the acyclic cover, this file includes various sanity checks, which should all return 'true'.\n";
load "setup_cover.m";

print "We next define and combine the Cech differential matrices out of various submatrices. These three files should produce no output.\n";

// Compute the differential relating the double and triple intersections.
load "cech_differential_submatrices_1.m";
load "cech_differential_assembly_1.m";

// Compute the differential relating the single and double intersections.
load "cech_differential_0.m";

// Compute two horizontal maps in the short exact sequence of Cech complexes.
print "We construct the horizontal maps required in the short exact sequence of Cech complexes, used to compute the tropical cycle.\n";
load "ses_maps.m"; 
print "We now process the results of our calculations.\n";

// First we compute the ranks of the matrices obtained so far.
print "Rank of first Cech differential (wrt fine cover, expect 1689):\n",Rank(full_1);
print "Rank of the second Cech differential (wrt fine cover) computing cohomology of the constant Z2 sheaf, (expect 1561):\n",Rank(incidence);
print "Rank of second Cech differential (wrt fine cover, expect 10820):\n",Rank(incidence_full);

// Next we check that the two differentials we have found form a complex.
print "Now check that the two cech differentials for the push-forward sheaf form a complex:\n", IsZero(incidence_full*full_1);

// Now we construct the connecting map beta in Cech cohomology.

// First we make an important sanity check, the image of the composition of the inclusion of C^1(R^1f_\star\ZZ_2)and
// the (second) Cech differential should have codimension 1 inside the intersection of the image of the Cech differential

smaller_space := RowSpace(Transpose(incidence_full*linear_maps_1));
larger_space := RowSpace(Transpose(incidence_full)) meet RowSpace(Transpose(linear_maps_2));
print "Dimension of the composition:\n", Rank(smaller_space);
print "Dimension of the intersection:\n", Dimension(larger_space); 
print "Is composition contained in the intersection:\n", smaller_space subset larger_space;

// Now we generate the tropical cycle via a complement of this space
tropical_cycle := Eltseq(Basis(Complement(larger_space,smaller_space))[1]);
non_zero_indices := [i : i in [1..#tropical_cycle] | tropical_cycle[i] ne 0];
support_indices := {((i-1) div 7)+1 : i in non_zero_indices};

print "The indices of the open sets contained in the support of the tropical cycle:\n",support_indices;