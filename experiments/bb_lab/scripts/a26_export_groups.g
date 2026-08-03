# A26: export multiplication tables for the mitten-code groups
# (arXiv:2607.28795 Table VII/XIII), in GAP's Elements(G) ordering — the
# paper's element-index convention. Regenerates instances/mitten_groups/*.txt.
#
# Run:  gap -q --nointeract a26_export_groups.g
#
# Output format per file:
#   line 1: n;  line 2: StructureDescription;  line 3: element orders;
#   line 4: inverse table (0-based);  lines 5..n+4: multiplication table rows
#   (0-based; row i, col j = index of e_i * e_j).
# NOTE: GAP wraps long lines with backslash continuations; the Python loader
# strips them ("\\\n" -> "").

ExportGroup := function(G, fname)
    local els, n, i;
    els := Elements(G);
    n := Size(G);
    PrintTo(fname, n, "\n");
    AppendTo(fname, StructureDescription(G), "\n");
    AppendTo(fname, JoinStringsWithSeparator(List([1..n], i -> String(Order(els[i]))), " "), "\n");
    AppendTo(fname, JoinStringsWithSeparator(List([1..n], i -> String(Position(els, els[i]^-1) - 1)), " "), "\n");
    for i in [1..n] do
        AppendTo(fname, JoinStringsWithSeparator(List([1..n], j -> String(Position(els, els[i]*els[j]) - 1)), " "), "\n");
    od;
    Print("wrote ", fname, "\n");
end;

for id in [ [30,1], [40,5], [60,11], [100,9], [108,9], [126,1], [156,13], [195,1] ] do
    ExportGroup(SmallGroup(id[1], id[2]),
                Concatenation("group_", String(id[1]), "_", String(id[2]), ".txt"));
od;

# The shipped processor_codes [[300,60,14]] lives in DirectProduct(C10,S3)
# Elements order (== row-major Kronecker order), not SmallGroup(60,11) order.
ExportGroup(DirectProduct(CyclicGroup(10), SymmetricGroup(3)), "group_dp_c10_s3.txt");

QUIT;
