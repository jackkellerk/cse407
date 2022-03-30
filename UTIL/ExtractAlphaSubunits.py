# TODO: KCN pdb files do not seem to have alpha subunits?
# TODO: When extracting subunits, are we extracting SEQRES or ATOM lines?
# I.e. SEQRES   1 ??A??   53  ASP ASN PHE ASN GLN GLN LYS LYS LYS PHE GLY GLY GLN    
# Or   ATOM     1  ??N??   GLN **A**1486      12.961  -6.234   7.725  1.00  0.00    ??N??  
# What am i looking for?

# Fifth column that starts with ATOMS (at least for SCN) (cna be like B12345 would be same as B)