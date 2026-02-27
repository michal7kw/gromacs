# Dihedral principal component analysis

`gmx angle `, `gmx covar `,
`gmx anaeig `  
Principal component analysis can be performed in dihedral
space `172 <refMu2005a>` using GROMACS. You
start by defining the dihedral angles of interest in an index file,
either using `gmx mk_angndx ` or otherwise. Then you use
the `gmx angle ` program with the `-or` flag to produce a new
`trr` file containing the cosine and sine of each dihedral angle in two
coordinates, respectively. That is, in the `trr` file you will have a
series of numbers corresponding to: cos($\phi_1$), sin($\phi_1$),
cos($\phi_2$), sin($\phi_2$), ..., cos($\phi_n$), sin($\phi_n$),
and the array is padded with zeros, if necessary. Then you can use this
`trr` file as input for the `gmx covar ` program and perform
principal component analysis as usual. For this to work you will need to
generate a reference file (`tpr`, `gro`, `pdb` etc.) containing the same
number of “atoms” as the new `trr` file, that is for $n$ dihedrals you
need 2$n$/3 atoms (rounded up if not an integer number). You should
use the `-nofit` option for `gmx covar ` since the
coordinates in the dummy reference file do not correspond in any way to
the information in the `trr` file. Analysis of the results is done using
`gmx anaeig `.
