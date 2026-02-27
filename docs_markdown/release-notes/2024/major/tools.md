# Improvements to GROMACS tools

## Improved Einstein viscosity calculation in gmx energy

Viscosity calculation using the Einstein formula is convenient as this
does not require extremely frequent pressure tensor data. However, the
implementation of the calculation was inconvenient for long simulations
and could take hours to complete. Improved stepping through the data
reduces the computational time to minutes and provides much clearer
output.

## XVG output from `gmx rdf` now uses 6 decimal places

The output from `gmx rdf` now uses more decimal places in order to avoid
rounding issues. These issues led to perceived erroneous shifts in the
results.

`4647`

## Handle CYX-CYX disulfide bonds in `gmx pdb2gmx`

Naming CYS residues as CYX shows that they should form a disulfide bond.
`gmx pdb2gmx` will now correctly interpret them as disulfide bond
forming residues.

`4929`
