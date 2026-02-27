# Answers to frequently asked questions (FAQs)

## Questions regarding GROMACS installation

1.  Do I need to compile all utilities with MPI?

    With one rarely-used exception (`pme_error `), only
    `mdrun ` is able to use the `MPI <mpi-support>`
    parallelism. So you only need to use the `-DGMX_MPI=on` flag when
    `configuring <configure-cmake>` for a build intended to run the main
    simulation engine `mdrun `. Generally that is desirable
    when running on a multi-node cluster, and necessary when using
    multi-simulation algorithms. Usually also installing a build of
    GROMACS configured without MPI is
    convenient for users.

2.  Should my version be compiled using double precision?

    In general, GROMACS only needs to be build
    in its default mixed-precision mode. For more details, see the
    discussion in Chapter 2 of the [reference manual](). Sometimes,
    usage may also depend on your target system, and should be decided
    upon according to the `individual instructions `.

## Questions concerning system preparation and preprocessing

1.  Where can I find a solvent `coordinate file `
    for use with `solvate `?

    Suitable equilibrated boxes of solvent
    `structure files ` can be found in the
    `$GMXDIR/share/gromacs/top` directory. That location will be
    searched by default by `solvate `, for example by using
    `-cs spc216.gro` as an argument. Other solvent boxes can be prepared
    by the user as described on the manual page for
    `solvate ` and elsewhere. Note that suitable topology
    files will be needed for the solvent boxes to be useful in
    `grompp `. These are available for some force fields,
    and may be found in the respective subfolder of
    `$GMXDIR/share/gromacs/top`.

2.  How to prevent `solvate ` from placing waters in
    undesired places?

    Water placement is generally well behaved when solvating proteins,
    but can be difficult when setting up membrane or micelle
    simulations. In those cases, waters may be placed in between the
    alkyl chains of the lipids, leading to problems later
    `during the simulation <blowing-up>`. You can either remove those
    waters by hand (and do the accounting for molecule types in the
    `topology ` file), or set up a local copy of the `vdwradii.dat`
    file from the `$GMXLIB` directory, specific for your project and
    located in your working directory. In it, you can increase the vdW
    radius of the atoms, to suppress such interstitial insertions.
    Recommended e.g. at a common
    [tutorial](http://www.mdtutorials.com/gmx/lysozyme/03_solvate.html)
    is the use of 0.375 instead of 0.15.

<!-- -->

1.  How do I provide multiple definitions of bonds / dihedrals in a
    topology?

    You can add additional bonded terms beyond those that are normally
    defined for a residue (e.g. when defining a special ligand) by
    including additional copies of the respective lines under the
    `[ bonds ]`, `[ pairs ]`, `[ angles ]` and `[ dihedrals ]` sections
    in the `[ moleculetype ]` section for your molecule, found either in
    the `itp` file or the `topology ` file. This will **add** those
    extra terms to the potential energy evaluation, but **will not**
    remove the previous ones. So be careful with duplicate entries. Also
    keep in mind that this **does not** apply to duplicated entries for
    `[ bondtypes ]`, `[ angletypes ]`, or `[ dihedraltypes ]`, in
    force-field definition files, where duplicates overwrite the
    previous values.

2.  Do I really need a `gro` file?

    The `gro` file is used in GROMACS as a
    unified `structure file ` format that can be
    read by all utilities. The large majority of
    GROMACS routines can also use other file
    types such as `pdb`, with the limitations that no velocities are
    available in `this case `. If you need a
    text-based format with more digits of precision, the `g96` format is
    suitable and supported.

3.  Do I always need to run `pdb2gmx ` when I already
    produced an `itp` file elsewhere?

    You don't need to prepare additional files if you already have all
    `itp` and `top` files prepared through other tools.

    Examples for those can be found in the
    `System Preparation section of this user guide <../user-guide/system-preparation>`.

4.  How can I build in missing atoms?

    GROMACS has no support for building
    coordinates of missing non-hydrogen atoms. If your system is missing
    some part, you will have to add the missing pieces using external
    programs to avoid the `missing atom ` error. This
    can be done using programs such as
    [Chimera](https://www.cgl.ucsf.edu/chimera/) in combination with
    [Modeller](https://salilab.org/modeller/), [Swiss PDB
    Viewer](https://spdbv.unil.ch/),
    [Maestro](https://www.schrodinger.com/maestro). **Do not run** a
    simulation that had missing atoms unless you know exactly why it
    will be stable.

5.  Why is the total charge of my system not an integer like it should
    be?

    In `floating point ` math, real numbers can not
    be displayed to arbitrary precision (for more on this, see e.g.
    [Wikipedia](https://en.wikipedia.org/wiki/Floating-point_arithmetic)).
    This means that very small differences to the final integer value
    will persist, and GROMACS will not lie to
    you and round those values up or down. If your charge differs from
    the integer value by a larger amount, e.g. at least 0.01, this
    usually means that something went wrong during your system
    preparation

## Questions regarding simulation methodology

1.  Should I couple a handful of ions to their own temperature-coupling
    bath?

    **No**. You need to consider the minimal size of your temperature
    coupling groups, as explained in `gmx-thermostats` and more
    specifically in `gmx-thermostats-dont`, as well as the
    implementation of your chosen thermostat as described in the
    [reference manual]().

2.  Why do my grompp restarts always start from time zero?

    You can choose different values for `tinit` and `init-step`.

    > **TODO:** Add "Continuing simulations" content (label: gmx-cont-simulation)
    and link.

    e.g. `` :ref:`Continuing simulations `. ``

3.  Why can't I do conjugate gradient minimization with constraints?

    Minimization with the conjugate gradient scheme can not be performed
    with constraints as described in the [reference manual](), and some
    additional information on
    [Wikipedia](https://en.wikipedia.org/wiki/Conjugate_gradient_method).

4.  How do I hold atoms in place in my energy minimization or
    simulation?

    Groups may be frozen in place using `freeze groups` (see the
    [reference manual]()). It is more common to use a set of position
    restraints, to place penalties on movement of the atoms. Files that
    control this kind of behaviour can be created using
    `genrestr `.

5.  How do I extend a completed a simulation to longer times?

    Please see the section on `managing long simulations`. You can
    either prepare a new `mdp` file, or extend the simulation time in
    the original `tpr` file using `convert-tpr `.

    1.  How do I complete a crashed simulation?

    Need gmx-cont-crash doc target.

    ``` none
    This can be easily achieved using the checkpoint reading
    :ref:`available ` in |Gromacs| versions newer than 4.
    ```

    1.  How can I do a simulation at constant pH?

    Need gmx-howto-cph doc target.

    > ``` none
    > This is a rather large topic, and you should at least read the short
    > :ref:`Constant pH How-To ` and all of the literature
    > included there to get an overview over the topic.
    > ```

6.  How should I compute a single-point energy?

    This is best achieved with the `-rerun` option to
    `mdrun `. See the `single-point energy` section.

## Parameterization and Force Fields

1.  I want to simulate a molecule (protein, DNA, etc.) which complexes
    with various transition metal ions, iron-sulfur clusters, or other
    exotic species. Parameters for these exotic species aren't available
    in force field X. What should I do?

    First, you should consider how well `MD ` will actually
    describe your system (e.g. see some of the [recent
    literature](https://dx.doi.org/10.1021%2Facs.chemrev.6b00440)). Many
    species are infeasible to model without either atomic
    polarizability, or QM treatments. Then you need to prepare your own
    set of parameters and add a new residue to your
    `force field ` of choice. Then you will have to
    validate that your system behaves in a physical way, before
    continuing your simulation studies. You could also try to build a
    more simplified model that does not rely on the complicated
    additions, as long as it still represents the correct *real* object
    in the laboratory.

2.  Should I take parameters from one force field and apply them inside
    another that is missing them?

    **NO**. Molecules parametrized for a given
    `force field ` will not behave in a physical manner
    when interacting with other molecules that have been parametrized
    according to different standards. If your required molecule is not
    included in the force field you need to use, you will have to
    parametrize it yourself according to the methodology of this force
    field.

## Analysis and Visualization

1.  How do I visualize a trajectory?

gmx-howto-visualize doc target:

``` none
Use one of the number of different programs that can visualize
coordinate :ref:`files and trajectories `.
```


1.  Why am I seeing bonds being created when I watch the trajectory?

    Most visualization softwares determine the bond status of atoms
    depending on a set of predefined distances. So the bonding pattern
    created by them might not be the one defined in your
    `topology ` file. What matters is the information encoded in
    there. If the software has read a `tpr <tpr>` file, then the
    information is in reliable agreement with the topology you supplied
    to `grompp `.

2.  When visualizing a trajectory from a simulation using PBC, why are
    there holes or my peptide leaving the simulation box?

    Those holes and molecules moving around are just a result of
    molecules ranging over the
    `box boundaries and wrapping around `, and are not a reason
    for concern. You can fix the visualization using
    `trjconv ` to prepare the structure for analysis.

3.  Why is my total simulation time not an integer like it should be?

    As the simulation time is calculated using
    `floating point arithmetic `, rounding errors
    can occur but are not of concern.
