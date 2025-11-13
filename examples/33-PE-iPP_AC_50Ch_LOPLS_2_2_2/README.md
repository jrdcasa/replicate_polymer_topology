# Typing molecular systems with different forcefields (OPLS and L-OPLS)

The current system is composed by polyethylene (PE) and isotactic polypropylene (iPP) chains. In this example, we want use the [LOPLS force field](https://pubs.acs.org/doi/10.1021/ct200908r) for PE and the [OPLS force field](https://pubs.acs.org/doi/10.1021/ja9621760) for iPP. Unfortunately, it is very difficult to generate different semantics based on SMART to type these systems using these force fields within the XML file used by [replicate_polymer](https://github.com/jrdcasa/replicate_polymer_topology). Essentially, these is because the LOPLS FF is a reparametrization for linear hydrocarbons, which do not define new types, but change the non-bonded parameters of some of them.

As the [replicate_polymer](https://github.com/jrdcasa/replicate_polymer_topology) software use the SMART semantic defined in [Foyer](https://github.com/mosdef-hub/foyer/issues/63), the direct use of the program does not assign correctly the parameters to iPP.

A shortcut to type the atoms of different polymer chains with distinct but related force fields, as L-OPLS and OPLS, is presentent in this tutorial.

## **1.** Initial structure

----

The initial structure is prepared using the amorphous builder module in Materials Studio. The system contains 25 PE chains with 40 monomers each (C $_{80}$ H $_{162}$) and 25 iPP chains with 40 monomers (C $_{120}$ H $_{242}$, 80 backbone carbons). The XSD structure can be found [here](./01-PREPARE/iPP-PE_40mon_25-25Ch.xsd).

The script [01-topology_script.py](./01-PREPARE/01-topology_script.py) along with the information of the head and tail atoms ([iPP-PE_40mon_25-25Ch_headtail.dat](./01-PREPARE/iPP-PE_40mon_25-25Ch_headtail.dat)) and residues ([iPP-PE_40mon_25-25Ch_residues.dat](./01-PREPARE/iPP-PE_40mon_25-25Ch_residues.dat)) can be used with the [topology library](https://github.com/jrdcasa/topology.git) to put order and assign residues to the whole system. Compare the [original numbering](./01-PREPARE/iPP-PE_40mon_25-25Ch_order.pdb) with [the ordered one](./01-PREPARE/iPP-PE_40mon_25-25Ch_residues.pdb). The commands to produce the files from the xsd structure are:

```bash
source ~/Programacion/sandbox/sandbox_common/bin/activate
python 01-topology_script.py
```

The first 25 chains in the pdb file ([iPP-PE_40mon_25-25Ch_residues.pdb](./01-PREPARE/iPP-PE_40mon_25-25Ch_residues.pdb)) corresponds to PE polymers whereas the 25 last chains are the isotactic PP polymers.

## **2.** Typing the molecule

----

The [file](./01-PREPARE/iPP-PE_40mon_25-25Ch_residues.pdb) is suitable to be typed using [replicate_polymer](https://github.com/jrdcasa/replicate_polymer_topology) software. If the program is applied directly, both PE and iPP topology files will be typed with L-OPLS type atoms, which is not valid for iPP polymer (see directory [99-DIRECT_NOVALID](./99-DIRECT_NOVALID/)). However, the [coordinate file](./99-DIRECT_NOVALID/iPP-PE_40mon_25-25Ch_residues_replicate.gro) is valid and it will be used in the next steps. In this case, we produce a 2x2x2 replica of the system:

```bash
replicate_polymer -p iPP-PE_40mon_25-25Ch_residues.pdb -f ../../../forcefields/oplsaa.xml --images 2 2 2 -e lammps
```

The trick to apply is as follows:

### **2.a.** Separate one polymer for each chain

----

Using VMD load the [01-PREPARE/iPP-PE_40mon_25-25Ch_residues.pdb](./01-PREPARE/iPP-PE_40mon_25-25Ch_residues.pdb). Then, save a PE polymer chain (i.e. selecting the fragment 0) and one iPP polymer chain (i.e. selecting the fragment 25). Move or copy the generate files to the folders [PE_SC.pdb](./02-PE_SC/PE_SC.pdb) and [iPP_SC.pdb](./03-iPP_SC/iPP_SC.pdb). These single chains will be used to create the topology files for each polymer chain separately by using [replicate_polymer](https://github.com/jrdcasa/replicate_polymer_topology).

```bash
cd 02-PE_SC
replicate_polymer -p PE_SC.pdb -f ../../../forcefields/loplsaa.xml --images 1 1 1
cd 03-iPP_SC
replicate_polymer -p iPP_SC.pdb -f ../../../forcefields/oplsaa.xml --images 1 1 1
```

### **2.b.** Combining files

----

Then, the topology files [PE_SC_con_replicate.top](./02-PE_SC/PE_SC_con_replicate.top) and [iPP_SC_con_replicate.top](./03-iPP_SC/iPP_SC_con_replicate.top) will be combined using a text editor.

Keep in mind the following things:

1. [atomtypes] section must only appear once in the top file.

2. In both topoogy files the [moleculetype] is defined as **system1**. Correct [moleculetype] sections:

    ```test
    [ moleculetype ]
    ; Name            nrexcl
    system1          3
    (...)
    [ moleculetype ]
    ; Name            nrexcl
    system2          3
    ```

3. Change the [molecules] section in this case:

    ```test
    [ molecules ]
    system1 25   # PE
    system2 25   # iPP
    system1 25   # PE
    system2 25   # iPP
    system1 25   # PE
    system2 25   # iPP
    system1 25   # PE
    system2 25   # iPP
    system1 25   # PE
    system2 25   # iPP
    system1 25   # PE
    system2 25   # iPP
    system1 25   # PE
    system2 25   # iPP
    system1 25   # PE
    system2 25   # iPP
    ```

The new topology file is called [iPP-PE_40mon_25-25Ch_residues_replicate.top](./iPP-PE_40mon_25-25Ch_residues_replicate.top) Finally, copy the gro file from [99-DIRECT_NOVALID](./99-DIRECT_NOVALID/).

Check the vality of the files running a GROMACS simulation [minim_template.mdp](../00-TEMPLATES/minim_template.mdp)

### **2.c.** Creating the lammps files

----

Last step is to create the LAMMPS files using the utility ``create_lammps_from_gromacs``

```bash
create_lammps_from_gromacs -g iPP-PE_40mon_25-25Ch_residues_replicate.gro -t iPP-PE_40mon_25-25Ch_residues_replicate.top -f ../../forcefields/oplsaa.xml
```
