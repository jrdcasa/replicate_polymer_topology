# Typing sytems with different forcefields.

The current system is composed by polyethylene (PE) and isotactic polypropylene (iPP) chains. In this example, we want use the [LOPLS force field](https://pubs.acs.org/doi/10.1021/ct200908r) for PE and the [OPLS force field](https://pubs.acs.org/doi/10.1021/ja9621760) for iPP. Unfortunately, it is very difficult to generate different semantics based on SMART to type these systems using these force fields within the XML file used by [replicate_polymer](https://github.com/jrdcasa/replicate_polymer_topology). Essensially, these is because the LOPLS FF is a reparametrization for linear hydrocarbons, which do not define new types.
The [replicate_polymer](https://github.com/jrdcasa/replicate_polymer_topology) software use the SMART semantic defined in [Foyer](https://github.com/mosdef-hub/foyer/issues/63). The direct use of the program does not assign correctly the parameters to iPP.

A shortcut to type the atoms of different polymer chains with distinct but related force fields, as L-OPLS and OPLS, is presentent in this tutorial.

## 1. Initial structure

----

The initial structure is prepared using the amorphous builder module in Materials Studio. The system contains 25 PE chains with 40 monomers each (C $_{80}$ H $_{162}$) and 25 iPP chains with 40 monomers (C $_{120}$ H $_{242}$, 80 backbone carbons). The XSD structure can be found [here](./01-PREPARE/iPP-PE_40mon_25-25Ch.xsd).

The script [01-topology_script.py](./01-PREPARE/01-topology_script.py) along with the information of the head and tail atoms ([iPP-PE_40mon_25-25Ch_headtail.dat](./01-PREPARE/iPP-PE_40mon_25-25Ch_headtail.dat)) and residues ([iPP-PE_40mon_25-25Ch_residues.dat](./01-PREPARE/iPP-PE_40mon_25-25Ch_residues.dat)) can be used with the topology library to put order and assign residues to the whole system. Compare the [original numbering](./01-PREPARE/iPP-PE_40mon_25-25Ch_order.pdb) with [the ordered one](./01-PREPARE/iPP-PE_40mon_25-25Ch_residues.pdb). The commands to produce the files from the xsd structure are:

```bash
source ~/Programacion/sandbox/sandbox_common/bin/activate
python 01-topology_script.py
```

The first 25 chains in the pdb file ([iPP-PE_40mon_25-25Ch_residues.pdb](./01-PREPARE/iPP-PE_40mon_25-25Ch_residues.pdb)) corresponds to PE polymers whereas the 25 last chains are the isotactic PP polymers.

## 2. Typing the molecule

----

The [file](./01-PREPARE/iPP-PE_40mon_25-25Ch_residues.pdb) is suitable to be typed using [replicate_polymer](https://github.com/jrdcasa/replicate_polymer_topology) software.