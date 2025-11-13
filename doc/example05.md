# Example 5: Replicate a trimer P4HB from a single chain (OPLS)

**The files to run this example are available in the [examples](../examples/05-Trimer_P4HB_TrappeUA/) folder inside the distribution**

```bash
cp ../04-Trimer_P4HB_OPLSAA/TrimerP4HB_residues.pdb .

replicate_polymer -p TrimerP4HB_residues.pdb -f ../../forcefields/trappe-ua.xml --images 1 1 1  -e lammps --boxlength 2.3 2.3 2.3 --noh --verbose
```

In this case, a lot of errors is obtained due to there is not parameters for some dihedral types.

```text
		 Checking angles    ... (verbose: True) (25-07-2022 18:09:53).
			 0 angles of 22
			 4 angles of 22
			 8 angles of 22
			 12 angles of 22
			 16 angles of 22
			 20 angles of 22
		 End checking angles ... (verbose: True) (25-07-2022 18:09:53).
		 Checking dihedrals    ... (verbose: True) (25-07-2022 18:09:53).
			 0 dihedral of 21
		ERROR: Missing dihedral with ids (2, 3, 4, 5) and types ['CH2_OH', 'CH2_sp3', 'CH2_sp3', 'COO'].
		ERROR: Missing dihedral with ids (3, 4, 5, 6) and types ['CH2_sp3', 'CH2_sp3', 'COO', 'OC'].
			 4 dihedral of 21
		ERROR: Missing dihedral with ids (3, 4, 5, 7) and types ['CH2_sp3', 'CH2_sp3', 'COO', 'OCO'].
		ERROR: Missing dihedral with ids (4, 5, 7, 8) and types ['CH2_sp3', 'COO', 'OCO', 'CH2_OC'].
		ERROR: Missing dihedral with ids (5, 7, 8, 9) and types ['COO', 'OCO', 'CH2_OC', 'CH2_sp3'].
		ERROR: Missing dihedral with ids (6, 5, 7, 8) and types ['OC', 'COO', 'OCO', 'CH2_OC'].
			 8 dihedral of 21
		ERROR: Missing dihedral with ids (7, 8, 9, 10) and types ['OCO', 'CH2_OC', 'CH2_sp3', 'CH2_sp3'].
		ERROR: Missing dihedral with ids (8, 9, 10, 11) and types ['CH2_OC', 'CH2_sp3', 'CH2_sp3', 'COO'].
		ERROR: Missing dihedral with ids (9, 10, 11, 12) and types ['CH2_sp3', 'CH2_sp3', 'COO', 'OC'].
		ERROR: Missing dihedral with ids (9, 10, 11, 13) and types ['CH2_sp3', 'CH2_sp3', 'COO', 'OCO'].
			 12 dihedral of 21
		ERROR: Missing dihedral with ids (10, 11, 13, 14) and types ['CH2_sp3', 'COO', 'OCO', 'CH2_OC'].
		ERROR: Missing dihedral with ids (11, 13, 14, 15) and types ['COO', 'OCO', 'CH2_OC', 'CH2_sp3'].
		ERROR: Missing dihedral with ids (12, 11, 13, 14) and types ['OC', 'COO', 'OCO', 'CH2_OC'].
		ERROR: Missing dihedral with ids (13, 14, 15, 16) and types ['OCO', 'CH2_OC', 'CH2_sp3', 'CH2_sp3'].
			 16 dihedral of 21
		ERROR: Missing dihedral with ids (14, 15, 16, 17) and types ['CH2_OC', 'CH2_sp3', 'CH2_sp3', 'COO'].
		ERROR: Missing dihedral with ids (15, 16, 17, 18) and types ['CH2_sp3', 'CH2_sp3', 'COO', 'OC'].
		ERROR: Missing dihedral with ids (15, 16, 17, 19) and types ['CH2_sp3', 'CH2_sp3', 'COO', 'OCO'].
		ERROR: Missing dihedral with ids (16, 17, 19, 20) and types ['CH2_sp3', 'COO', 'OCO', 'CH3_OC'].
			 20 dihedral of 21

		*********************************************************************************************************************************
		ERROR: Parameters have not been assigned to all proper dihedrals.
		Total system dihedrals: 21, Parameterized dihedrals: 3.
		Note that if your system contains torsions of Ryckaert-Bellemans functional form, all of these torsions are processed as propers.
		Molecule cannot be typed!!!!!. Exiting....
		*********************************************************************************************************************************

```

This information can be used to insert (add) new types in the xml file.
