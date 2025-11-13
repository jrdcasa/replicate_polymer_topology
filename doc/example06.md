# Example 6: A system with two different molecules

**The files to run this example are available in the [examples](../examples/06-pdb_two_diferent_molecules/) folder inside the distribution**

```bash
cp ../04-Trimer_P4HB_OPLSAA/TrimerP4HB_residues.pdb .

replicate_polymer -p two_molecules.pdb -f ../../forcefields/oplsaa.xml --images 2 2 2 -e lammps --boxlength 2.3 2.3 2.3
```