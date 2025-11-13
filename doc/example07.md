# Example 7: Polymer blends

**The files to run this example are available in the [examples](../examples/07-blend_xsd/) folder inside the distribution**

Five inputs files are required to run the program: 

* A pdb file containing the molecule to be replicated
* A xml file with force field parameters. The basic format is described [here](https://foyer.mosdef.org/en/stable/topic_guides/smarts.html)
* To order the atoms and type the residues, two files more are required: a file with information of the head and tail atoms [blend_headtail.dat](../examples/07-blend_xsd/blend_headtail.dat) and other with the residue information [blend_residues.dat](../examples/07-blend_xsd/blend_residues.dat).
* A script to order atoms and assign residues [01-topology_script.py](../examples/07-blend_xsd/01-topology_script.py)


```bash
source ~/Programacion/sandboxes/sandbox_common/bin/activate

python 01-topology_script.py 


```