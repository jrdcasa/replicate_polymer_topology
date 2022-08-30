# How to prepare from gabedit

1. Sketch with Gabedit two molecules: A linear PE with 5 monomers and a SCB-PE with 5 monomers and one C4 branch ([pdb](./00-INPUT_GABEDIT/PE-PESCB_gabedit.pdb)). Edit the occupancy and beta fields to have a value of 0.00.

2. Apply topology tools.

    * Create a file with the head and tail info. Use VMD to find the indices [Toxwaerd_twokinds_headtail.dat](./00-INPUT_GABEDIT/Toxwaerd_twokinds_headtail.dat).

        * Mol A: Head: 27 Tail: 0
        * Mol B: Head: 32 Tail: 58

    * Create a file with the residue info. Use VMD to find the residues [Toxwaerd_twokinds_residues.dat](./00-INPUT_GABEDIT/Toxwaerd_twokinds_residues.dat).

        * Mol A: Three type of residues, PES (terminal CH3-CH2-), PEE (middle -CH2-CH2-) and HEX (branched residue, -CH2-CH(C4H9)-).
        * Mol B: Two type of residues, PES (terminal CH3-CH2-) and PEE (middle -CH2-CH2-).  

    * cd [01-PREPARE](./01-PREPARE/)
    * python [01-topology_script.py](./01-PREPARE/01-topology_script.py)
    * The ouput files are [Toxwaerd_twokinds_residues.gro](./01-PREPARE/Toxwaerd_twokinds_residues.gro) and [Toxwaerd_twokinds_residues.pdb](./01-PREPARE/Toxwaerd_twokinds_residues.pdb).
    * Check the new files with VMD.
    * CAUTION: The name of the output files have the same root adding "_residues" tag that the residue file.