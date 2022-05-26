import unittest
import datetime
import os
import shutil
import replicate_polymer.replicate_polymer as replicate_polymer
from replicate_polymer_lib.remove_hydrogens import remove_hydrogens
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import mdtraj
warnings.filterwarnings("ignore")

class TestRemoveHydrogen(unittest.TestCase):

    # ===============================================================
    @classmethod
    def setUpClass(cls):

        cls.filelog = "./test/test02_remove_hydrogen.log"
        cls.log = replicate_polymer.init_logger("Output", fileoutput=cls.filelog, append=False, inscreen=False)
        m = "\n\t***************** START Remove_Hydrogen TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tStarting: \t {}\n".format(now))

    # ===============================================================
    @staticmethod
    def movefile(filenamepdb, extrafolder="test"):

        # Check if file already exists and move the file
        fullpath = os.path.abspath(filenamepdb)
        destdir = os.path.join(os.path.dirname(fullpath), extrafolder)
        filename = os.path.basename(fullpath)
        destfullpath = os.path.join(destdir,filename)
        if os.path.isfile(destfullpath):
            os.remove(destfullpath)

        shutil.move(fullpath, destdir)

    # ===============================================================
    def test_01_remove_hydrogen_alkane(self):

        m = "\tTest_01: CONECT section exists in the PDB file."
        print(m) if self.log is None else self.log.info(m)
        m2 = "\t"+len(m)*"="+"\n"
        print(m2) if self.log is None else self.log.info(m2)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        # CONECT exists in the pdb file
        fnamepdb1 = "./data/103191_noctane.order_cryst.pdb"
        fnamepdb_noh_1, fnamegro_noh_1 = remove_hydrogens(fnamepdb1, removeallh=False)
        fnamepdb_noh_2, fnamegro_noh_2 = remove_hydrogens(fnamepdb1, removeallh=True)

        mdtrajpdb1 = mdtraj.load_pdb(fnamepdb_noh_1)
        mdtrajgro1 = mdtraj.formats.GroTrajectoryFile(fnamegro_noh_1)
        mdtrajpdb2 = mdtraj.load_pdb(fnamepdb_noh_2)
        mdtrajgro2 = mdtraj.formats.GroTrajectoryFile(fnamegro_noh_2)

        # Check natoms
        self.assertEqual(mdtrajpdb1.n_atoms, 8)
        self.assertEqual(mdtrajgro1.n_atoms, 8)
        m = "\tAll atoms are equal to 8 with removeallh=False"
        print(m) if self.log is None else self.log.info(m)

        self.assertEqual(mdtrajpdb2.n_atoms, 8)
        self.assertEqual(mdtrajgro2.n_atoms, 8)
        m = "\tAll atoms are equal to 8 with removeallh=True"
        print(m) if self.log is None else self.log.info(m)

        # Are there hydrogen atoms?
        nhs = len(mdtrajpdb1.topology.select("name H"))
        self.assertEqual(nhs, 0)
        m = "\tHydrogens are not present removeallh=False"
        print(m) if self.log is None else self.log.info(m)
        nhs = len(mdtrajpdb2.topology.select("name H"))
        self.assertEqual(nhs, 0)
        m = "\tHydrogens are not present removeallh=True"
        print(m) if self.log is None else self.log.info(m)

        self.movefile(fnamepdb_noh_1)
        self.movefile(fnamegro_noh_1)

        m = "\tPASSED\n"
        print(m) if self.log is None else self.log.info(m)

    # ===============================================================
    def test_02_remove_hydrogen_heavy_atom(self):

        m = "\tTest_02: CONECT section exists in the PDB file."
        print(m) if self.log is None else self.log.info(m)
        m2 = "\t"+len(m)*"="+"\n"
        print(m2) if self.log is None else self.log.info(m2)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        # CONECT exists in the pdb file
        fnamepdb1 = "./data/00-Inital_trimerP4HB.pdb"
        fnamepdb_noh_1, fnamegro_noh_1 = remove_hydrogens(fnamepdb1, removeallh=False)
        mdtrajpdb1 = mdtraj.load_pdb(fnamepdb_noh_1)
        mdtrajgro1 = mdtraj.formats.GroTrajectoryFile(fnamegro_noh_1)

        fnamepdb_noh_2, fnamegro_noh_2 = remove_hydrogens(fnamepdb1, removeallh=True)
        mdtrajpdb2 = mdtraj.load_pdb(fnamepdb_noh_2)
        mdtrajgro2 = mdtraj.formats.GroTrajectoryFile(fnamegro_noh_2)

        # Check natoms
        self.assertEqual(mdtrajpdb1.n_atoms, 21)
        self.assertEqual(mdtrajgro1.n_atoms, 21)
        m = "\tAll atoms are equal to 20 with removeallh=False"
        print(m) if self.log is None else self.log.info(m)

        self.assertEqual(mdtrajpdb2.n_atoms, 20)
        self.assertEqual(mdtrajgro2.n_atoms, 20)
        m = "\tAll atoms are equal to 8 with removeallh=True"
        print(m) if self.log is None else self.log.info(m)

        # Are there hydrogen atoms?
        nhs = len(mdtrajpdb1.topology.select("name H"))
        self.assertEqual(nhs, 1)
        m = "\t1 Hydrogen is present removeallh=False"
        print(m) if self.log is None else self.log.info(m)
        nhs = len(mdtrajpdb2.topology.select("name H"))
        self.assertEqual(nhs, 0)
        m = "\tHydrogens are not present removeallh=True"
        print(m) if self.log is None else self.log.info(m)

        self.movefile(fnamepdb_noh_1)
        self.movefile(fnamegro_noh_1)

        m = "\tPASSED\n"
        print(m) if self.log is None else self.log.info(m)


    # ===============================================================
    @classmethod
    def tearDownClass(cls):

        m = "\n\t***************** END Remove_Hydrogen TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tFinishing: \t {}\n".format(now))