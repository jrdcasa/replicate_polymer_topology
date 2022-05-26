import unittest
import datetime
import io
import os
import shutil
import replicate_polymer.replicate_polymer as replicate_polymer
from replicate_polymer_lib.check_connect_pdb import check_conect_pdb


class TestCheckConnectPdb(unittest.TestCase):

    # ===============================================================
    @classmethod
    def setUpClass(cls):

        cls.filelog = "./test/test01_check_connect_pdb.log"
        cls.log = replicate_polymer.init_logger("Output", fileoutput=cls.filelog, append=False, inscreen=False)
        m = "\n\t***************** START Check_connect_pdb TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tStarting: \t {}\n".format(now))

    # ===============================================================
    def test_01_pdb_connect(self):

        m = "\tTest_01: CONECT section exists in the PDB file."
        print(m) if self.log is None else self.log.info(m)
        m2 = "\t"+len(m)*"="+"\n"
        print(m2) if self.log is None else self.log.info(m2)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        # CONECT exists in the pdb file
        fnamepdb = "./data/103191_noctane.order_cryst.pdb"
        filenamepdb = check_conect_pdb(fnamepdb)

        self.assertListEqual(
            list(io.open(fnamepdb)),
            list(io.open(filenamepdb)))
        m = "\tBoth files are the same\n"
        m += "\tPASSED\n"
        print(m) if self.log is None else self.log.info(m)

    # ===============================================================
    def test_02_pdb_noconnect(self):

        m = "\tTest_02: CONECT section does not exist in the PDB file."
        print(m) if self.log is None else self.log.info(m)
        m2 = "\t"+len(m)*"="+"\n"
        print(m2) if self.log is None else self.log.info(m2)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        # CONECT does not exist in the pdb file
        fnamepdb = "./data/103191_noctane.order_cryst_noconect.pdb"
        filenamepdb = check_conect_pdb(fnamepdb)

        # Check if file already exists and move the file
        fullpath = os.path.abspath(filenamepdb)
        destdir = os.path.join(os.path.dirname(fullpath), "test")
        filename = os.path.basename(fullpath)
        destfullpath = os.path.join(destdir,filename)
        if os.path.isfile(destfullpath):
            os.remove(destfullpath)
        with open(filenamepdb, 'r') as f:
            lines = f.readlines()
            ncon = sum('CONECT' in s for s in lines)

        shutil.move(fullpath, destdir)

        self.assertGreater(ncon, 0)

        m = "\tCONECT entries have been created.\n"
        m += "\tPASSED\n"
        print(m) if self.log is None else self.log.info(m)

    # ===============================================================
    def test_03_pdb_noconnect_noresname(self):

        m = "\tTest_03: CONECT section does not exist in the PDB file.\n"
        m += "\t         RESNAME in PDB are not defined. This raises an error\n"
        m += "\t         in parmed package used by check_conect_pdb"
        print(m) if self.log is None else self.log.info(m)
        m2 = "\t"+len(m)*"="+"\n"
        print(m2) if self.log is None else self.log.info(m2)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        # CONECT does not exist in the pdb file
        fnamepdb = "./data/103191_noctane.order_cryst_noconect_noresname.pdb"
        filenamepdb = check_conect_pdb(fnamepdb)

        self.assertIsNone(filenamepdb)

        m = "\tBad format for PDB. No residues name in PDB file\n"
        m += "\tNone object is capture!!!\n"
        m += "\tPASSED\n"
        print(m) if self.log is None else self.log.info(m)


    # ===============================================================
    @classmethod
    def tearDownClass(cls):

        m = "\n\t***************** END Check_connect_pdb TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tFinishing: \t {}\n".format(now))