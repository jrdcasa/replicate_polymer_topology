import unittest
import datetime
import os
import shutil
import replicate_polymer.replicate_polymer as replicate_polymer
from replicate_polymer_lib.remove_hydrogens import remove_hydrogens
from replicate_polymer_lib.setup_chiral_impropers import setup_chiral_impropers

class TestRemoveHydrogen(unittest.TestCase):

    # ===============================================================
    @classmethod
    def setUpClass(cls):

        cls.filelog = "./test/test03_setup_chiral_impropers.log"
        cls.log = replicate_polymer.init_logger("Output", fileoutput=cls.filelog, append=False, inscreen=False)
        m = "\n\t***************** START Impropers TEST *****************"
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
    def test_01_improper01(self):

        m = "\tTest_01. Automatic impropers in an United Atom model."
        print(m) if self.log is None else self.log.info(m)
        m2 = "\t"+len(m)*"="+"\n"
        print(m2) if self.log is None else self.log.info(m2)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        fnamepdb1 = "./data/C08_02br01_UA.pdb"
        nimpr, imp_lines = setup_chiral_impropers(fnamepdb1, improper_filename=None)

        self.assertEqual(nimpr, 2)
        string = "; improper angles\n;  ai    aj    ak    al funct  c0  c1  c2  c3  c4  c5\n" \
                 "       3        2        4        9  2    30.500 0.5178E+03\n" \
                 "       6        5        7       10  2    30.500 0.5178E+03\n"
        self.assertEqual(imp_lines, string)

        m = "\tPASSED\n"
        print(m) if self.log is None else self.log.info(m)

    # ===============================================================
    def test_02_improper02(self):

        m = "\tTest_02. Initial model PDB with Hs. Remove Hs and Trappe force field."
        print(m) if self.log is None else self.log.info(m)
        m2 = "\t"+len(m)*"="+"\n"
        print(m2) if self.log is None else self.log.info(m2)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        fnamepdb1 = "./data/C08_02br01.pdb"

        # Remove Hs ========
        filenamepdb_new, filenamegro_new = remove_hydrogens(fnamepdb1)
        # Guess impropers ========
        nimpr, imp_lines = setup_chiral_impropers(filenamepdb_new, improper_filename=None)

        self.assertEqual(nimpr, 2)
        string = "; improper angles\n;  ai    aj    ak    al funct  c0  c1  c2  c3  c4  c5\n" \
                 "       3        2        4        9  2    30.500 0.5178E+03\n" \
                 "       6        5        7       10  2    30.500 0.5178E+03\n"
        self.assertEqual(imp_lines, string)

        self.movefile(filenamepdb_new)
        self.movefile(filenamegro_new)

        m = "\tPASSED\n"
        print(m) if self.log is None else self.log.info(m)

    # ===============================================================
    def test_03_improper03(self):

        m = "\tTest_03. "
        print(m) if self.log is None else self.log.info(m)
        m2 = "\t"+len(m)*"="+"\n"
        print(m2) if self.log is None else self.log.info(m2)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        fnamepdb1 = "./data/C08_02br01.pdb"

        # Guess impropers ========
        nimpr, imp_lines = setup_chiral_impropers(fnamepdb1, improper_filename=None)

        self.assertEqual(nimpr, 2)
        string = "; improper angles\n;  ai    aj    ak    al funct  c0  c1  c2  c3  c4  c5\n" \
                 "       8        5       10       25  2    30.500 0.5178E+03\n" \
                 "      16       13       18       29  2    30.500 0.5178E+03\n"
        self.assertEqual(imp_lines, string)

        m = "\tPASSED\n"
        print(m) if self.log is None else self.log.info(m)

    # ===============================================================
    def test_04_improper04(self):

        m = "\tTest_04. A P4HB oligomer."
        print(m) if self.log is None else self.log.info(m)
        m2 = "\t"+len(m)*"="+"\n"
        print(m2) if self.log is None else self.log.info(m2)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        fnamepdb1 = "./data/00-Inital_trimerP4HB.pdb"

        # Guess impropers ========
        nimpr, imp_lines = setup_chiral_impropers(fnamepdb1, improper_filename=None)

        self.assertEqual(nimpr, 0)

        m = "\tPASSED\n"
        print(m) if self.log is None else self.log.info(m)

    # ===============================================================
    @classmethod
    def tearDownClass(cls):

        m = "\n\t***************** END Impropers TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tFinishing: \t {}\n".format(now))