import unittest
import datetime
#import replicate_polymer_lib as replicate_polymer
import replicate_polymer_lib.logger
import replicate_polymer_lib.replicate_pdb
import numpy as np


class TopocordpolymerTests(unittest.TestCase):

    # ##################################################################################################################
    @classmethod
    def setUpClass(cls):

        cls.filelog = "./test05/test01_replicate_polymer.log"
        cls.log = replicate_polymer_lib.logger.init_logger("Output", fileoutput=cls.filelog, append=False, inscreen=False)
        m = "\n\t***************** START Gecos rdkit TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tStarting: \t {}\n".format(now))

    # ##################################################################################################################
    def test_01_simulation_box(self):

        """
        Test if the simulation box is correctly handle

        """

        m = "\tTest_01: Test box information."
        print(m) if self.log is None else self.log.info(m)
        m2 = "\t"+len(m)*"="+"\n"
        print(m2) if self.log is None else self.log.info(m2)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        # Test 01
        f1 = "./data/103191_noctane.order_nocryst.pdb"
        m = "\t\tPDB without CRYST section ({}): ".format(f1)
        trj_single, _, _, _ = replicate_polymer_lib.replicate_pdb.replicate_pdb(f1, [1, 1, 1])
        bl = np.array([0.4194, 0.2949, 0.9888])
        al = np.array([90., 90., 90.])
        np.testing.assert_almost_equal(trj_single.unitcell_lengths[0, :], bl)
        np.testing.assert_almost_equal(trj_single.unitcell_angles[0, :], al)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)

        # Test 02
        m = "\t\tPDB without CRYST section ({}) boxlength through interface: ".format(f1)
        trj_single, _, _, _ = replicate_polymer_lib.replicate_pdb.replicate_pdb(f1, [1, 1, 1], boxlength=[10.1, 10.2, 10.3])
        bl = np.array([10.1, 10.2, 10.3], dtype=np.float32)
        al = np.array([90., 90., 90.], dtype=np.float32)
        np.testing.assert_almost_equal(trj_single.unitcell_lengths[0, :], bl)
        np.testing.assert_almost_equal(trj_single.unitcell_angles[0, :], al)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)

        # Test 03
        m = "\t\tPDB without CRYST section ({}) boxangles through interface: ".format(f1)
        trj_single, _, _, _ = replicate_polymer_lib.replicate_pdb.replicate_pdb(f1, [1, 1, 1], boxangles=[88.0, 78.9, 90.0])
        bl = np.array([0.4194, 0.2949, 0.9888], dtype=np.float32)
        al = np.array([88.0, 78.9, 90.0], dtype=np.float32)
        np.testing.assert_almost_equal(trj_single.unitcell_lengths[0, :], bl)
        np.testing.assert_almost_equal(trj_single.unitcell_angles[0, :], al)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)

        # Test 04
        m = "\t\tPDB without CRYST section ({}) boxlengths and boxangles through interface: ".format(f1)
        trj_single, _, _, _ = replicate_polymer_lib.replicate_pdb.replicate_pdb(f1, [1, 1, 1], boxlength=[10.1, 10.2, 10.3],
                                                        boxangles=[88.0, 78.9, 90.0])
        bl = np.array([10.1, 10.2, 10.3], dtype=np.float32)
        al = np.array([88.0, 78.9, 90.0], dtype=np.float32)
        np.testing.assert_almost_equal(trj_single.unitcell_lengths[0, :], bl)
        np.testing.assert_almost_equal(trj_single.unitcell_angles[0, :], al)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)

        # Test 05
        f2 = "./data/103191_noctane.order_cryst.pdb"
        m = "\t\tPDB with CRYST section ({}): ".format(f2)
        trj_single, _, _, _ = replicate_polymer_lib.replicate_pdb.replicate_pdb(f2, [1, 1, 1])
        bl = np.array([0.4123, 0.4686, 1.0974], dtype=np.float32)
        al = np.array([85.06, 83.72, 75.10], dtype=np.float32)
        np.testing.assert_almost_equal(trj_single.unitcell_lengths[0, :], bl)
        np.testing.assert_almost_equal(trj_single.unitcell_angles[0, :], al)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)

        # Test 06
        m = "\t\tPDB with CRYST section ({}) boxlength through interface: ".format(f2)
        trj_single, _, _, _ = replicate_polymer_lib.replicate_pdb.replicate_pdb(f2, [1, 1, 1], boxlength=[10.1, 10.2, 10.3])
        bl = np.array([10.1, 10.2, 10.3], dtype=np.float32)
        al = np.array([85.06, 83.72, 75.10], dtype=np.float32)
        np.testing.assert_almost_equal(trj_single.unitcell_lengths[0, :], bl)
        np.testing.assert_almost_equal(trj_single.unitcell_angles[0, :], al)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)

        # Test 07
        m = "\t\tPDB with CRYST section ({}) boxangles through interface: ".format(f2)
        trj_single, _, _, _ = replicate_polymer_lib.replicate_pdb.replicate_pdb(f2, [1, 1, 1], boxangles=[88.0, 78.9, 90.0])
        bl = np.array([0.4123, 0.4686, 1.0974], dtype=np.float32)
        al = np.array([88.0, 78.9, 90.0], dtype=np.float32)
        np.testing.assert_almost_equal(trj_single.unitcell_lengths[0, :], bl)
        np.testing.assert_almost_equal(trj_single.unitcell_angles[0, :], al)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)

        # Test 08
        m = "\t\tPDB with CRYST section ({}) boxlength and boxangles through interface: ".format(f2)
        trj_single, _, _, _ = replicate_polymer_lib.replicate_pdb.replicate_pdb(f2, [1, 1, 1], boxlength=[10.1, 10.2, 10.3],
                                                        boxangles=[88.0, 78.9, 90.0])
        bl = np.array([10.1, 10.2, 10.3], dtype=np.float32)
        al = np.array([88.0, 78.9, 90.0], dtype=np.float32)
        np.testing.assert_almost_equal(trj_single.unitcell_lengths[0, :], bl)
        np.testing.assert_almost_equal(trj_single.unitcell_angles[0, :], al)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)
    #
    # ##################################################################################################################
    def test_02_shift_positions(self):

        """
        Test if the coordinates in the replicas are correctly calculated.
        """

        m = "\tTest_02: Test atom shifts in a replica."
        print(m) if self.log is None else self.log.info(m)
        m2 = "\t"+len(m)*"="+"\n"
        print(m2) if self.log is None else self.log.info(m2)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        # Test 01
        f2 = "./data/103191_noctane.order_cryst.pdb"
        m = "\t\tReplica [0,0,1]. PDB with CRYST section ({}): ".format(f2)
        trj_single, _, _, _ = replicate_polymer_lib.replicate_pdb.replicate_pdb(f2, [1, 1, 1])
        c0 = trj_single.xyz[0, 0, :]
        v1 = np.array([0, 0, 1])
        newc = c0 + np.dot(v1, trj_single.unitcell_vectors)
        np.testing.assert_almost_equal(newc[0], np.array([0.3067, 0.4108, 1.2180]), decimal=4)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)

        # Test 02
        m = "\t\tReplica [1,0,1]. PDB with CRYST section ({}): ".format(f2)
        v1 = np.array([1, 0, 1])
        newc = c0 + np.dot(v1, trj_single.unitcell_vectors)
        np.testing.assert_almost_equal(newc[0], np.array([0.7190, 0.4108, 1.2180]), decimal=4)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)

        # Test 03
        m = "\t\tReplica [1,0,1]. PDB with CRYST section ({}): ".format(f2)
        v1 = np.array([3, 2, 1])
        newc = c0 + np.dot(v1, trj_single.unitcell_vectors)
        np.testing.assert_almost_equal(newc[0], np.array([1.7846, 1.3165, 1.21800]), decimal=4)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)

        # Test 04
        m = "\t\tReplica [0,0,0]. PDB with CRYST section ({}): ".format(f2)
        v1 = np.array([0, 0, 0])
        newc = c0 + np.dot(v1, trj_single.unitcell_vectors)
        np.testing.assert_almost_equal(newc[0], np.array([0.1866, 0.3450, 0.12920]), decimal=4)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)

        # Test 05
        f1 = "./data/103191_noctane.order_nocryst.pdb"
        m = "\t\tPDB without CRYST section ({}): ".format(f1)
        trj_single, _, _, _ = replicate_polymer_lib.replicate_pdb.replicate_pdb(f1, [1, 1, 1], boxlength=[0.4123, 0.4686, 1.0974],
                                                        boxangles=[85.06, 83.72, 75.10])
        m = "\t\tReplica [1,0,1]. PDB without CRYST section ({}): ".format(f2)
        v1 = np.array([1, 0, 1])
        newc = c0 + np.dot(v1, trj_single.unitcell_vectors)
        np.testing.assert_almost_equal(newc[0], np.array([0.7190, 0.4108, 1.2180]), decimal=4)
        m += "PASSED\n"
        print(m) if self.log is None else self.log.info(m)

        print(c0, trj_single.unitcell_vectors, newc)
        pass

    # ##################################################################################################################
    def test_03_inputmore100Katoms_noconnect(self):

        f1 = "./data/103191_noctane.order_cryst_noconect.pdb"
        m = "\t\tPDB without CRYST section ({}): ".format(f1)
        trj_single, _, _, _ = replicate_polymer_lib.replicate_pdb.replicate_pdb(f1, [1, 1, 1])



    # ##################################################################################################################
    @classmethod
    def tearDownClass(cls):

        m = "\n\t***************** END Gecos rdkit TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tFinishing: \t {}\n".format(now))