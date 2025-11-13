import unittest
import datetime
import subprocess
import os
import shutil
import create_library_monomer.create_library_monomer as create_library_monomer
import replicate_polymer.replicate_polymer as replicate_polymer


class TestPolyply(unittest.TestCase):

    # ===============================================================
    @classmethod
    def setUpClass(cls):

        dir_target = "./test05"
        if os.path.exists("./test05"):
            shutil.rmtree(dir_target)
        os.mkdir(dir_target)
        os.chdir(dir_target)

        cls.filelog = "test05_test_examples.log"
        cls.log = replicate_polymer.init_logger("Output", fileoutput=cls.filelog, append=False, inscreen=False)
        m = "\n\t***************** START Test generate systems with polyply TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tStarting: \t {}\n".format(now))

    # ===============================================================
    def test01_PE_OPLS(self):

        m = "\n\t ============================== START 05 OPLS PE ==============================\n"
        print(m) if self.log is None else self.log.info(m)

        exe = "../../create_library_monomer/create_library_monomer.py"

        bash_command = "python {} -p ../00-MONOMERS_MODEL/PE_model_residues.pdb " \
                       "-f /home/jramos/Programacion/GITHUB_REPO_DIR/replicate_polymer_topology/forcefields/oplsaa.xml " \
                       "--pattern ff_PE_opls_cJ".format(exe)
        process = subprocess.Popen(bash_command.split(),
                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\tTest 05: Output from PE_OPLS.({})".format(now)
        m += "\n\t\t{}\n".format(len(m) * "*")
        m += output.decode()
        m += error.decode()
        print(m) if self.log is None else self.log.info(m)

        m = "\n\t ============================== END TEST 05 OPLS PE ==============================\n"
        print(m) if self.log is None else self.log.info(m)


    # ===============================================================
    @classmethod
    def tearDownClass(cls):

        m = "\n\t***************** END Test generate systems with polyply END *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tFinishing: \t {}\n".format(now))