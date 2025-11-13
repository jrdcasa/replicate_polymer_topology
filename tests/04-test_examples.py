import unittest
import datetime
import subprocess
import replicate_polymer.replicate_polymer as replicate_polymer

class TestExamples(unittest.TestCase):

    # ===============================================================
    @classmethod
    def setUpClass(cls):

        cls.filelog = "./test/test04_test_examples.log"
        cls.log = replicate_polymer.init_logger("Output", fileoutput=cls.filelog, append=False, inscreen=False)
        m = "\n\t***************** START Check_connect_pdb TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tStarting: \t {}\n".format(now))

    # ===============================================================
    def test_01_octane(self):


        try:
            result = subprocess.run(['python', 'replicate_polymer.py'], check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as error:
            print(error.stdout)
            print(error.stderr)
            raise error