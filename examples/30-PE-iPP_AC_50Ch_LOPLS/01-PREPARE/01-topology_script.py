import datetime
from topology.readmol.readxsdformat import ReadXsdFormat
from topology.readmol.readpdbformat import ReadPdbFormat
import utils

filename = "./iPP-PE_40mon_25-25Ch.xsd"
pattern = "iPP-PE_40mon_25-25Ch_order"
headinfo_file = "./iPP-PE_40mon_25-25Ch_headtail.dat"
residueinfo_file = "./iPP-PE_40mon_25-25Ch_residues.dat"

# Logger
filelog = pattern+".log"
log = utils.init_logger("Output", fileoutput=filelog, append=False, inscreen=True)
m = "\n\t***************** BUTANE 1 chain *****************"
print(m) if log is None else log.info(m)
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
log.info("\t\tStarting: \t {}\n".format(now))

# Create the xsd object
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Reading {}...({})".format(filename, now)
obj = ReadXsdFormat(filename)
print(m) if log is None else log.info(m)

# Write pdb file
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Writing pdb file...({})".format(now)
print(m) if log is None else log.info(m)
filenamepdb = "{}.pdb".format(pattern)
obj.write_pdb(filename_pdb=filenamepdb, separate_chains=False)

# Renumber PDB
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Read PDB and renumber pdb file...({})".format(now)
print(m) if log is None else log.info(m)
pdb = ReadPdbFormat(filenamepdb)
head_atoms, tail_atoms = pdb.read_head_tail_info(headinfo_file)
test = pdb.write_renumber_pdb(head_idx_atom=head_atoms, tail_idx_atom=tail_atoms)

# Assign residues PDB
filenamepdb=pattern+"_renumber.pdb"
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Read PDB and assign residues...({})".format(now)
print(m) if log is None else log.info(m)
pdb_new = ReadPdbFormat(filenamepdb)
pdb_new.assign_residues_chains(residueinfo_file)


# Logger
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
log.info("\n\t\tFinishing: \t {}\n".format(now))
m = "\t============== END   ==============================="
print(m) if log is None else log.info(m)


