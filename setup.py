import sys
import logging
import subprocess
import os
import tempfile
import shutil
import logging
import numpy
from datetime import datetime
from setuptools import setup, Extension
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler


# Formatter for the logger
class CustomFormatter(logging.Formatter):

    """Logging Formatter to add colors and count warning / errors"""
    FORMATS = {
        logging.ERROR: "\n\tERROR: %(asctime)s: %(msg)s",
        logging.WARNING: "\n\tWARNING: %(msg)s",
        logging.DEBUG: "%(asctime)s: %(msg)s",
        "DEFAULT": "%(msg)s",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        date_fmt = '%d-%m-%Y %d %H:%M:%S'
        formatter = logging.Formatter(log_fmt, date_fmt)
        return formatter.format(record)


# Install packages from pip ==============================================================
def install_with_pip(pack, vers=None, log=None):

    # sys.executable gives the path of the python interpreter
    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    if vers is None:
        m = "{}: ** REPLICATE: Installing {}".format(now, pack)
        print(m) if log is None else log.info(m)
        # subprocess.call([sys.executable, "-m", "pip", "install", "{0}".format(pack)])
        p = subprocess.Popen([sys.executable, "-m", "pip", "install", "{0}".format(pack)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.communicate()
    else:
        m = "{}: ** REPLICATE: Installing {}=={}".format(now, pack, vers)
        print(m) if log is None else log.info(m)
        # subprocess.call([sys.executable, "-m", "pip", "install", "{0}=={1}".format(pack, vers), " &>install.log"])
        p = subprocess.Popen([sys.executable, "-m", "pip", "install", "{0}=={1}".format(pack, vers)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.communicate()

# ================================================================================================
def install_topology_library(log=None, namepkg=None):
    """
    Installing the python topology library if is not present in the python environment.
    """

    import git

    giturl = "https://github.com/jrdcasa/topology.git"
    install_dir = 'topology'
    install_dir = './thirdparty/topology'

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    m = "\n\t\t INSTALLING TOPOLOGY LIBRARY INTO THE PYTHON ENVIRONMENT FROM GITHUB\n\n"

    if os.path.isdir(os.path.join(os.getcwd(), install_dir)):
        m += "{}: ** {}: topology library is already installed in your system. {}".format(now, namepkg, giturl)
        print(m) if log is None else log.info(m)
    else:
        m += "{}: ** {}: topology library is not installed in your system\n".format(now, namepkg)
        m += "{}: ** {}: Installing from git... {}\n".format(now, namepkg, giturl)
        print(m) if log is None else log.info(m)

        fullpath_install = os.path.abspath(install_dir)

        # Look at thirdparty directory
        if os.path.isdir(fullpath_install):
            shutil.rmtree(fullpath_install)
        os.makedirs(fullpath_install)

        try:
            git.Repo.clone_from(giturl, fullpath_install)
        except git.GitCommandError as e:
            m = "================= ERROR INSTALL ================\n"
            m += "** {}: The github repository for topology is not valid or not exists.!!!\n".format(namepkg)
            m += "** {}: giturl     : {}\n".format(namepkg, giturl)
            m += "** {}: install_dir: {}\n".format(namepkg, fullpath_install)
            m += "** {}: Topology library cannot be installed\n".format(namepkg)
            m += "** {}: The installation is aborted\n".format(namepkg)
            m += "** Error: {}\n".format(e.stderr)
            m += "\n================= ERROR INSTALL ================"
            print(m) if log is None else log.info(m)
            exit()

        os.chdir(fullpath_install)
        #print(os.getcwd())
        subprocess.call(["python", "setup.py", "install"])
        os.chdir("..")

    # Check if topology can be imported
    try:
        import topology
        import topology.readmol
        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        m = "topology installed in the python environment ({})".format(now)
        print(m) if log is None else log.info(m)
        return True
    except (ModuleNotFoundError, ImportError):
        m = "================= ERROR INSTALL ================\n"
        m += "** REPLICATE: Automatic installation of topology python is unsuccessful.\n"
        m += "** REPLICATE: Install the topology library manually using the repi:\n"
        m += "** REPLICATE: https://github.com/jrdcasa/topology.git\n"
        m += "** REPLICATE: and then try to re-install\n"
        m += "================= ERROR INSTALL ================"
        print(m) if log is None else log.info(m)
        exit()


# Install foyer ==============================================================
def install_foyer(log=None):

    import git

    giturl = 'https://:@github.com/jrdcasa/foyer_cj.git'
    # giturl = 'https://:@github.com/mosdef-hub/foyer.git'
    install_dir = './thirdparty/foyer'

    # Check if exists a distribution of foyer in the thirdparty directory
    # git clone https://github.com/jrdcasa/foyer_cj.git
    if os.path.isdir("thirdparty/foyer"):
        m = "{} is already cloned in {}".format(giturl, install_dir)
        print(m) if log is None else log.info(m)
        pass
    else:
        try:
            m = "Cloning from {} in {}".format(giturl, install_dir)
            print(m) if log is None else log.info(m)
            git.Repo.clone_from(giturl, install_dir)  # progress=CloneProgress())
            m = "Repository cloned\n"
            print(m) if log is None else log.info(m)
        except git.GitCommandError:
            if not os.path.isdir(install_dir):
                m = "================= ERROR INSTALL ================\n"
                m += "** REPLICATE: The github repository for foyer is not valid or not exist.!!!\n"
                m += "** REPLICATE: giturl     : {}\n".format(giturl)
                m += "** REPLICATE: install_dir: {}\n".format(install_dir)
                m += "** REPLICATE: foyer cannot be installed\n"
                m += "** REPLICATE: The installation is aborted\n"
                m += "================= ERROR INSTALL ================"
                print(m) if log is None else log.info(m)
                exit()
            else:
                pass

    # Install into the python environment
    wk = os.getcwd()
    os.chdir(install_dir)
    subprocess.call(["python", "setup.py", "install"])
    m = "foyer installed in the python environment"
    print(m) if log is None else log.info(m)
    os.chdir(wk)


# Install openmm ==============================================================
def install_openmm(log=None):

    import git
    giturl = 'https://:@github.com/openmm/openmm.git'
    install_dir = './thirdparty/openmm'

    try:
        import openmm
        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        m = "openmm installed in the python environment ({})".format(now)
        print(m) if log is None else log.info(m)

        return True
    except (ModuleNotFoundError, ImportError):
        m = "Trying to install openmm python in an automatic way."
        print(m) if log is None else log.info(m)

    # Check if exists a distribution of openmm in the thirdparty directory
    # git clone https://:@github.com/openmm/openmm.git
    if os.path.isdir("thirdparty/openmm"):
        m = "{} is already cloned in {}".format(giturl, install_dir)
        print(m) if log is None else log.info(m)
        pass
    else:
        try:
            m = "Cloning from {} in {}".format(giturl, install_dir)
            print(m) if log is None else log.info(m)
            git.Repo.clone_from(giturl, install_dir)  # progress=CloneProgress())
            m = "Repository cloned\n"
            print(m) if log is None else log.info(m)
        except git.GitCommandError:
            if not os.path.isdir(install_dir):
                m = "================= ERROR INSTALL ================\n"
                m += "** REPLICATE: The github repository for openmm is not valid or not exist.!!!\n"
                m += "** REPLICATE: giturl     : {}\n".format(giturl)
                m += "** REPLICATE: install_dir: {}\n".format(install_dir)
                m += "** REPLICATE: foyer cannot be installed\n"
                m += "** REPLICATE: The installation is aborted\n"
                m += "================= ERROR INSTALL ================"
                print(m) if log is None else log.info(m)
                exit()
            else:
                pass

    # Install into the python environment
    wk = os.getcwd()
    os.chdir(install_dir)
    subprocess.call(["mkdir", "build"])
    subprocess.call(["mkdir", "openmm_bin"])
    os.chdir("./build")
    print(os.getcwd())
    m = "Making cmake...\n"
    print(m) if log is None else log.info(m)
    subprocess.call(["cmake", "..", "-DCMAKE_INSTALL_PREFIX={}".format("../openmm_bin")])
    m = "Make install (compile)...\n"
    print(m) if log is None else log.info(m)
    subprocess.call(["make", "install", "-j4"])
    m = "Making PythonInstall...\n"
    print(m) if log is None else log.info(m)
    subprocess.call(["make", "PythonInstall"])
    os.chdir(wk)

    try:
        import openmm
        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        m = "openmm installed in the python environment ({})".format(now)
        print(m) if log is None else log.info(m)
        return True
    except (ModuleNotFoundError, ImportError):
        m = "================= ERROR INSTALL ================\n"
        m += "** REPLICATE: Automatic installation of openmm python is unsuccessful.\n"
        m += "** REPLICATE: Install the openmm library manually following the instructions in:\n"
        m += "** REPLICATE: http://docs.openmm.org/latest/userguide/library/02_compiling.html\n"
        m += "** REPLICATE: and then try to re-install\n"
        m += "================= ERROR INSTALL ================"
        print(m) if log is None else log.info(m)
        exit()


# Install gmso ==============================================================
def install_gmso(log=None):

    import git

    giturl = 'https://:@github.com/jrdcasa/gmso_cj.git'
    # giturl = 'https://:@github.com/mosdef-hub/gmso.git'
    install_dir = './thirdparty/gmso'

    # Check if exists a distribution of foyer in the thirdparty directory
    # git clone https://github.com/jrdcasa/gmso_cj.git
    if os.path.isdir("thirdparty/gmso"):
        m = "{} is already cloned in {}".format(giturl, install_dir)
        print(m) if log is None else log.info(m)
        pass
    else:
        try:
            m = "Cloning from {} in {}".format(giturl, install_dir)
            print(m) if log is None else log.info(m)
            git.Repo.clone_from(giturl, install_dir)  # progress=CloneProgress())
            m = "Repository cloned\n"
            print(m) if log is None else log.info(m)
        except git.GitCommandError:
            if not os.path.isdir(install_dir):
                m = "================= ERROR INSTALL ================\n"
                m += "** REPLICATE: The github repository for gmso is not valid or not exist.!!!\n"
                m += "** REPLICATE: giturl     : {}\n".format(giturl)
                m += "** REPLICATE: install_dir: {}\n".format(install_dir)
                m += "** REPLICATE: gmso cannot be installed\n"
                m += "** REPLICATE: The installation is aborted\n"
                m += "================= ERROR INSTALL ================"
                print(m) if log is None else log.info(m)
                exit()
            else:
                pass

    # Install into the python environment
    wk = os.getcwd()
    os.chdir(install_dir)
    subprocess.call(["python", "setup.py", "install"])
    m = "gmso installed in the python environment"
    print(m) if log is None else log.info(m)
    os.chdir(wk)


# Install gmso ==============================================================
def install_mbuild(log=None):

    import git

    giturl = 'https://:@github.com/jrdcasa/mbuild.git'
    # giturl = 'https://:@github.com/mosdef-hub/mbuild.git'
    install_dir = './thirdparty/mbuild'

    # Check if exists a distribution of foyer in the thirdparty directory
    # git clone https://github.com/jrdcasa/mbuild.git
    if os.path.isdir("thirdparty/mbuild"):
        m = "{} is already cloned in {}".format(giturl, install_dir)
        print(m) if log is None else log.info(m)
        pass
    else:
        try:
            m = "Cloning from {} in {}".format(giturl, install_dir)
            print(m) if log is None else log.info(m)
            git.Repo.clone_from(giturl, install_dir)  # progress=CloneProgress())
            m = "Repository cloned\n"
            print(m) if log is None else log.info(m)
        except git.GitCommandError:
            if not os.path.isdir(install_dir):
                m = "================= ERROR INSTALL ================\n"
                m += "** REPLICATE: The github repository for mbuild is not valid or not exist.!!!\n"
                m += "** REPLICATE: giturl     : {}\n".format(giturl)
                m += "** REPLICATE: install_dir: {}\n".format(install_dir)
                m += "** REPLICATE: mbuild cannot be installed\n"
                m += "** REPLICATE: The installation is aborted\n"
                m += "================= ERROR INSTALL ================"
                print(m) if log is None else log.info(m)
                exit()
            else:
                pass

    # Install into the python environment
    wk = os.getcwd()
    os.chdir(install_dir)
    subprocess.call(["python", "setup.py", "install"])
    m = "mbuild installed in the python environment"
    print(m) if log is None else log.info(m)
    os.chdir(wk)


# Install intermol ==============================================================
def install_intermol(log=None):

    import git

    giturl = 'https://:@github.com/jrdcasa/intermol_cj.git'
    # giturl = 'https://:@github.com/shirtsgroup/InterMol'
    install_dir = './thirdparty/intermol'

    # Check if exists a distribution of foyer in the thirdparty directory
    # git clone https://github.com/jrdcasa/intermol_cj.git
    if os.path.isdir("thirdparty/intermol"):
        m = "{} is already cloned in {}".format(giturl, install_dir)
        print(m) if log is None else log.info(m)
        pass
    else:
        try:
            m = "Cloning from {} in {}".format(giturl, install_dir)
            print(m) if log is None else log.info(m)
            git.Repo.clone_from(giturl, install_dir)  # progress=CloneProgress())
            m = "Repository cloned\n"
            print(m) if log is None else log.info(m)
        except git.GitCommandError:
            if not os.path.isdir(install_dir):
                m = "================= ERROR INSTALL ================\n"
                m += "** REPLICATE: The github repository for intermol is not valid or not exist.!!!\n"
                m += "** REPLICATE: giturl     : {}\n".format(giturl)
                m += "** REPLICATE: install_dir: {}\n".format(install_dir)
                m += "** REPLICATE: foyer cannot be installed\n"
                m += "** REPLICATE: The installation is aborted\n"
                m += "================= ERROR INSTALL ================"
                print(m) if log is None else log.info(m)
                exit()
            else:
                pass

    # Install into the python environment
    wk = os.getcwd()
    os.chdir(install_dir)
    subprocess.call(["python", "setup.py", "install"])
    m = "intermol installed in the python environment"
    print(m) if log is None else log.info(m)
    os.chdir(wk)


# Remove warnings in intermol, parmed and mdtraj
def remove_warnings():

    import sysconfig
    path = sysconfig.get_paths()["purelib"]

    # mdtraj --> simtk.openmm warning
    wk = os.getcwd()
    os.chdir(os.path.join(path, 'mdtraj'))
    bashCommand = "grep -rl 'simtk.openmm' . |xargs sed -i 's/simtk\.openmm/openmm/g'"
    os.system(bashCommand)
   
    # parmed --> simtk.openmm warning
    os.chdir(os.path.join(path, 'parmed'))
    os.system(bashCommand)

    # intermol --> collections to collections.abc
    os.chdir(os.path.join(path, 'intermol-0.1.0-py3.8.egg'))
    bashCommand = "grep -rl 'collections.Set' intermol/orderedset.py |xargs sed -i 's/collections\.Set/collections\.abc\.Set/g'"
    os.system(bashCommand)

    os.chdir(os.path.join(wk))


# Disabling-output-when-compiling-with-distutil =================================================
def hasfunction(cc, funcname, include=None, extra_postargs=None):
    # From http://stackoverflow.com/questions/
    #            7018879/disabling-output-when-compiling-with-distutils
    tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
    devnull = oldstderr = None
    try:
        try:
            fname = os.path.join(tmpdir, 'funcname.c')
            with open(fname, 'w') as fout:
                if include is not None:
                    fout.write('#include {0!s}\n'.format(include))
                fout.write('int main(void) {\n')
                fout.write('    {0!s};\n'.format(funcname))
                fout.write('}\n')
            # Redirect stderr to /dev/null to hide any error messages
            # from the compiler.
            # This will have to be changed if we ever have to check
            # for a function on Windows.
            devnull = open('/dev/null', 'w')
            oldstderr = os.dup(sys.stderr.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())
            objects = cc.compile([fname], output_dir=tmpdir,
                                 extra_postargs=extra_postargs)
            cc.link_executable(objects, os.path.join(tmpdir, "a.out"))
        except Exception:
            return False
        return True
    finally:
        if oldstderr is not None:
            os.dup2(oldstderr, sys.stderr.fileno())
        if devnull is not None:
            devnull.close()
        shutil.rmtree(tmpdir)


# Does this compiler support OpenMP parallelization?""" ==============================================================
def detect_openmp():
    print("REPLICATE: Attempting to autodetect OpenMP support... ")
    compiler = new_compiler()
    customize_compiler(compiler)
    compiler.add_library('gomp')
    include = '<omp.h>'
    extra_postargs = ['-fopenmp']
    hasopenmp = hasfunction(compiler, 'omp_get_num_threads()', include=include,
                            extra_postargs=extra_postargs)
    if hasopenmp:
        print("REPLICATE: Compiler supports OpenMP")
    else:
        print("REPLICATE: Did not detect OpenMP support.")

    return hasopenmp


# Setup external extensions ==============================================================
def setup_external_extensions(debug_cflags=False, use_openmp=True):

    has_openmp = detect_openmp()

    # parallel_libraries = ['gomp'] if has_openmp and use_openmp else []
    mathlib = ['m']
    define_macros = []
    extra_compile_args = ['-std=c99', '-ffast-math', '-O3', '-funroll-loops', '-Wno-cpp']
    if debug_cflags:
        extra_compile_args.extend(['-Wall', '-pedantic'])
        define_macros.extend([('DEBUG', '1')])

    parallel_args = ['-fopenmp'] if has_openmp and use_openmp else []
    parallel_libraries = ['gomp'] if has_openmp and use_openmp else []
    parallel_macros = [('PARALLEL', None)] if has_openmp and use_openmp else []

    extensions_install = [
        Extension("ext_libc.c_distC", ["replicate_polymer_lib/ext_libc/c_distC.pyx"],
                  libraries=mathlib + parallel_libraries,
                  define_macros=define_macros + parallel_macros,
                  extra_compile_args=parallel_args + extra_compile_args,
                  extra_link_args=['-lgomp']),

    ]

    return extensions_install


# Main setup
if __name__ == '__main__':

    # Creating the logger to install.log file ===================================
    logger = logging.getLogger(name="INSTALL_LOG")
    logger.setLevel(logging.DEBUG)
    h1 = logging.FileHandler("install.log", 'w')
    h1.setFormatter(CustomFormatter())
    # Output also in the screen
    logger.addHandler(h1)
    f1 = logging.StreamHandler()
    f1.setFormatter(CustomFormatter())
    logger.addHandler(f1)

    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t Starting installation!!!! at {}\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)

    # Print sys path ===================================
    m1 = "\t\t SYS PATH\n"
    for item in sys.path:
        m1 += item + "\n"
    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 += "\n\t\t INSTALLING PIP PACKAGES ({})\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)

    # Install requirements ===================================
    with open('requirements.txt') as f:
        required = f.read().splitlines()
    for ipack in required:
        try:
            pkg, version = ipack.split(">=")[0:2]
            if pkg[0] == "#":
                continue
            install_with_pip(pkg, vers=version, log=logger)
        except ValueError:
            pkg = ipack
            if pkg[0] == "#":
                continue
            install_with_pip(pkg, log=logger)

    # Install Topology library from github =======================================
    install_topology_library(log=logger)
    # Check for topology installation
    try:
        import topology
    except ImportError:
        m = "ERROR. Topology library is not correctly installed\n"
        m +="ERROR. Please check ./topology/install.log"
        print(m) if logger is None else logger.info(m)
        exit()

    ## Install foyer ===================================
    #nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    #m1 = "\n\t\t INSTALLING FOYER (https://github.com/mosdef-hub/foyer) ({})\n".format(nowm)
    #print(m1) if logger is None else logger.info(m1)
    #install_foyer(log=logger)

    ## Install openmm ===================================
    #nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    #m1 = "\n\t\t INSTALLING OPENMM (https://openmm.org/) ({})\n".format(nowm)
    #print(m1) if logger is None else logger.info(m1)
    #install_openmm(log=logger)

    ## Install gmso ===================================
    #nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    #m1 = "\n\t\t INSTALLING GMSO (https://github.com/mosdef-hub/gmso) ({})\n".format(nowm)
    #print(m1) if logger is None else logger.info(m1)
    #install_gmso(log=logger)

    ## Install mbuild ===================================
    #nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    #m1 = "\n\t\t INSTALLING MBUILD (https://github.com/mosdef-hub/mbuild) ({})\n".format(nowm)
    #print(m1) if logger is None else logger.info(m1)
    #install_mbuild(log=logger)

    # Install intermol ===================================
    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t INSTALLING InterMol (https://github.com/shirtsgroup/InterMol) ({})\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)
    install_intermol(log=logger)

    # Setup Replicate ===========================================
    from Cython.Build import cythonize
    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    extensions = setup_external_extensions()
    m1 = "\n\t\t RUNNING SETUP FROM SETUPTOOLS {}\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)
    setup(ext_modules=cythonize(extensions,
          compiler_directives={'language_level': sys.version_info[0]}),
          include_dirs=[numpy.get_include()])

    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t Installation Done!!!! at {}\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)

    # Remove warnings for simtk
    remove_warnings()
