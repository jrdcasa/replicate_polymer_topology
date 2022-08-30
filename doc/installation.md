# 3. Installation instructions

## Pre-installation

---

Previous to install the program  the following software have to be installed in the system:

* `git` (sudo apt-get install git)
* `python3-venv` (sudo apt update; sudo apt-get install python3-venv)
* `cmake` >=3.17
* `doxygen` (sudo apt-get install doxygen)
* `swig` (sudo apt-get install swig)
* `python3-dev` (sudo apt-get install python3-dev)

The version of cmake request is greater than 3.17. By default in Ubuntu 20.04 the version installed from the repositores is 3.16. Thus, the way to install cmake is the following:

```bash
sudo apt remove --purge --auto-remove cmake
sudo apt update
sudo apt install -y software-properties-common lsb-release 
sudo apt clean all

wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | sudo tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null
sudo apt-add-repository "deb https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main"

sudo apt update
sudo apt install kitware-archive-keyring
sudo rm /etc/apt/trusted.gpg.d/kitware.gpg

sudo apt update
sudo apt-get install cmake
cmake --version
```

## Create a virtualenv environment

---

First, we create an isolated Python environment to install the required packages (see dependencies below). Then, activate the virtual environment.

```bash
python3 -m venv <name_of_env>
source <name_of_env>/bin/activate
pip list
pip install --upgrade pip
```

This virtual environment must be activate in order to use the program.

## Prerequisites for topology library

---

This program makes use of the [``topology`` library](<https://github.com/jrdcasa/topology.git>). Thus, some packages or programs must be installed at hand before to start the ``replicate_polymer`` installation. These requisites are needed to properly install the ``topology`` library.

```bash
sudo apt get python3-dev
sudo apt get libgraphviz-dev
python -m pip install wheel
python -m pip install pygraphviz
```

The topology library can be installed using the following [instructions](https://github.com/jrdcasa/topology/blob/main/docs/02-installation.md):

```bash
git clone https://github.com/jrdcasa/topology.git
```
```

### Install from source

---

The setup.py script should install all dependencies automatically

```bash
git clone https://github.com/jrdcasa/replicate_polymer
cd replicate_polymer
python setup.py install
```

## Dependencies

---

To use **replicate_polymer**, the following libraries and software will need to be installed.

### Pip packages

---

Several packages from **pip** need to be installed. The **requirements.txt** contains the list of packages to be installed

* Cython>=0.29.24
* argcomplete>=1.12.3
* numpy>=1.21.2
* GitPython>=3.1.1
* ParmEd>=3.4.3
* lark-parser>=0.12.0
* ele>=0.2.0
* pydantic>=1.8.2
* unyt>=2.8.0
* lxml>=4.6.3
* boltons>=21.0.0
* networkx>=2.6.3
* requests>=2.26.0
* mdtraj>=1.9.6

### Package installed from the source

---

Packages **openmm**, **mbuild**, **foyer**, **gmso** and **intermol** need to be installed from the source code. In the best scenario these programs will be installed by running **setup.py**.

* **openmm**: <https://github.com/openmm/openmm.git>. No modifications from the original
* **intermol_cj**: <https://github.com/openmm/intermol_cj.git>
* **topology**: <https://github.com/jrdcasa/topology.git>. Library to implement topology from xsd, pdb and gro files

After installation you can check if all depencies are correctly installed:

 pip install
 replicate_polymer -h
