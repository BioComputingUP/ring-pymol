# RING-PyMOL
RING-PyMOL is a plugin for [PyMOL](https://pymol.org/) providing a set of analysis tools for structural ensembles and molecular dynamic (MD) simulations. RING-PyMOL combines residue interaction networks, as provided by the **RING software** ([RING source](https://biocomputingup.it/download), [RING webserver](https://ring.biocomputingup.it/)), with structural clustering to enhance the analysis and visualization of the conformational complexity. RING software has been improved extensively. It is ten times faster, can process mmCIF files and it identifies typed interactions also for nucleic acids.

Go to the [Wiki](https://github.com/BioComputingUP/ring-pymol/wiki) for the full documentation.


RING-PyMOL features
- Precise calculation of non-covalent interactions 
- Identifies and highlights correlating contacts and interaction patterns that can explain structural allostery, active sites and structural heterogeneity connected with molecular function
- Easy to use and extremely fast, processing and rendering hundreds of models and long trajectories in seconds
- Generates a number of interactive plots and output files for use with external tools

Reference 

RING 3.0: fast generation of probabilistic residue interaction networks from structural ensembles /
Clementel D, Del Conte A, Monzon AM, Camagni GF, Minervini G, Piovesan D and Tosatto SCE /
(2022) Nucleic Acids Research [https://doi.org/10.1093/nar/gkac365]

# Install
In order to work, Ring-PyMOL requires PyMOL and some Python packages. To simplify 
installation you can use the Singularity container provided below and skip completely
the manual install section.

## Manual install (HARD)
### Dependencies 
We provide three diffent solutions to insatll PyMOL and RING-PyMOL dependecies. The two using Conda gives the same result,
the one using Apt works only for Linux users and might give different results depending on the operating system.

**NOTE**
Sometimes PyMOL comes with a custom Python, please make sure you install the
dependencies for the correct PyMOL executable. 
In Linux you can type `which pymol` to see the path of the PyMOL executable. 
If installed with Conda, the command should return something like `/opt/miniconda3/envs/myenv/bin/pymol`

#### Conda from yml file (RECOMMENDED)

- Install [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) or [miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Download the environment.yml file from this repository
- Create the environment `conda env create -f environment.yml`
- Activate the environment `conda activate ring-pymol-plugin`

#### Conda commands
Same as before but with all commands issued explicitly

- Install [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) or [miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Create a new environment and activate it
    - `conda create -n myenv`
    - `conda activate myenv`
- Install PyMOL
    - `conda install -c conda-forge -c schrodinger pymol-bundle` (shrodinger version)
    - `conda install -c conda-forge pymol-open-source` (open-source version)
- Install Python dependencies
    - `conda install networkx numpy scipy seaborn pandas requests biopython`
    - `pip install qt-material` (This will install qt-material in the conda environment)

#### Apt commands (Linux only)
This will use system environment. 
The version of PyMOL and Python packages depend on the version and distribution of your operating system (OS).
RING-PyMOL might not work with an obsolete OS.

- `sudo apt install pymol python3-pip python3-tk`
- `pip install pmw networkx numpy~=1.20 scipy seaborn pandas qt-material biopython requests`

### RING-PyMOL plugin
Once you have installed PyMOL and all the RING-PyMOL dependencies you have to install the RING-PyMOL plugin. 

- Open PyMOL and go to Plugin > Plugin Manager > Install New Plugin > Install from Repository > Add..
    - Add https://ring.biocomputingup.it/plugin/
- Click on ring-plugin.zip in the right panel and then Install
- Set the installation directory
- The plugin should now appear on the Plugin menu of PyMOL

## Install container (EASY)

Another option for installing the plugin is to use the [singularity](https://docs.sylabs.io/guides/latest/admin-guide/) container definition file provided in this
repository.
To create the image file you can follow these steps:

- `sudo singularity build -F ring-pymol-plugin.sif singularity.def` (this will create the image file)
- `singularity shell --cleanenv --writable-tmpfs -B ~/.Xauthority ring-pymol-plugin.sif` (this will open a shell in the
  container).
  Note that the -B option is needed to allow the container to access the X server of the host machine for displaying the
  GUI.
- Start PyMOL with `pymol`
- Add a new directory where to find new plugins
    - Plugin > Plugin Manager > Settings > Add new directory...
    - Add `/opt`
- Restart PyMOL
- The plugin should now appear on the Plugin menu of PyMOL

# Usage

Go to the [Wiki](https://github.com/BioComputingUP/ring-pymol/wiki) for the full documentation.
