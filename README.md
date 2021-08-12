# How to install
 
## Installation of pymol with Conda:

- First install dependencies 
  - `sudo apt-get install git build-essential python3-dev libglew-dev libpng-dev libfreetype6-dev libxml2-dev libmsgpack-dev python3-pyqt5.qtopengl libglm-dev libnetcdf-dev`
- Install conda, follow the instructions on their website
- Create a new environment and switch to it
  - `conda crete -n myenv`
  - `conda activate myenv`
- Install pymol in the new environment
  - `conda install -c schrodinger pymol`
- Install python dependencies for the plugin
  - `conda install networkx numpy scipy seaborn pandas`
- Open pymol and go to Plugin > Plugin Manager > Settings
  - Add a new directory (/home/user/.pymol/startup)
- Copy the plugin folder in the specified directory
- Restart Pymol
- The plugin should be showing on the Plugin menu of Pymol
