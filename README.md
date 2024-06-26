# Magnetic critical temperature Calculator :magnet:
![Algo Details](/images/algo.png)
## :warning: Warning
For all the calculations you need VASP ab-initio program. 
Project is no more supported.
---
## :computer: Installation
Code is currently under development. Now it working fine for few materials (e.g EuO).


```
git clone --recursive https://github.com/Volkov-da/curie_calculator.git
```

1. `./install.sh`

2. `python -m venv .venv`

3. `source .venv/bin/activate`

4. `pip install -r requirements.txt`

if needed:

```
chmod +x install.sh
```
```
chmod +x src/curie_calculator.py
```

Export path to _enumlib_ and _curie_calculator_ executables (you might need to specify absolute path):
Add this to your `.bashrc` or `.zshrc` file. Also, it might be useful to create an alias for running a python script.

```
export PATH="home/username/curie_calculator/enumlib/src:$PATH"
export PATH="/home/username/curie_calculator/src:$PATH"
alias curie_calculator='python ~/curie_calculator/src/stat_file_builder.py'
```
---

## :compass: How to run examples

To run any example, you simply need a `POSCAR` file in the folder. Also, it is important to have automated access to pseudopotential (`POTCAR`) files used in VASP (i.e., `POT_GGA_PAW_PBE` or `POT_LDA_PAW` folders). For these purposes please check how to use `.pmgrc.yaml` file.

```
cd examples/EuO_test
curie_calculator
```
