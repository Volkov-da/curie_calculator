# Magnetic critical temperature Calculator
---
## Installation
Code is currently under development. Now it working fine for few materials (e.g EuO).


```
git clone --recursive https://github.com/Volkov-da/curie_calculator.git

pip3 install siman

chmod +x install.sh (if needed)

./install.sh
```

Export path to enumlib executables (you might need to specify absolute path):
Add this to your `.bashrc` or `.zshrc` etc. file.

```
5. export PATH="home/username/curie_calculator/enumlib/src:$PATH"
```
---

## How to run examples

```
cd examples/EuS
python ../../new_src/file_builder.py
python ../../new_src/swaper.py
python ../../new_src/solver.py
```
