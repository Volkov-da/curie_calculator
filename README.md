# Magnetic critical temperature Calculator

Code is currently under development. Now it working fine for few materials (e.g EuO).


```
1. git clone --recursive https://github.com/Volkov-da/curie_calculator.git

2. pip3 install siman

3. chmod +x install.sh (if needed)

4. ./install.sh
```

Export path to enumlib executables (you might need to specify absolute path):
Add this to your `.bashrc` or `.zshrc` etc. file.

```
5. export PATH="home/username/curie_calculator/enumlib/src:$PATH"
```


```
cd examples/EuS
python ../../new_src/curie_calculator.py
python ../../new_src/swaper.py
```