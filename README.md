# Magnetic critical temperature Calculator :magnet:
![Algo Details](/images/algo.png)

---
## :computer: Installation
Code is currently under development. Now it working fine for few materials (e.g EuO).


```
git clone --recursive https://github.com/Volkov-da/curie_calculator.git

pip3 install siman

chmod +x install.sh (if needed)
chmod +x src/curie_calculator.py

./install.sh
```

Export path to enumlib executables (you might need to specify absolute path):
Add this to your `.bashrc` or `.zshrc` file.

```
export PATH="home/username/curie_calculator/enumlib/src:$PATH"
export PATH="/Users/dmitry.volkov/curie_calculator/src:$PATH"
```
---

## :compass: How to run examples

```
cd examples/EuO_test
curie_calculator.py
```
