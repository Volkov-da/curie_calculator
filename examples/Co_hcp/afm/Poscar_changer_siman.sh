#!/bin/bash
sed -i.bak $'s/  2  /  Co Fe\\\n  1 1/g' en*/POSCAR
sed -i.bak $'s/  4  /  Co Fe\\\n  2 2/g' en*/POSCAR
sed -i.bak $'s/  6  /  Co Fe\\\n  3 3/g' en*/POSCAR
sed -i.bak $'s/  8  /  Co Fe\\\n  4 4/g' en*/POSCAR
sed -i.bak -e 's/^D/Dir/g' en*/POSCAR
