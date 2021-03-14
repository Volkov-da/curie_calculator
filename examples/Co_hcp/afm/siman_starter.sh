#!/bin/bash
for dir in */; do
	cd "$dir"
	python3.8 nn_scrypt.py > log_magnetic
	cd ..
done
