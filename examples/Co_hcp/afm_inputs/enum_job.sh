#!/bin/bash
/Users/dmitry.volkov/enumlib/src/enum.x
rm VERSION.enum symops_enum_parent_lattice.out readcheck_enum.out debug_*
s=$(line=$(tail -1 struct_enum.out); echo ${line:0:11})
python3.8 /Users/dmitry.volkov/enumlib/aux_src/makeStr.py 1 $s
echo $s >> log
for ((i=1; i<=s; i++))
do
mkdir en"$i"
mv vasp.$i en$i
cp INCAR_en KPOINTS POTCAR jobscript nn_scrypt.py en$i
cd en$i
mv vasp.$i POSCAR
mv INCAR_en INCAR
cd ..
sed -i'' -e 's/job-name=Co-hcp/job-name=Co-hcp'$i'/g' en$i/jobscript
done
sed -i'' -e 's/^  1   1 /  2/g' */POSCAR
sed -i'' -e 's/^  2   2 /  4/g' */POSCAR
sed -i'' -e 's/^  3   3 /  6/g' */POSCAR
sed -i'' -e 's/^  4   4 /  8/g' */POSCAR

for ((i=1; i<=s; i++))
do
cd en$i
p=$(head -n 6 POSCAR | tail -1)
if [ ${p:0:3} = 4 ];  then
  sed -i'' -e 's/MAGMOM = 5 -5/MAGMOM = 5 5 -5 -5/g' INCAR
elif [ ${p:0:3} = 6 ]; then
  sed -i'' -e 's/MAGMOM = 5 -5/MAGMOM = 5 5 5 -5 -5 -5/g' INCAR
elif [ ${p:0:3} = 8 ]; then
  sed -i'' -e 's/MAGMOM = 5 -5/MAGMOM = 5 5 5 5 -5 -5 -5 -5/g' INCAR
fi
cd ..
done
#all copies of changed files, which was created(with "-e" index) will be removed:
rm */jobscript-e
rm */POSCAR-e
rm */INCAR-e
