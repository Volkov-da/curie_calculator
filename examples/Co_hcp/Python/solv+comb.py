import itertools
import numpy as np


comb = 0  # counter for the total number of combinations

fm = [2, 6, 6, 6] #-14.070633
afm1 = [2, 6, -6, 6] #-13.611630
afm2 = [4, 0, -6, 0] #-27.651858
afm3 = [4, 0, 2, 0] #-27.228849
afm4 = [6, 0, -6, 0] #-41.666455
afm5 = [6, 6, -6, 6] #-41.666457
afm6 = [6, 2, -2, 4] #-41.047883
afm7 = [6, 4, -2, 0] #-40.916701
afm8 = [6, 0, -2, 2] #-41.008527

values = [-14.070633, -13.611630, -27.651858, -27.228849, -41.666455, -41.666457, -41.047883, -40.916701, -41.008527]
data = [fm, afm1, afm2, afm3, afm4, afm4, afm5, afm6, afm7, afm8]
number_singula = 0
no = []
index = []
Eg = []
J1 = []
J2 = []
J3 = []

print('Enter the number of configurations:')
configurations_number = int(input()) #number of FM configurations
print('Enter the number of equations in one system:')
equations_number = int(input()) #number of equations in one linear system

for i in itertools.combinations(range(1, configurations_number + 1, 1), equations_number):
    n = list(i)
    no.append(n)
    comb += 1

counter = 0
print('Combinations:')
print(no)
print('Total number of combinations:')
print(comb)
print('#################################')
print('#################################')
print('#################################')
while counter < comb:
    s = no[counter]
    eq2 = s[0]
    eq3 = s[1]
    eq4 = s[2]
    A = [data[0], data[eq2], data[eq3], data[eq4]]
    b = [values[0], values[eq2], values[eq3], values[eq4]]
    print('Calculation #', counter + 1)
    print('Inputs:')
    print(A[0], b[0], 'fm')
    print(A[1], b[1], 'afm', eq2)
    print(A[2], b[2], 'afm', eq3)
    print(A[3], b[3], 'afm', eq4)
    if int(np.linalg.det(A)) == 0:
        print('_________________________________')
        print('Outs:')
        print('Singular matrix')
        number_singula += 1
        counter += 1
        print('   #     #   ')
        print('    #   #    ')
        print('     # #     ')
        print('    #   #    ')
        print('   #     #   ')
        print('_________________________________')
        print('#################################')
        print('_________________________________')
    else:
        index.append(no[counter])
        x = np.linalg.solve(A, b)
        print('_________________________________')
        print('Outs:')
        print('Eg =', "%.4f" % x[0], 'eV')
        Eg.append("%.4f" % x[0])
        print('J1 =', "%.4f" % x[1], 'eV')
        J1.append("%.4f" % x[1])
        print('J2 =', "%.4f" % x[2], 'eV')
        J2.append("%.4f" % x[2])
        print('J3 =', "%.4f" % x[3], 'eV')
        J3.append("%.4f" % x[3])
        print('_________________________________')
        print('#################################')
        print('_________________________________')
        counter += 1
print('Number of singular matrix:', number_singula)
print('Indexes:', index)
print('Eg=', Eg)
print('J1=', J1)
print('J2=', J2)
print('J3=', J3)
