import itertools
import numpy as np

comb = 0  # counter for the total number of combinations
fm = [1, 8, 6] #-8.245
afm1 = [2, 0, -4] #-16.030
afm2 = [2, -16, 12] #-15.565
afm3 = [4, 16, 8] #-32.615
afm4 = [4, 0, -8] #-32.170
afm5 = [4, 0, 8] #-32.228
afm6 = [4, 0, -24] #-32.200
values = [-8.245, -16.030, -15.565, -32.615, -32.170, -32.228, -32.200]
data = [fm, afm1, afm2, afm3, afm4, afm4, afm5, afm6]
no = []
Eg = []
J1 = []
J2 = []

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
print('O', no)
print('Total number of combinations:')
print(comb)
print('#################################')
print('#################################')
print('#################################')
while counter < comb:
    s = no[counter]
    eq2 = s[0]
    eq3 = s[1]
    A = [data[0], data[eq2], data[eq3]]
    print(A)
    b = [values[0], values[eq2], values[eq3]]
    print(b)
    print('Calculation #', counter + 1)
    print('Inputs:')
    print(A[0], b[0], 'fm')
    print(A[1], b[1], 'afm', eq2)
    print(A[2], b[2], 'afm', eq3)
    if int(np.linalg.det(A)) == 0:
        print('_________________________________')
        print('Singular matrix')
        counter += 1
        print('_________________________________')
        print('#################################')
        print('_________________________________')
    else:
        x = np.linalg.solve(A, b)
        print('_________________________________')
        print('Outs:')
        print('Eg =', "%.4f" % x[0], 'eV')
        Eg.append("%.4f" % x[0])
        print('J1 =', "%.4f" % x[1], 'eV')
        J1.append("%.4f" % x[1])
        print('J2 =', "%.4f" % x[2], 'eV')
        J2.append("%.4f" % x[2])
        print('_________________________________')
        print('#################################')
        print('_________________________________')
        counter += 1
print('Eg=', Eg)
print('J1=', J1)
print('J2=', J2)
