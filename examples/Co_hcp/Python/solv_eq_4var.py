import numpy as np

fm = [2, 6, 6, 6]
afm1 = [2, 6, -6, 6]
afm2 = [4, 0, -6, 0]
afm3 = [4, 0, 2, 0]
afm4 = [6, 6, -6, 6]
afm5 = [6, 2, -2, 4]
afm6 = [6, 4, -2, 0]
afm7 = [6, 0, -2, 2]
data = [fm, afm1, afm2, afm3, afm4, afm4, afm5, afm6, afm7]
values = [-14.0706, -13.61163, -27.651858, -27.228849, -41.666455, -41.666457, -41.047883, -40.916701, -41.008527]


eq2 = 1
eq3 = 5
eq4 = 6


A = [data[0], data[eq2], data[eq3], data[eq4]]
b = [values[0], values[eq2], values[eq3], values[eq4]]
print('_________________________________')
print('Inputs:')
print(A[0], b[0], 'fm')
print(A[1], b[1], 'afm', eq2)
print(A[2], b[2], 'afm', eq3)
print(A[3], b[3], 'afm', eq4)
x = np.linalg.solve(A, b)
print('_________________________________')
print('Outs:')
print('Eg =', "%.4f" % x[0], 'eV')
print('J1 =', "%.4f" % x[1], 'eV')
print('J2 =', "%.4f" % x[2], 'eV')
print('J3 =', "%.4f" % x[3], 'eV')
print('_________________________________')