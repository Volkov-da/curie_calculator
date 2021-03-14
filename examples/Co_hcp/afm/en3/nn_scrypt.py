import siman
from siman.calc_manage import smart_structure_read
st = smart_structure_read('POSCAR')
st.printme()
st.magmom = [5, 5, -5, -5]
st.name = 'Fe'
st = st.replic([3, 3, 3])
out = st.nn(1, 30)
# out = st.nn(2, 30)
# out = st.nn(3, 30)
# out = st.nn(4, 30)

unique = []
# print(out['dist'])
print(out['el'])
nn = 0
neighbors = []
# neighbors{1} = 0

unique_dists = []
number_of_neig = 0
nJ = 0
nJ_list = []
dnn = out['dist'][1]
for el, d in zip(out['el'][1:], out['dist'][1:]):
    # print(el, d)
    # unique_dists.append(d)

    if abs(d - dnn) < 0.01:
        number_of_neig += 1
        if el == 'Fe':
            nJ += 1
        elif el == 'Co':
            nJ -= 1
    else:
        neighbors.append(number_of_neig)
        nJ_list.append(nJ)
        dnn = d
        number_of_neig = 1
        if el == 'Fe':
            nJ = 1
        elif el == 'Co':
            nJ = - 1
    # neighbors{nn} += 1
print('Number of neighbours:')
print(neighbors)
print('Coefficients:')
print(nJ_list)


# ss = [ ]
# ar = np.array(out['el'])

# print(ar[])
