import itertools
comb = 0 #counter for the total number of combinations
print('Number of configurations:')
configurations_number = int(input()) #number of FM configurations
print('Number of equations in one system:')
equation_number = int(input()) #number of equations in one linear system
no = []
for i in itertools.combinations(range(1, configurations_number + 1, 1), equation_number):
    n = list(i)
    no.append(n)
    comb += 1
print(no)

counter = 0
print(comb)
while counter < comb:
    print(no[counter])
    s = no[counter]
    print(s)
    eq1 = s[1]
    print(eq1)
    counter += 1


