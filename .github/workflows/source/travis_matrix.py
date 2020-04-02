#!/usr/bin/env python

# Concatenation of tests in each of the 6 elements in Travis matrix
f = open('./travis_matrix_elements.txt') ; matrix_elements = f.readlines() ; f.close()
# All tests read by prepare_travis_tests.py
f = open('./travis_all_tests.txt') ; all_tests = f.readlines() ; f.close()

# Now let's make sure these two are equal

# Remove these elements from both lists, as they are are not test names
elements_to_remove = ['[main]\n', '[AMReX]\n', '[source]\n', '[extra-PICSAR]\n']
for element in elements_to_remove:
    for x in range(matrix_elements.count(element)):
        matrix_elements.remove(element)
    for x in range(all_tests.count(element)):
        all_tests.remove(element)

# Sort lists, and make sure they are equal
matrix_elements.sort()
all_tests.sort()
print("Tests in matrix, but not in initial list (typically if a test is done twice):")
print(list(set(matrix_elements) - set(all_tests)))
print("Tests in initial list but not in the matrix:")
print(list(set(all_tests) - set(matrix_elements)))

assert( matrix_elements == all_tests )
