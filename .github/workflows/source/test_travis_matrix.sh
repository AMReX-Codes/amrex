#! /usr/bin/env sh

set -e

cp .github/workflows/source/travis_matrix.py Regression/
cd Regression/

# Put the name of all travis tests into a text file
python prepare_file_travis.py
grep "\[" travis-tests.ini > travis_all_tests.txt

# Concatenate the names of all elements in Travis matrix into another test file
WARPX_CI_REGULAR_CARTESIAN=TRUE python prepare_file_travis.py
grep "\[" travis-tests.ini >  travis_matrix_elements.txt
WARPX_CI_PSATD=TRUE             python prepare_file_travis.py
grep "\[" travis-tests.ini >> travis_matrix_elements.txt
WARPX_CI_PYTHON_MAIN=TRUE       python prepare_file_travis.py
grep "\[" travis-tests.ini >> travis_matrix_elements.txt
WARPX_CI_SINGLE_PRECISION=TRUE  python prepare_file_travis.py
grep "\[" travis-tests.ini >> travis_matrix_elements.txt
WARPX_CI_RZ_OR_NOMPI=TRUE      python prepare_file_travis.py
grep "\[" travis-tests.ini >> travis_matrix_elements.txt
WARPX_CI_QED=TRUE               python prepare_file_travis.py
grep "\[" travis-tests.ini >> travis_matrix_elements.txt

# Check that the resulting lists are equal
{
    python travis_matrix.py &&
        rm travis_all_tests.txt travis_matrix_elements.txt travis_matrix.py &&
        echo "test passed" &&
        exit 0
} || {
    rm travis_all_tests.txt travis_matrix_elements.txt travis_matrix.py &&
        echo "tests failed" &&
        exit 1
}
