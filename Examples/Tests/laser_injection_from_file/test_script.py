#!/usr/bin/python3

#ADD COMMENT
#ADD COMMENT
#ADD COMMENT

import sys
import glob
import os

def launch_analysis(file):
    print(file)
    os.system("./test_injection.py --generate_txye_files")
    os.system("./" + file + " inputs.2d_test_txye")
    os.system("./test_injection.py --check ")

def main() :
    files = glob.glob("main2d*")
    if len(files) == 1 :
        launch_analysis(files[0])
        return True
    else :
        return False


if __name__ == "__main__":
    main()
