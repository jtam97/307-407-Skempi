import getopt
import math
import sys
import argparse

"""
GreedySearch main function
- Takes 2 PDB files, optional distance variable, optional output file NameError
"""

# Parse arguments from command line


def parseArg():

    # parse output to take two inputs -i
    parser = argparse.ArgumentParser(
        description='Identify all points between two proteins that are within a certain distance of each other.')
    parser.add_argument('-i', nargs=2, metavar='InputPDB',
                        help='Input PDB file to be compared.')
    parser.add_argument('-d', nargs='?', metavar='distance',
                        type=float, help='Resolution for distance checking.')
    parser.add_argument('-m', nargs='?', metavar='mode',
                        help='Search mode. Contact check = [], ionic bond check = i, hydrogen bond check = h, cation pi check = c')
    parser.add_argument('-o', nargs='?', metavar='OutputPDB',
                        help='Output PDB file name')

    # parse list of points from inputs
    args = parser.parse_args()
    args = vars(args)
    i1, i2 = args['i'][0], args['i'][1]
    d = args['d']
    m = args['m']
    o = args['o']
    return i1, i2, d, m, o

def main():
    parseArg()

if __name__ == "__main__":
    main()
