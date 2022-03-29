import pandas as pd
import getopt
import math
import sys
import argparse

def read_PDB():
    combined = pd.read_csv("combined.csv", index_col=0)

def main():
    read_PDB()


if __name__ == "__main__":
    main()