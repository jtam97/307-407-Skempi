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

# Parsing PDB here gives several useful outputs. First, all lines that do not contain ATOMS are completely ignored. Then, it breaks up both PDB files into :
# - a list of coordinates of every atom
# - a list of lines as they are in the PDB file. We can write back to PDB format easier by just referencing the line index.
# - a list of lines that have been split by space so we can use specific parts of each line for different parsing methods in the future.


def parsePDB(PDB1, PDB2):
    f = open(PDB1, "r")
    templines1 = f.readlines()
    f.close()

    splitLines_PDB1 = []
    unsplitLines_PDB1 = []
    # X,Y,Z	coordinates	found in lines[6],lines[7],lines[8]
    for x in range(len(templines1)):
        if templines1[x][0] == 'A':
            unsplitLines_PDB1.append(templines1[x])
            splitLines_PDB1.append(templines1[x].split())

    f = open(PDB2)
    templines2 = f.readlines()
    f.close()

    splitLines_PDB2 = []
    unsplitLines_PDB2 = []
    # X,Y,Z	coordinates	found in lines[6],lines[7],lines[8]
    for x in range(len(templines2)):
        if templines2[x][0] == 'A':
            unsplitLines_PDB2.append(templines2[x])
            splitLines_PDB2.append(templines2[x].split())

    # coordinates are added to a list and then appended to make a list of points
    points1 = []
    for x in splitLines_PDB1:
        if x[0] == "ATOM":
            temp = []
            temp.append(float(x[6]))
            temp.append(float(x[7]))
            temp.append(float(x[8]))
            points1.append(temp)

    points2 = []
    for x in splitLines_PDB2:
        if x[0] == "ATOM":
            temp = []
            temp.append(float(x[6]))
            temp.append(float(x[7]))
            temp.append(float(x[8]))
            points2.append(temp)

    return points1, points2, unsplitLines_PDB1, splitLines_PDB1, unsplitLines_PDB2, splitLines_PDB2


# Checks a connected graph from all points in 1st PDB file to 2nd PDB file.
# If within dist, then append index of point to output
# This point index can be referenced using splitLines (for further parsing using other parameters) or unsplitLines (to write back to PDB format easier)
def compareDist(points1, points2, dist):
    output1 = []
    for i in range(len(points1)):
        temp = [i]
        temp.append([])
        for j in range(len(points2)):
            d1 = (points1[i][0]-points2[j][0])**2
            d2 = (points1[i][1]-points2[j][1])**2
            d3 = (points1[i][2]-points2[j][2])**2
            d = math.sqrt(d1+d2+d3)
            if d < dist:
                temp[1].append(str(j))
        if temp[1]:
            output1.append(temp)
    return output1
    # return output1, output2

# Output format is the same as above but adds an extra parameter of checking if amino acids are oppositely charged AND checks for a distance of 3


def compareDistIonic(output, points1, points2, dist, splitLines_PDB1, splitLines_PDB2):
    output1 = []
    for i in range(len(output)):
        pts1_idx = int(output[i][0])
        temp_AminoAcid = splitLines_PDB1[pts1_idx][3]
        charge1 = 0
        temp = [pts1_idx]
        temp.append([])
        if temp_AminoAcid == 'HIS' or temp_AminoAcid == 'ARG' or temp_AminoAcid == 'LYS':
            if ('NZ' in splitLines_PDB1[pts1_idx][2] or 'NE' in splitLines_PDB1[pts1_idx][2] or 'ND' in splitLines_PDB1[pts1_idx][2] or 'NH' in splitLines_PDB1[pts1_idx][2]):
                charge1 = 1
        elif temp_AminoAcid == 'ASP' or temp_AminoAcid == 'GLU':
            if ('OD' in splitLines_PDB1[pts1_idx][2] or 'OE' in splitLines_PDB1[pts1_idx][2]):
                charge1 = -1
        if charge1 != 0:
            for j in range(len(output[i][1])):
                charge2 = 0
                pts2_idx = int(output[i][1][j])
                temp_AminoAcid = splitLines_PDB2[pts2_idx][3]
                if temp_AminoAcid == 'HIS' or temp_AminoAcid == 'ARG' or temp_AminoAcid == 'LYS':
                    if ('NZ' in splitLines_PDB2[pts2_idx][2] or 'NE' in splitLines_PDB2[pts2_idx][2] or 'ND' in splitLines_PDB2[pts2_idx][2] or 'NH' in splitLines_PDB2[pts2_idx][2]):
                        charge2 = 1
                elif temp_AminoAcid == 'ASP' or temp_AminoAcid == 'GLU':
                    if ('OD' in splitLines_PDB2[pts2_idx][2] or 'OE' in splitLines_PDB2[pts2_idx][2]):
                        charge2 = -1
                total_charge = charge1 + charge2
                if total_charge == 0:
                    # d1 = (points1[pts1_idx][0]-points2[pts2_idx][0])**2
                    # d2 = (points1[pts1_idx][1]-points2[pts2_idx][1])**2
                    # d3 = (points1[pts1_idx][2]-points2[pts2_idx][2])**2
                    # d	=	math.sqrt(d1+d2+d3)
                    # if d < dist:
                    temp[1].append(str(pts2_idx))
        if temp[1]:
            output1.append(temp)
    return output1


def compareDistCatPi(output, points1, points2, dist, splitLines_PDB1, splitLines_PDB2):
    output1 = []
    for i in range(len(output)):
        pts1_idx = int(output[i][0])
        temp_AminoAcid = splitLines_PDB1[pts1_idx][3]
        charge1 = 0
        temp = [pts1_idx]
        temp.append([])
        if temp_AminoAcid == 'ARG' or temp_AminoAcid == 'LYS':
            charge1 = 1
        elif temp_AminoAcid == 'PHE' or temp_AminoAcid == 'TYR' or temp_AminoAcid == 'TRP':
            charge1 = -1
        if charge1 != 0:
            for j in range(len(output[i][1])):
                charge2 = 0
                pts2_idx = int(output[i][1][j])
                temp_AminoAcid = splitLines_PDB2[pts2_idx][3]
                if temp_AminoAcid == 'ARG' or temp_AminoAcid == 'LYS':
                    charge2 = 1
                elif temp_AminoAcid == 'PHE' or temp_AminoAcid == 'TYR' or temp_AminoAcid == 'TRP':
                    charge2 = -1
                total_charge = charge1 + charge2
                if total_charge == 0:
                    # d1 = (points1[pts1_idx][0]-points2[pts2_idx][0])**2
                    # d2 = (points1[pts1_idx][1]-points2[pts2_idx][1])**2
                    # d3 = (points1[pts1_idx][2]-points2[pts2_idx][2])**2
                    # d	=	math.sqrt(d1+d2+d3)
                    # if d < dist:
                    temp[1].append(str(pts2_idx))
        if temp[1]:
            output1.append(temp)
    return output1


# Removes duplicates and also parses out any atoms with # = -1 or 0
# Since changing the code to use only line indexes, duplicates will not occur so this function is unused
def removeDupe(output):
    output = list(dict.fromkeys(output))
    for i in output:
        if i == -1:
            output.remove(-1)
        if i == 0:
            output.remove(0)
    return output


def numToAA(PDB_split_lines):
    list = {}
    for i in PDB_split_lines:
        if i[5] not in list.keys():
            list[i[5]] = i[3]
    return list


# Output is structured so that we have a list of lists. Each list contains:
# - 0th index = line index of atom in PDB1 that is within dist of another atom in PDB2
# - 1st index = [] list of indexes of atoms in PDB2 in contact with above atom
# Writes back to PDB by getting the 0th index of each list item in output.
# These are all atoms that have a contact, but does not tell you what atoms in PDB2 are in contact.
def writePDB(PDBlines, output, outputFileName):
    f = open("Contact_" + outputFileName + ".pdb", "w+")
    for i in output:
        f.write(PDBlines[i[0]])
    f.close
    print("Contact PDB written to Contact_" + outputFileName + ".pdb")


# ##This will give a list in CSV format of all atoms in PDB1 with a contact, and a list following that of atoms in PDB2 in contact with the PDB1 atom.
# def writeContactListTemp(unsplitLines_PDB1, splitLines_PDB1, unsplitLines_PDB2, splitLines_PDB2, output, outputFileName, list1, list2):
# 	AAindex1 = -1
# 	write_list = []
# 	count = -1
# 	dupe_list = []
# 	for i in output:
# 		#check if first amino acid is the same as previous in list. If it is not, then we just go through the list and add amino acids to list as usual.
# 		if splitLines_PDB1[i[0]][5] != AAindex1:
# 			dupe_list = []
# 			AAindex1 = splitLines_PDB1[i[0]][5]
# 			AAindex2 = -1
# 			for j in i[1]:
# 				AAindex2 = splitLines_PDB2[int(j)][5]
# 				if (AAindex2 not in dupe_list):
# 					dupe_list.append(AAindex2)
# 			temp = [AAindex1]
# 			temp.append(dupe_list)
# 			write_list.append(dupe_list)
# 			count = count + 1
# 		#If first amino acid is the same as previous, that means we want to merge the two lists without duplicates.
# 		#We take the most recent list from write_list and will merge into this.
# 		else:
# 			AAindex2 = -1
# 			for j in i[1]:
# 				AAindex2 = splitLines_PDB2[int(j)][5]
# 				# if splitLines_PDB2[int(j)][5] != AAindex2:
# 				if (AAindex2 not in dupe_list):
# 					dupe_list.append(AAindex2)
# 					# temp.append(AAindex2)
# 			temp = [AAindex1]
# 			temp.append(dupe_list)
# 			write_list[-1] = temp
# 	f=open("ContactList_" + outputFileName + ".pdb", "w+")
# 	#reference dictionary from AA index to AA name for writing to file
# 	for i in write_list:
# 		temp = list1[i[0]] + " "+ str(i[0])
# 		for j in i[1]:
# 			temp = temp + ", " + list2[j] + " " + str(j)
# 		temp = temp + "\n"
# 		f.write(temp)
# 	f.close
# 	print("Contact list written to ContactList_" + outputFileName + ".pdb")


# This will give a list in CSV format of all atoms in PDB1 with a contact, and a list following that of atoms in PDB2 in contact with the PDB1 atom.
def writeContactList(unsplitLines_PDB1, splitLines_PDB1, unsplitLines_PDB2, splitLines_PDB2, output, outputFileName, list1, list2, dist):
    AAindex1 = -1
    write_list = []
    count = -1
    dupe_list = []
    for i in output:
        # check if first amino acid is the same as previous in list. If it is not, then we just go through the list and add amino acids to list as usual.
        if splitLines_PDB1[i[0]][5] != AAindex1:
            AAindex1 = splitLines_PDB1[i[0]][5]
            # dupe_list = []
            dupe_list = i[1]
            temp = [i[0]]
            temp.append(dupe_list)
            write_list.append(temp)
        # If first amino acid is the same as previous, that means we want to merge the two lists without duplicates.
        # We take the most recent list from write_list and will merge into this.
        else:
            for j in i[1]:
                if j not in dupe_list:
                    dupe_list.append(j)
            write_list[-1][1] = dupe_list

    # print(write_list)
    f = open("ContactList_" + outputFileName + "_" + str(dist) + ".csv", "w+")

    for i in write_list:
        dupe_list = []
        temp = "PDB1 AA: " + \
            str(splitLines_PDB1[i[0]][3]) + \
            str(splitLines_PDB1[i[0]][5]) + " | "
        i[1].sort()
        for j in i[1]:
            a = str(splitLines_PDB2[int(
                j)][4]) + " " + str(splitLines_PDB2[int(j)][3]) + str(splitLines_PDB2[int(j)][5])
            dupe_list.append(a)
            # temp = temp + str(splitLines_PDB2[int(j)][4]) + " " + str(splitLines_PDB2[int(j)][3]) + str(splitLines_PDB2[int(j)][5]) + ", "
            # dupe_list.append(tempAA)
        # print(dupe_list)
        final = []
        for k in dupe_list:
            if k not in final:
                final.append(k)
                temp = temp + str(k) + ", "
        temp = temp + "\n"
        # print(temp)
        f.write(temp)
    f.close
    print("Contact list written to ContactList_" +
          outputFileName + "_" + str(dist) + ".pdb")


def main():
    input1, input2, dist, mode, outputFileName = parseArg()
    if dist == None:
        dist = 5

    if mode == None:
        mode = "contact"

    if outputFileName == None:
        outputFileName = input1.split(".")[0] + "_" + input2.split(".")[0]

    # list of points from first input PDB or second input PDB
    # lines1, lines2 are original parsed lines from PDB files used for writing back to output files
    points1, points2, unsplitLines_PDB1, splitLines_PDB1, unsplitLines_PDB2, splitLines_PDB2 = parsePDB(
        input1, input2)

    # parse pdb file for all atoms
    # atoms1,atoms2	=	parsePDB(PDB.inputs[0],PDB.inputs[1])
    #output1,output2	=	compareDist(points1, points2, dist)

    # creates a dictionary matching each AA index number to the AA name
    list1 = numToAA(splitLines_PDB1)
    list2 = numToAA(splitLines_PDB2)

    output = compareDist(points1, points2, dist)
    # print(output)
    output_ionic = compareDistIonic(
        output, points1, points2, dist, splitLines_PDB1, splitLines_PDB2)
    # print(output_ionic)
    # output_catpi = compareDistCatPi(
    #     output, points1, points2, dist, splitLines_PDB1, splitLines_PDB2)
    # print(output_catpi)

    # write to PDB file
    writePDB(unsplitLines_PDB1, output, outputFileName)
    # write list of contacts
    # writeContactList(unsplitLines_PDB1, splitLines_PDB1, unsplitLines_PDB2,
    #                  splitLines_PDB2, output, outputFileName, list1, list2, dist)

    # write list of ionic bonds in contacts
    outputFileName_csv = outputFileName + "_ionic"
    writeContactList(unsplitLines_PDB1, splitLines_PDB1, unsplitLines_PDB2,
                     splitLines_PDB2, output_ionic, outputFileName_csv, list1, list2, dist)

    # outputFileName_csv = outputFileName + "_CationPi"
    # writeContactList(unsplitLines_PDB1, splitLines_PDB1, unsplitLines_PDB2,
    #                  splitLines_PDB2, output_catpi, outputFileName_csv, list1, list2, dist)


if __name__ == "__main__":
    main()
