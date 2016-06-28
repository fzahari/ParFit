#!/usr/bin/env python
#
# This a ParFit input file generating program.
#
#  To Do:
#  1. Error out if dihedral is given incorrectly (such as two numbers are the
#     same or not all numbers are given)
#  2. For values that are not known, the program will leave a place holder that a
#     user can then change using a text editor.
#

####################################################################################

# --- Functions ---
def getPESCoords():
    print("The following prompts define the variables used to generate" +
          "the QM PESs.")
    atomIndices = raw_input("Enter atom indices separated by a space.\n")
    iCoord = raw_input("Enter the initial length or angle.\n")
    fCoord = raw_input("Enter the final length or angle.\n")
    sCoord = raw_input("Enter the length or angle step size.\n")
    return "{0}, {1} {2} {3}".format(atomIndices, iCoord, fCoord, sCoord)

def getQMFileFormat():
    print("<< QM data file format >>\n" +
          "(a) Compact: QM data is contained in one file that includes fixed\n" +
          "    bond lengths, bond angles, or torsion angle geometries.\n" +
          "(b) Series: GAMESS log files, one for each fixed bond length,\n" +
          "    bond angle, or torsion angle geometry.\n" +
          "Enter: a or b. [a]")
    fileFormat = raw_input()
    if fileFormat != "b" :
        fileFormat = "comp"
    else:
        fileFormat = "full"
    return fileFormat

def getQMFilePropertyLines(variedCoord):
    filePropertyLines = []
    if (variedCoord == "bond" or variedCoord == "angl"):
        QMFilenameRoot = raw_input("Enter the root filename. ")
        noOfVariedCoords = 1
        QMFormatChoice = getQMFileFormat()
        if QMFormatChoice == "comp":
            filePropertyLines.extend((QMFormatChoice, QMFilenameRoot, variedCoord))
        elif QMFormatChoice == "full":
            filePropertyLines.extend((QMFormatChoice, QMFilenameRoot, getPESCoords(), variedCoord))
    elif (variedCoord == "diha"):
        noOfVariedCoords = int(raw_input("Enter the number of PESs to be fit.\n"))
        for n in range(0, noOfVariedCoords):
            QMFormatChoice = getQMFileFormat()
            QMFilenameRoot = raw_input("Enter the root filename. ")
            if QMFormatChoice == "comp":
                filePropertyLines.extend((QMFormatChoice, QMFilenameRoot, variedCoord))
            elif QMFormatChoice == "full":
                filePropertyLines.extend((QMFormatChoice, QMFilenameRoot, getPESCoords(), variedCoord))
    return filePropertyLines, noOfVariedCoords

def getParameterLines(variedCoord):
    paramList = []
    if variedCoord == "diha":
        noOfParameterLines = int(raw_input("Enter number of parameter lines. "))
        for m in range(noOfParameterLines):
            paramList.append(raw_input("Enter parameter line number for the "
                + "dihedral to be fit. "))
            for i in range(1, 4):
                paramOrConst = raw_input("Enter 'p' if Line " + paramList[m] + " V"
                    + str(i) + " should be varied during ParFit run. ")
                if paramOrConst != "p":
                    paramOrConst = "c"
                paramList.append(paramOrConst)
        i = raw_input("How many pairs of parameters are coupled? [0] ")
        if i == "":
            i = 0
        else:
            i = int(i)
            print("Please identify the coupled parameters by giving the line numbers \n"
                + "and parameter (1, 2 or 3 for V1, V2 and V3) "
                + "in the following format: \n"
                + "\t[line number] [line number] [parameter number] ")
        for n in range(i):
            q, r, s = str.split(raw_input("Enter the line and parameter numbers. "))
            qIndex = int(paramList.index(q))
            rIndex = int(paramList.index(r))
            s = int(s)
            paramList[qIndex + s] = paramList[rIndex + s] = "p" + str(n + 1)
    else:
        paramList.append(raw_input("Enter parameter line number of the "
            + "parameters to be fit. ") + " p p")
#        print paramList  #print for debugging
    return paramList
###################################################################################
###################################################################################

# --- Determine ParFit input file name ---
print("Enter the name of ParFit input file to create, if blank, the file " +
      "name will be PFinput.\n")
outFileName = raw_input()
if (outFileName == ""):
    outFileName = "PFinput"
    print "[PFinp]: Input filename:", outFileName, "\n"
else:
    outFileName == outFileName
    print "[PFinp]: Input filename:", outFileName, "\n"

# --- Open the file for writing ---
f = open(outFileName,'w')

# --- Select the parameter type, bond length, bond angle, torsion. ---
property_type = raw_input('''Choose from the properties below:
(a) bond length
(b) bond angle
(c) torsion (default)
    \n
Enter: a, b, or c.\n''')

if (property_type == "a"):
    variedCoord = 'bond'
    parameterize = "bond length"
elif (property_type == "b"):
    variedCoord = 'angl'
    parameterize = "bond angle"
elif (property_type == "c"):
    variedCoord = 'diha'
    parameterize = "torsion angle"
else:
    print "The default, torsion angle, was chosen."
    variedCoord = 'diha'
    parameterize = "torsion angle"

# --- Multiple dihedral angle file fitting. ---
lines, noOfVariedCoords = getQMFilePropertyLines(variedCoord)


# --- Determine which parameters will be changed by ParFit ---
paramList = getParameterLines(variedCoord)

# --- Get engine path ---
enginePath = raw_input("\nWhat is the full engine.exe path?\n")

# --- Determine the type of MM file that is to be modified ---
MMTypeChoice = raw_input("\nChoose the MM type (mm3 or mmff94) parameters to be fit\n(a) MM3 - default\n(b) MMFF94\nChoose a or b.\n")
if (MMTypeChoice == 'a' or MMTypeChoice == ""):
    carbontyp = 50
    MMType = 'mm3'
elif (MMTypeChoice == 'b'):
    carbontyp = 37
    MMType = 'mmff94'
else:
    MMType = 'mm3'
    print("Warning: Check the MM type you entered, the only options are a and b.\n" +
          "Default will be chosen.")

# --- Choose the algorithm used to fit parameters. ---
print("<< Fitting alorithms >>\n" +
      "(a) genetic algorithm\n" +
      "(b) Nedler-Mead algorithm\n" +
      "(c) hybrid: genetic followed by Nedler-Mead algorithm\n")
alg = raw_input("Select a, b, or c. [c] ")
if (alg == 'a'):
    alg = 'ga'
elif (alg == 'b'):
    alg = 'fmin'
elif (alg == 'c'):
    alg = 'hybr'
else:
   alg = 'hybr'
   print "Default was chosen."

# --- Printing csv file option ---
printcsv = raw_input("\nEnter \"n\" if you do NOT want ParFit to print a csv format file\ncontaining the angles, QM energy, and the optmized MM energies.\n")
if (printcsv == 'n'):
    csv = "csv_off"
else:
    csv = "csv_on"

# --- Print out input file ---
if (variedCoord == "diha"):
    print >> f, "mult, " + str(noOfVariedCoords)
for n in range(0, len(lines)):
    if lines[n] == "full":
        print >> f, ", ".join(lines[n: n+4])
for s in range(0, len(lines)):
    if lines[s] == "comp":
        print >> f, ", ".join(lines[s: s+3])
print >> f, enginePath
print >> f, MMType
print >> f, alg
if variedCoord == "diha":
    for i in range(0, len(paramList), 4):
        n = i + 4
        print >> f, " ".join(paramList[i:n])
else:
    print >> f, paramList
#else:                 #print for debugging
#    print paramList  #print for debugging
print >> f, csv

print "\nYour ParFit input file name {0} has been generated.\n".format(outFileName)

exit()
