#!/usr/bin/env python

#    error out if dihedral is given incorrectly (such as two numbers are the same or not all numbers are given)

# --- Determine ParFit input file name ---

pyout = raw_input( "Enter the name of ParFit input file to create, if blank, the file name will be PFinput.\n" )           # take file name form stdin
if ( pyout == "" ) :
    pyout = 'PFinput'
else :
    pyout == pyout

# --- Open the file for writing ---

f = open(pyout,'w')   # ................................................  open a file for writing ('w')

# --- Create GAMESS input files or use existing energy/geometry data. ---

runtyp = raw_input( '''Choose from the scenarios below:
(a) I have compact file that includes all of the geometry and energy information
    for the torsion angles described above.
(b) I have a GAMESS output file for each torsion angles in the range described above.
(c) I need GAMESS input files to run a series of constrained optimizations with the
    torsion angles described above.
    \n
Enter: a, b, or c. Default is a.\n''' )  # Prompt for defining run type options.

if ( runtyp == "a" ) :
    runtyp = 'comp'
elif ( runtyp == "b" ) :
    runtyp = 'full'
elif ( runtyp == "c" ) :
    runtyp = 'ginp'
else :
    runtyp = 'comp'
#dict = {}   #........................................................ Run type definitions.
#dict['a'] = 'comp'
#dict['b'] = 'full'
#dict['c'] = 'ginp'

# --- Description of Molecule and Rotation used for the Fit ---

torsion   = raw_input( "What are the indices of the four atoms creating the dihedral angle to be fit?\n" )
TorInit = raw_input( "What is the initial torsion angle?\n" )
TorFin  = raw_input( "What is the final torsion angle?\n" )
TorStep = raw_input( "What is the angle step size?\n" )

# --- Create short form input file ---

if ( runtyp == 'ginp' ) :
    filenameroot = raw_input( "Enter the filename root of the series of GAMESS input files to be generated.\nFormat: String with no spaces.\n" )
    onlyline = '{0}, {1}, {2}, {3} {4} {5}'.format( runtyp , filenameroot , torsion , TorInit , TorFin , TorStep ) # formats only line in short form input file.
    print >> f,onlyline
    print "\nYour ParFit input file name {0} has been generated.\n".format( pyout )
    exit()

# --- Create long form input file ---

elif ( runtyp == 'comp' or 'full' ) :

# --- Get engine path ---

    engine_path = raw_input( "\nWhat is the full engine.exe path?\n" )

# --- Identify and format the double bonds found in the molecule ---
#
#    n= int( raw_input( "\nNumber of double bonds in the molecule:\n" ) ) # take input and convert to an integer
#    print "\nEnter the pair of atoms making up the double bonds.\nFor multiple double bonds, enter the atom numbers pairwise, #pressing return after each pair of atoms.\n"
#    formatDblBndStr = ""
#    for i in range( n ) :     # ........................................... loop prints comma separated pairs of integers designating #atom pairs forming double bonds.
#           DblBndA,DblBndB = raw_input().split()     # .................... split input (pair of integers) into two strings
#           DblBndA,DblBndB = int( DblBndA ),int( DblBndB )    # ........... convert the string values to integers
#           formatDblBndStrElem="{0} {1}, "              # ................. format integers into a comma separated values
#           formatDblBndStr += formatDblBndStrElem.format( DblBndA , DblBndB )
#    Double_bonds = formatDblBndStr[:-2]

# --- Determine the type of MM file that is to be modified ---

    mmtyp = raw_input( "\nChoose the MM type (mm3 or mmff94) parameters to be fit\n(a) MM3\n(b) MMFF94\nChoose a or b.\n" ) # a = mm3 and b = mmff94 MM type
    if ( mmtyp == 'a' ) :
        carbontyp = 50
        mmtyp = 'mm3'
    elif ( mmtyp == 'b' ) :
        carbontyp = 37
        mmtyp = 'mmff94'
    else :
        print "\nWarning: Check the MM type you entered, the only options are a for mm3 and b for mmff94\n"

# --- Choose the algorithm used to fit parameters. ---

    print "\nPlease choose the fitting algorithm."
    alg = raw_input("Enter ga for genetic algorithm, or fmin for simplex algorithm.\n")

# --- Determine which parameters will be changed by ParFit ---

    print "\nNow you will be prompted to enter the line numbers that contain the parameters to be fit.\n"
    m = int( raw_input( "\nHow many parameters in add_{0}.prm are to be fit?\n".format( mmtyp ) ) )
    print "\nYou have {0} parameters to fit. When prompted, please enter each line number followed by the parameter designation.\n".format( m )
    prm_lines = ""
    for i in range( m ) :
        line_no = raw_input( "\nLine number:\n" )
        var_param = raw_input( "\nWhich parameter in line {0} is to be fit?\n\t(a) first\n\t(b) second\n\t(c) third\n".format( line_no ) )
        if ( var_param == 'a' ) :
            param = "p c c"
        elif ( var_param == 'b' ) :
            param = "c p c"
        elif ( var_param == 'c' ) :
            param = "c c p"
        else :
            param = "p c c"
            print "\nWARNING: check the parameter in line {0} that should be fit. Input file will be printed with default values: p c c.\n\
open input file to modify.".format( line_no )
        formatedline = "{0} {1}\n".format( line_no , param )
        prm_lines += formatedline

# --- Obtain file name root ---

    if ( runtyp == 'comp' ) :
        filenameroot = raw_input("\nEnter the root file name. It should match the name of the\ncompact file containing energies and geometries, minus the word 'scan'.\n" )
    elif ( runtyp == 'full') :
        filenameroot = raw_input("\nEnter the root file name. It should match the root file\nname of your GAMESS log files minus '***.log' where *** is an angle.\n" )

# --- Print csv file option ---

    printcsv = raw_input( "\nEnter \"n\" if you do NOT want ParFit to print a csv format file\ncontaining the angles, QM energy, and the optimized MM energies.\n" )
    if ( printcsv == 'n') :
        csv = "csv_off"
    else :
        csv = "csv_on"
#        print "\nSorry, I didn't understand your input. The default, 'yes' will be set and a csv file will be printed.\n"

# --- Format and print the ParFit input file ---

    inputfile = '''{0}, {1}, {2}, {3} {4} {5}\n{6}\n{7}\n{8}\n{9}{10}'''\
            .format( runtyp , filenameroot , torsion , TorInit , TorFin , TorStep , \
            engine_path , \
            mmtyp , \
            alg , \
#           carbonlist , \
            prm_lines , \
            csv )

    print >> f,inputfile
    print "\nYour ParFit input file name {0} has been generated.\n".format( pyout )
    exit()

else :
    print "\nError: the only options are a, b, or c. Please start over.\n"
    exit()

