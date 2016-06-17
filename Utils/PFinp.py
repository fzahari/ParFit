#!/usr/bin/env python
#
# This a ParFit input file generating program.
#
#  To Do:
#  1. Error out if dihedral is given incorrectly (such as two numbers are the
#     same or not all numbers are given)
#  2. For values that are not known, the program will leave a place holder that a
#     user can then change using a text editor.
#  3. Add option to change two or three parameters per line.
#  4. Add new parameter fitting options, bond length, bond angle.
#

# --- Functions ---
def quantumdata() :
    print "Choose the QM data format:"
    print "\t(a) Compact: one file that includes fixed torsion angle geometries."
    print "\t(b) Series: GAMESS log files, one for each fixed torsion angle geometry."
    qmdatachoice = raw_input( "Enter: a or b. Default is a.\n" )
    if ( qmdatachoice == "a" ) :
        qmdata = 'comp'
        filenameroot = raw_input( '''\nEnter the root file name. It should match the name of the compact file containing energies and geometries, minus the word 'scan'.\n''' )
    elif ( qmdatachoice == "b" ) :
        qmdata = 'full'
        filenameroot = raw_input('''\nEnter the root file name. It should match the root file name of your GAMESS log files minus '***.log' where *** is an angle.\n''' )
    else :
        qmdata = 'comp'
        filenameroot = raw_input('''\nEnter the root file name. It should match the name of
the compact file containing energies and geometries, minus the word 'scan'.\n''' )

# --- Determine ParFit input file name ---

pyout = raw_input( "Enter the name of ParFit input file to create, if blank, the file name will be PFinput.\n" )
if ( pyout == "" ) :
    pyout = "PFinput"
    print "[PFinp]: Input filename:", pyout, "\n"
else :
    pyout == pyout
    print "[PFinp]: Input filename:", pyout, "\n"

# --- Open the file for writing ---
f = open(pyout,'w')

# Select the parameter type, bond length, bond angle, torsion.
property_type = raw_input( '''Choose from the properties below:
(a) bond length
(b) bond angle
(c) torsion (default)
    \n
Enter: a, b, or c.\n''' )

if ( property_type == "a" ) :
    property_type = 'bond'
    parameterize = "bond length"
elif ( property_type == "b" ) :
    property_type = 'angl'
    parameterize = "bond angle"
elif ( property_type == "c" ) :
    property_type = 'diha'
    parameterize = "torsion angle"
else :
    print "The default, torsion angle, was chosen."
    property_type = 'diha'
    parameterize = "torsion angle"

# Multiple dihedral angle file fitting.
if ( property_type == "bond" ) :
############# insert making an input file for a bond parameter fitting
    print "Place holder for bond parameter fitting."
elif ( property_type == "angl" ):
############# insert making an input file for a bond parameter fitting
    print "Place holder for bond parameter fitting."
elif ( property_type == "diha" ) :
    no_torsions = int( raw_input( ''' How many torisions are to be fit?\n''' ) )
    print >> f, "mult, ", no_torsions
    for n in range( 0, no_torsions ) :
# --- Create GAMESS input files or use existing energy/geometry data. ---
        qmdatachoice = raw_input( '''Choose from the scenarios below:
        (a) I have compact file that includes all of the geometry and energy information
            for the torsion angles described above.
        (b) I have a GAMESS output file for each torsion angles in the range described above.
            \n
        Enter: a or b. Default is a.\n''' )
        if ( qmdatachoice == "a" ) :
            qmdata = 'comp'
            filenameroot = raw_input("\nEnter the root file name. It should match the name of the\ncompact file containing energies and geometries, minus the word 'scan'.\n" )
        elif ( qmdatachoice == "b" ) :
            qmdata = 'full'
            filenameroot = raw_input("\nEnter the root file name. It should match the root file\nname of your GAMESS log files minus '***.log' where *** is an angle.\n" )
        else :
                qmdata = 'comp'
                filenameroot = raw_input("\nEnter the root file name. It should match the name of the\ncompact file containing energies and geometries, minus the word 'scan'.\n" )

# --- Description of Molecule and Rotation used for the Fit ---

        torsion = raw_input( "What are the indices of the four atoms creating the dihedral angle to be fit?\n" )
        TorInit = raw_input( "What is the initial torsion angle? \nFor the default of 0 degrees, press enter.\n" )
        if ( TorInit == "" ) :
            TorInit = "0"
        else :
            TorInit = TorInit
        TorFin  = raw_input( "What is the final torsion angle?\n" )
        TorStep = raw_input( "What is the angle step size?\nFor the default of 5 degrees, press enter.\n" )
        if ( TorStep == "" ) :
            TorStep = "5"
        else :
            TorStep = TorStep
        print >> f, "{0}, {1}, {2}, {3} {4} {5}".format( qmdata , filenameroot , torsion , TorInit , TorFin , TorStep )

# if not running a bond length/angle run, use this path.
else :
    print "[PFinp] Error: You have not properly chosen a property to parameterize."

# --- Get engine path ---

engine_path = raw_input( "\nWhat is the full engine.exe path?\n" )

# --- Determine the type of MM file that is to be modified ---

mmtypchoice = raw_input( "\nChoose the MM type (mm3 or mmff94) parameters to be fit\n(a) MM3 - default\n(b) MMFF94\nChoose a or b.\n" )
if ( mmtypchoice == 'a' ) :
    carbontyp = 50
    mmtyp = 'mm3'
elif ( mmtypchoice == 'b' ) :
    carbontyp = 37
    mmtyp = 'mmff94'
else :
    mmtyp = 'mm3'
    print "\nWarning: Check the MM type you entered, the only options are a and b. Default will be chosen.\n"
print >> f, mmtyp

# --- Choose the algorithm used to fit parameters. ---

print "\nPlease choose from the following options for algorithm to be used."
print "(a) genetic algorithm"
print "(b) Nedler-Mead algorithm"
print "(c) hybrid: genetic followed by Nedler-Mead algorithm - default"
alg = raw_input("For default, just press enter.\n")
if ( alg == 'a' ) :
    alg = 'ga'
elif ( alg == 'b' ) :
    alg = 'fmin'
elif ( alg == 'c' ) :
    alg = 'hybr'
else :
   alg = 'hybr'
   print "Default was chosen."
print >> f, alg

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
        print "\nWarning: check the parameter in line {0} that should be fit.\n".format( line_no )
    print >> f, "{0} {1}".format( line_no , param )

# --- Printing csv file option ---

printcsv = raw_input( "\nEnter \"n\" if you do NOT want ParFit to print a csv format file\ncontaining the angles, QM energy, and the optmized MM energies.\n" )
if ( printcsv == 'n') :
    csv = "csv_off"
else :
    csv = "csv_on"
print >> f, csv

print "\nYour ParFit input file name {0} has been generated.\n".format( pyout )
exit()
