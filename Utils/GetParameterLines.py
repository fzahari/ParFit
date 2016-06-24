#!/usr/bin/env python
# change the way ParamList is created and populated by using python list methods.

def GetParameterLines(VariedCoord, NoOfVariedCoords):
    ParamList = []
    if VariedCoord == "diha":
        for m in range(NoOfVariedCoords):
            ParamList.append(raw_input("Enter parameter line number for the dihedral to be fit. "))
            for i in range(1, 4):
                ParamOrConst = raw_input("Enter 'p' if Line " + ParamList[m] + " V"
                    + str(i) + " should be varied during ParFit run. ")
                if ParamOrConst != "p":
                    ParamOrConst = "c"
                ParamList.append(ParamOrConst)
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
            q_index = int(ParamList.index(q))
            r_index = int(ParamList.index(r))
            s = int(s)
            ParamList[q_index + s] = ParamList[r_index + s] = "p" + str(n + 1)
    else:
        ParamList.append(raw_input("Enter parameter line number of the "
            + "parameters to be fit. ") + " p p")
#        print ParamList  #print for debugging
    return ParamList

print GetParameterLines("angl", 2)
