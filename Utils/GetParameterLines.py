#!/usr/bin/env python
# change the way ParamList is created and populated by using python list methods.


def GetParameterLines(VariedCoord, NoOfVariedCoords):
    if VariedCoord == "diha":
        NoOfParamLines = NoOfVariedCoords
        ParamListLen = NoOfParamLines*4
        ParamList = [None] * ParamListLen
        for m in range(0, ParamListLen, 4):
            ParamList[m]= raw_input("Enter parameter line number for the dihedral to be fit. ")
            for i in range(1, 4):
                c_or_p = raw_input("Enter 'p' if Line " + ParamList[m] + " V"
                    + str(i) + " should be varied during ParFit run. ")
                if c_or_p == "p":
                    c_or_p = "p"
                else:
                    c_or_p = "c"
                ParamList[m + i] = c_or_p
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
        m = raw_input("Enter parameter line number of the parameters to be fit. ")
        ParamList = m + " p p"
#        print ParamList  #print for debugging
    return ParamList

print GetParameterLines("diha", 2)
