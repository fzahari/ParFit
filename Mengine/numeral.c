# Copyright (C) 2004-2016 Kevin E. Gilbert, dba Serena Software
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Kevin E. Gilbert
# Serena Software
# Box 3076
# Bloomington, IN 47402-3076
#
# gilbert@serenasoft.com

#include "pcwin.h"
#include "pcmod.h"

void numeral(int number, char *astring, int size)
{
    int multi,length,minsize,million,hthou,tthou;
    int thou, hun,tens,ones;
    int right, negative;
    char bstring[4];

    if (size == 0)
    {
        right = FALSE;
        size = 1;
    } else
        right = TRUE;

    minsize = size;
    length = size;

    if (number > 0)
       negative = FALSE;
    else
    {
        negative = TRUE;
        number = -number;
    }
    million = number/1000000;
    multi = 1000000*million;
    hthou = (number-multi)/100000;
    multi = multi + 100000*hthou;
    tthou = (number-multi)/10000;
    multi = multi + 10000*tthou;
    thou = (number-multi)/1000;
    multi = multi + 1000*thou;
    hun = (number-multi)/100;
    multi = multi + 100*hun;
    tens = (number-multi)/10;
    multi = multi + 10*tens;
    ones = number - multi;

    if (million != 0)
       size = 7;
    else if (hthou != 0)
       size = 6;
    else if (tthou != 0)
       size = 5;
    else if (thou != 0)
       size = 4;
    else if (hun != 0)
       size = 3;
    else if (tens != 0)
       size = 2;
    else
       size = 1;

    if ( length < size)
       size = length;
    if ( minsize > size)
       size = minsize;

    sprintf(bstring,"%d",number);
    size = strlen(bstring);

    if (size == 3)
    {
      astring[0] = bstring[0];
      astring[1] = bstring[1];
      astring[2] = bstring[2];
      astring[3] = '\0';
    } else if (size == 2)
    {
      bstring[2] = ' ';
      astring[0] = bstring[2];
      astring[1] = bstring[0];
      astring[2] = bstring[1];
      astring[3] = '\0';
    } else if (size == 1)
    {
      bstring[2] = ' ';
      astring[0] = bstring[2];
      astring[1] = bstring[2];
      astring[2] = bstring[0];
      astring[3] = '\0';
    }
        
}       

