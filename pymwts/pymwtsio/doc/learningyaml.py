#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      IPS User
#
# Created:     05/10/2011
# Copyright:   (c) IPS User 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import yaml
from numpy import *

def learningyaml():
    fin = open("infiles/yaml_new_01.mix","r")
    mixspec = yaml.load(fin)
    print mixspec
    print mixspec['mixes']
    print mixspec['mixes'][0]
    print mixspec['mixes'][0]['min_days_week']

    lenset = set([])
    for m in mixspec['mixes']:
        for s in m['shiftlengths']:
            lenset.add(s['numbins'])
    print lenset
    lengths = list(lenset)
    lengths.sort()
    print lengths
    param = list_to_param("lengths",lengths)
    print param

def list_to_param(pname,plist):
    """
    Convert a list to a GMPL representation of a parameter.
    """
    print 'param ' + pname + ':='

    # Convert parameter as list to an ndarray
    parray = array(plist)

    # Denumerate the array to get at the index tuple and array value
    paramrows = ndenumerate(parray)
    param = ''
    for pos, val in paramrows:
        poslist = [str(p + 1) for p in pos]
        datarow = ' '.join(poslist) + ' ' + str(val) + '\n'
        param += datarow

    param += ";\n"
    return param



def main():
    learningyaml()

if __name__ == '__main__':
    main()
