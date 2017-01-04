#!/usr/bin/env python

"""
Provide a command line interface to solvemwts procedure. 

User can specify inputs 
such as scenario, phase 1 dat file, output path, solver, mipgap, debug mode,
phase 1 model, and phase 2 model. It appears that there were plans to have
the ability to read a YAML formatted input config file - not implemented yet.
"""

import sys
import argparse

import pymwts.solvemwts as solve

def process_command_line(argv):
    """
    Return a Namespace representing the argument list.

    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    
    # If argv is empty, get the argument list from sys.argv.
    if argv is None:
        argv = sys.argv[1:]
        
    # Initialize parser object
    parser = argparse.ArgumentParser(
        description='Solve a multi-week tour scheduling problem.',
        epilog='May the force be with you.')
        
    # Add arguments
    parser.add_argument('--version', action='version', 
                        version='%(prog)s 1.0')
                        
    parser.add_argument('scenario',
                        help='Short string to be used in output filenames')
                        
    parser.add_argument('phase1dat',
        help='DAT file for phase 1')
        
    parser.add_argument('-p', '--path',
                        default='./',
        help='Relative path to output file directory. Terminate with /')
    
    parser.add_argument('-s', '--solver', choices=['cbc','glpk','gurobi'], 
                        default='cbc',
                        help='cbc or glpk for now')    
                        
    parser.add_argument('-t', '--timelimit', type=int, default=1e+6,
                        help='seconds')
                        
    parser.add_argument('-g', '--mipGap', type=float, default=1e-6,
                        help='Can prevent really long run times.')
                        
    parser.add_argument('-w', '--windebug', action='store_true',
                        help='Write out start window debug info.')
    
    parser.add_argument('-p1', '--phase1model', default='mwts_phase1.py',
                        help='Model for phase 1 problem') 
                        
    parser.add_argument('-p2', '--phase2model', default='mwts_phase2.py',
                        help='Model for phase 2 problem') 
                        
    parser.add_argument('-y', '--yaml',  
                        default=None,
                        help='YAML input config filename. NOT IMPLEMENTED.')    

    args = parser.parse_args()
    return args    

def main(argv=None):
    args = process_command_line(argv)
    print args
    # application code here, like:
    # run(settings, args)
    solve.solvemwts(args.scenario,args.phase1dat,args.path,
                        args.solver,args.timelimit,args.mipGap,
                        args.windebug)
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
