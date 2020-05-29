"""
Read input files for mwts problems and create a GMPL data file.
"""

# Author: misken
# License: TBD

import re
import json

def create_mwt(filenameInput,stubOutput,output_path):
    pattLineType = re.compile('^(TTS|PP4|n_weeks|n_tours)')    # The regex to match the file section headers

    #pattTourShift = re.compile('^tourshift\[([0-9]+),([0-9]+),([0-9]+),([0-9]+),([0-9]+),([0-9]+),([0-9]+)')
    outFiles = {}
    whichFile = None


#===============================================================================
# PP4
#  1  2  1 11  4 16 11  4
#  1  2  1 11  6 16 11  6
#  1  2  2 11  3 16 11  3
#  1  2  2 11  4 16 11  4
#  1  2  2 11  5 16 11  5
#  1  2  2 11  6 16 11  6
#  2  3  1 11  1 16 11  1
#  2  3  1 11  3 24 11  3
#  2  3  1 11  5 24 11  5
#  2  3  1 11  6 24 11  6
#  2  3  2 11  1 24 11  1
#  2  3  2 11  2 24 11  2
#  2  3  2 11  6 24 11  6
#===============================================================================


    # outFiles['TTS'] = open(output_path+stubOutput+'.tts','w')
    outFiles['MWT'] = open(output_path+stubOutput+'.mwt','w')
    # pp2key = []



    #outFiles['PP2'] = open(stubOutput+'.pp2','w')
    #outFiles['PP3'] = open(stubOutput+'.pp3','w')

    inFile = open(filenameInput)           # Open the file
    pp4mode = False
    toursInitialized = False
    num_weeks = 0
    num_tours = 0

    while 1:
        line = inFile.readline()           # Read a single line
        if not line: break                 # Bail if EOF

        matchobj = pattLineType.match(line)    # See if this line is start of new section
        if matchobj:
            if matchobj.group(1) == 'PP4':
                pp4mode = True
            else:
                pp4mode = False

            # If this is the num_weeks section, create the keys and files
            # for each week
            if matchobj.group(1) == 'n_weeks':
                num_weeks = int(line.split()[1])
                # for w in range(1,num_weeks+1):
                #     pp2key.append("pp2_%d" % (w))
                #     key = pp2key[w-1]
                #     outFiles[key] = open("%s%s_w%d.pp2" % (output_path, stubOutput, w),'w')

            if matchobj.group(1) == 'n_tours':
                num_tours = int(line.split()[1])

            # Initialize the tour data structure if not done yet
            if num_weeks>0 and num_tours>0 and not toursInitialized:
                tours = []
                mwtours = []
                tourtypes = []
                for w in range(1,num_weeks+1):
                    weeklytours = []
                    for t in range(1,num_tours+1):
                        tourspec = []

                        for _ in range(1,8):
                            tourspec.append(0)
                        tourspec.append(1)
                        for _ in range(1,8):
                            tourspec.append(0)
                        weeklytours.append(tourspec)
                    tours.append(weeklytours)

                for t in range(1,num_tours+1):
                    mwtourspec = []
                    for w in range(1,num_weeks+1):
                        for _ in range(1,8):
                            mwtourspec.append('x')
                    mwtours.append(mwtourspec)

                toursInitialized = True
                
        else:
            if pp4mode:
                tour = line.split()
                if len(tour) > 5:
                    tournum = int(tour[0])
                    ttype = int(tour[1])
                    week = int(tour[2])
                    startprd = int(tour[3])
                    dow = int(tour[4])
                    shiftlen = int(tour[5])


                    tours[week-1][tournum-1][dow-1] = startprd
                    tours[week-1][tournum-1][dow-1+8] = shiftlen

                    shift = '(' + str(startprd) + ':' + str(shiftlen) + ')'
                    mwtours[tournum-1][(week-1)*7+dow-1] = shift


                #whichFile.write(line)         # No match, just write the line to active output file
            else:
                # outFiles['TTS'].write(line)   # Write TTS line
                tourtypes.append(int(line))


    # for w in range(1,num_weeks+1):
    #     key = pp2key[w-1]
    #     for t in range(1,num_tours+1):
    #         for i in tours[w-1][t-1]:
    #              outFiles[key].write("%5d" % i)
    #         outFiles[key].write('\n')



    for t in range(1,num_tours+1):
        print (mwtours[t-1])
        # outFiles['MWT'].write(mwtours[t])
        # outFiles['MWT'].write(str(tourtypes[t-1]))
        json.dump(mwtours[t-1],outFiles['MWT'])
        outFiles['MWT'].write('\n')

    # Close all the files
    # for k in pp2key:
    #     outFiles[k].close

    outFiles['MWT'].close
    inFile.close


def main():
    #extract_outfiles(sys.argv[1],sys.argv[2])
    pass

if __name__ == '__main__':
    #extract_outfiles(sys.argv[1],sys.argv[2])
    pass

