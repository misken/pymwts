#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      IPS User
#
# Created:     28/06/2011
# Copyright:   (c) IPS User 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

def main():
    pass

if __name__ == '__main__':
    main()


def filterpatterns(x):
    """
    Creates a sequence of binary values to be used for list filtering. This
    function will contain the various rules used to filter out weekend days
    worked patterns that we don't want to allow.
    """
    keep = True

    if ((x[0]==1 and x[27]==1 and x[13]==1 and x[14]==1) or (x[6]==1 and x[7]==1 and x[20]==1 and x[21]==1)):
        if not (sum(x) == 4):
            keep = False # Making sure no more than 4 weekend days
        else:                                 # worked over the scheduling horizon
            keep = True
    else:
        keep = False

    return keep


# Assumptions for this version
# (1) 4 week horizon
# (2) Only dealing with (Sun,Sat) weekends - weekend type 1
# (3) All tour types have same set of weekend patterns

outfilename = 'wkends_everyother.txt'
numweeks = 4
numttypes = 6
patterns = []
goodpatterns = []

outfile = open(outfilename,'w')

for sat4 in range(0,2):
    for sun4 in range(0,2):
        for sat3 in range(0,2):
            for sun3 in range(0,2):
                for sat2 in range(0,2):
                    for sun2 in range(0,2):
                        for sat1 in range(0,2):
                            for sun1 in range(0,2):
                                #print sun1, sat1, sun2, sat2, sun3, sat3, sun4, sat4
                                patterns.append([sun1, 0,0,0,0,0, sat1, sun2, 0,0,0,0,0, sat2, sun3, 0,0,0,0,0, sat3, sun4, 0,0,0,0,0, sat4])

# Diagnostic printing
for p in patterns:
    if sum(p)==4 and p[0]==1 and p[27]==1:
        print p,sum(p)
goodpatterns.extend(filter(filterpatterns,patterns))

print goodpatterns
print len(goodpatterns)

# Print GMPL dat matrix

# param A{i in 1..max_weekend_patterns,j in DAYS,w in WEEKS,t in TTYPES,e in WEEKENDS: i <= num_weekend_patterns[t]} default 0;

# [*,*,*,1,1]
# 1	1	1	0

for t in range(numttypes):
    outfile.write('\n[*,*,*,%i,1]\n' % (t+1))
    patnum = 0
    for p in goodpatterns:
        patnum = patnum + 1
        for w in range(1,numweeks+1):
            for j in range(1,8):
                position = (w-1)*7 + j
                val = p[position-1]
                line = '%i %i %i %i\n' % (patnum,j,w,val)
                outfile.write(line)

outfile.close()
