# getResiduals.py
# A.WARD 2016 (original AWK script by M. Coombes)

"""
  Retrieve and save the time history of the residuals in a friendly format 
  for later plotting.
    
Usage
    python getResiduals.py e3mpi.xxxx.log
"""

import sys
import os

def parse(parseFilePath):
    outputFilePath = parseFilePath.split('.')[0] + "_residuals.dat"
    if os.path.isfile(outputFilePath):
        os.remove(outputFilePath)
    with open(outputFilePath, 'a') as out:
        with open(parseFilePath, 'r') as f:
            lines = f.readlines()
            massResiduals = []
            massTimes = []
            energyResiduals = []
            energyTimes = []
            for line in lines:
                if 'mass global' in line:
                    words = line.split()
                    massResiduals.append(float(words[4]))
                    massTimes.append(float(words[8]))
                if 'energy global' in line:
                    words = line.split()
                    energyResiduals.append(float(words[4]))
                    energyTimes.append(float(words[8]))
        out.write(str(massResiduals))
        out.write('\n')
        out.write(str(massTimes))
        out.write('\n')
        out.write(str(energyResiduals))
        out.write('\n')
        out.write(str(energyTimes))
    f.close()
    
if __name__ == '__main__':
    argv = sys.argv
    parse(argv[1])

