# AMARAsurf.py
# A.WARD & D.WARD, 2016
"""
This file iterates through the block looking for the boundary 
condition keyword specified on line 45. It then returns a file, 
AMARAblock_surfaceList.dat giving a list of block numbers and 
face name (NORTH, EAST etc) for all surfaces using that 
boundary condition.
Note it currently DOES NOT work correctly for SuperBlock3D objects.


Usage
# Run AMARAsurf.py, giving it your file specifying the blocks.
$ python AMARAsurf.py AMARAblock.py

# Paste content of AMARAblock_surfaceList.dat into the e3post command:
$ e3post.py --job=AMARAgeom --vtk-xml --surface-list="<output>"
"""

import sys
import os

def parse(parseFilePath):
    #Lines to search after occurence of Block3D or SuperBlock3D
    searchDist = 15
    outputFilePath = parseFilePath.split('.')[0] + "_surfaceList.dat"
    if os.path.isfile(outputFilePath):
        os.remove(outputFilePath)
    with open(outputFilePath, 'a') as out:
        with open(parseFilePath, 'r') as f:
            lines = f.readlines()
            b3dCount = 0
            counts = []
            words = []
            linenums = []
            for i, line in enumerate(lines):
                if 'Block3D' in line:
                    if '#' in line: continue
                    startTBCSearch = i+1
                    foundTBC = False
                    for j in range(startTBCSearch, startTBCSearch + searchDist):
                        if j > len(lines) - 1 or 'Block3D' in lines[j]:
                            break
                        if '#' in lines[j]: continue
                        if 'KEYWORD' in lines[j]:
                            foundTBC = True
                            try:
                                w = lines[j].split('[')[1].split(']')[0]
                            except:
                                w = ''
                                print "Bad Bracket Split / KeyWord"
                            if w in ['BOTTOM', 'WEST', 'SOUTH', 'TOP', 'NORTH', 'EAST']: 
                                counts.append(b3dCount)
                                linenums.append(j)
                                words.append(w)
                            else:
                                print "bad line:", j
                    b3dCount += 1
                    if not foundTBC:
                        continue
                        print "No 'FixedTBC' in", str(searchDist), "lines following 'Block3D' found at: line:", i
            print "Lines searched / Processed: ", len(lines)
        writeStr = ''
        for i,count in enumerate(counts):
            #Output file
            writeStr = writeStr + str(count) + ',' + str(words[i]) + ';'
            #Linenum debug
            #writeStr = writeStr + str(count) + ',' + str(words[i]) + ',' + str(linenums[i] + 1) + ';\n'
        #if you dont want the ';' on the end of the last occurence uncomment the last line
        writeStr = writeStr[:-1]
        out.write(writeStr)

if __name__ == '__main__':
    argv = sys.argv
    parse(argv[1])
