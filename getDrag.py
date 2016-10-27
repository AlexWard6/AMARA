# getDrag.py
# A.WARD, D.KING, D.WARD 2016

"""
This is the start of a script to pull out lift and drag forces 
for a three dimensional VTK grid. Harder than it seems at first!
"""

import os
import numpy as np

# SET DIRECTORY WHERE YOUR GRID FILES ARE
directory = ''

os.chdir(directory)

'list of unknown length currently, summate at end to get total force'
xforce_block_list = []
x = np.array([1.0, 0.0, 0.0], dtype=float)

for filename in os.listdir(directory):
    
    if filename.endswith('.vtu'): 

        with open(filename, 'r') as f:
            lines = f.readlines()
            
            '''points on line 2 = line 1 in python'''
            points_position_start = lines[1].index('NumberOfPoints="') + len('NumberOfPoints="')
            points_position_end = lines[1].index('" NumberOfCells')
            blocks_positon_start = lines[1].index('NumberOfCells="') + len('NumberOfCells="')
            blocks_positon_end = lines[1].index('">')
            points = int( lines[1][points_position_start : points_position_end ] )
            blocks = int( lines[1][blocks_positon_start : blocks_positon_end ] )
            '''
            create empty matricies
            3d array block length with 4 points stacked vertically
            '''
            blockData = np.zeros( (blocks, 4, 3), dtype=float )
            pointData = np.zeros( (points, 3), dtype=float )
            connectionData = np.zeros( (blocks, 4), dtype=int )
            directionVectors = np.zeros( (blocks, 3), dtype=float )
            rhos = np.zeros( (blocks), dtype=float)
            pressures = np.zeros( (blocks), dtype=float)
            velos = np.zeros( (blocks, 3), dtype=float)
            xforces = np.zeros( (blocks), dtype=float)
        
            '''get the  data on the point xyz position for corners of blocjs'''
            i = 4        
            while i < (4 + points):
                data = lines[i].split()
                pointData[ (i - 4): ] = [ float(data[0]), float(data[1]), float(data[2]) ]
                i += 1
            
            '''add 4 to get to connectivity'''
            i += 4
            nLocal = i
            while i < ( nLocal + blocks ):
                data = lines[i].split()
                connectionData[ (i - nLocal): ] = [ int(data[0]), int(data[1]), int(data[2]), int(data[3]) ]
                i += 1
                
            '''yeah, go ahead laugh but it works so f u brah'''
            while i < 20000:
                if 'Name="rho"' in lines[i]:
                    nLocal = i
                    i += 1
                    break
                i += 1
            
            while i < ( nLocal + blocks ):
                rhos[ (i - nLocal ) ] = float(lines[i])
                i += 1
            
        
            while i < 20000:
                if 'Name="p"' in lines[i]:
                    nLocal = i
                    i += 1
                    break
                i += 1
                
            while i < ( nLocal + blocks ):
                pressures[ (i - nLocal ) ] = float(lines[i])
                i += 1
                
                
            
            while i < 20000:
                if 'Name="vel.vector' in lines[i]:
                    nLocal = i
                    i += 1
                    break
                i += 1
                
            while i < ( nLocal + blocks ):
                data = lines[i].split() 
                velos[ (i - nLocal ): ] = [ float(data[0]), float(data[1]), float(data[2]) ]
                i += 1
                

                
            '''create these stoopid blocks
            why can't you just fucking print them in the first place PJ??'''            
            for j in xrange(0, len(connectionData)):
                data = connectionData[j]
                blockData[j][0] = pointData[ data[0] ]
                blockData[j][1] = pointData[ data[1] ]
                blockData[j][2] = pointData[ data[2] ]
                blockData[j][3] = pointData[ data[3] ]
                
                
                ''' This is critical, I have no idea which way the vectors are!!!'''
                '''going to assume ADxAB would be the outwards direction'''
                AB = pointData[ data[1] ] - pointData[ data[0] ]
                AD = pointData[ data[3] ] - pointData[ data[0] ]                
                
                vector = np.cross(AD, AB)

                Area = np.linalg.norm(vector)

                '''
                use the direction vector (divided by the magnitude of original vector <- just incas u is n00b)
                also cross product happen to be area of parallelogram thing so we'll use that in calcs m8
                '''
                directionVectors[j] = vector/np.linalg.norm(vector)
                
                xforces[j] = Area*pressures[j]*-1*np.dot( directionVectors[j], x ) + rhos[j]*Area*( ( np.linalg.norm(velos[j]) )**2 )*-1*np.dot( directionVectors[j], x )
                
                j += 1
        
    xforce_block_list.append( np.sum(xforces) )
        
XFORCE = sum(xforce_block_list)            
print XFORCE