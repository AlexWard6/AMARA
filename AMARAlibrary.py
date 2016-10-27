# AMARAlibrary.py
# A. WARD 2016
"""
This file containts the generic functions for all vehicle parts. See individual
descriptions for an explanation.

This file gets executed by AMARAgeom.py
"""

import numpy as np

#---------------------------------  Generic ----------------------------------#



def sphericalpolar2car(R,theta,phi):
    x = R * np.sin(phi) * np.cos(theta) 
    y = R * np.sin(phi) * np.sin(theta)
    z = R * np.cos(phi)
    return x, y, z

def cylindricalpolar2car(R,theta,z):
    x = R * np.cos(theta) 
    y = R * np.sin(theta)
    z = z
    return x, y, z

def cart_to_2D_polar(x,y):
    r = math.sqrt( x*x + y*y)   
    theta = np.arctan2(y,x)

    return r, theta

def cart_to_cylindrical(x,y):
    r = math.sqrt( x*x + y*y)   
    theta = np.arctan2(y,x)

    return r, theta


def surf2surf(r, s, t, surfTop, surfBottom):
    """
    Function to create a block by interpolating between the two given surfaces.
	Inputs:
    r,s,t   - parametric positons within the volume
    surf1   - top surface
    surf2   - bottom surface
    
	Returns:
	x,y,z	- cartesian coordinates defining the block
    """
	
    # Use r and s to find coordinates along two surfaces
    topCoords = surfTop.eval(r,s)
    bottomCoords = surfBottom.eval(r,s)
    
    # Interpolate between the two surfaces with t
    x = (1 - t)*topCoords.x + t*bottomCoords.x
    y = (1 - t)*topCoords.y + t*bottomCoords.y
    z = (1 - t)*topCoords.z + t*bottomCoords.z
    
    return x,y,z
    
    
def evalFace(volume, faceSelection, rReversal=0, sReversal=0, rsReversal=0):
    if faceSelection == 'N':
        if rsReversal == 1:
            f = lambda r, s: (volume.eval((1.-r), 1., (1.-s)).x, 
                              volume.eval((1.-r), 1., (1.-s)).y, 
                              volume.eval((1.-r), 1., (1.-s)).z)
            return PyFunctionSurface(f, 'bottom')
        else:
            f = lambda r, s: (volume.eval(r, 1., s).x, 
                              volume.eval(r, 1., s).y, 
                              volume.eval(r, 1., s).z)
            return PyFunctionSurface(f, 'bottom')        
        
        f = lambda r, s: (volume.eval(r, 1., s).x, 
                          volume.eval(r, 1., s).y, 
                          volume.eval(r, 1., s).z)
        return PyFunctionSurface(f, 'north')
    
    elif faceSelection == 'E':
        f = lambda r, s: (volume.eval(1., r, s).x, 
                          volume.eval(1., r, s).y, 
                          volume.eval(1., r, s).z)
        return PyFunctionSurface(f, 'east')
    
    elif faceSelection == 'S':
        if rsReversal == 1:
            f = lambda r, s: (volume.eval((1.-r), 0., (1.-s)).x, 
                              volume.eval((1.-r), 0., (1.-s)).y, 
                              volume.eval((1.-r), 0., (1.-s)).z)
        else:        
            f = lambda r, s: (volume.eval(r, 0., s).x, 
                              volume.eval(r, 0., s).y, 
                              volume.eval(r, 0., s).z)
        return PyFunctionSurface(f, 'south')
    
    elif faceSelection == 'W':
        f = lambda r, s: (volume.eval(0., r, s).x, 
                          volume.eval(0., r, s).y, 
                          volume.eval(0., r, s).z)
        return PyFunctionSurface(f, 'west')

    elif faceSelection == 'T':
        f = lambda r, s: (volume.eval(r, s, 1.).x, 
                          volume.eval(r, s, 1.).y, 
                          volume.eval(r, s, 1.).z)
        return PyFunctionSurface(f, 'top')
        
    elif faceSelection == 'B':
        if rReversal == 1:
            f = lambda r, s: (volume.eval((1.-r), s, 0.).x, 
                              volume.eval((1.-r), s, 0.).y, 
                              volume.eval((1.-r), s, 0.).z)
            return PyFunctionSurface(f, 'bottom')
        else:
            f = lambda r, s: (volume.eval(r, s, 0.).x, 
                              volume.eval(r, s, 0.).y, 
                              volume.eval(r, s, 0.).z)
            return PyFunctionSurface(f, 'bottom')

    print "evalFace: face selection not recognised"
    return None

def blockSubdivide(r, s, t, block, division, index, direction):
    """ This function takes a block and subdivides in i, j, and k 
    according to iDivision, jDivision, kDivision.
    These define the number of resultant blocks in i, j and k directions 
    respectively (MAX 3)
    """
    
    if direction == 'i':
        coords = block.eval(index/division+r/division, s, t)
        x, y, z = coords.x, coords.y, coords.z
        return x, y, z

    elif direction == 'j':
        coords = block.eval(r, index/division+s/division, t)
        x, y, z = coords.x, coords.y, coords.z
        return x, y, z
        
    elif direction == 'k':
        coords = block.eval(r, s, index/division+t/division)
        x, y, z = coords.x, coords.y, coords.z
        return x, y, z
    
    print "Oops, block division didn't go as planned..."    
    return None

def distance_Node2Node(A, B):
    """ This function returns the straight line distance between two nodes, 
    A and B using the Pythagorean Theorem.
    """
    
    x = B.x - A.x; y = B.y - A.y; z = B.z - A.z  
    distance_2D = (x**2 + y**2)**0.5; distance = (distance_2D**2 + z**2)**0.5
    
    return distance
    
def get_Face(paths, face):
    """ This function builds a specific face given a set of 12 paths defining 
    a volume.
    """
    
    if face == 'N':
        # p76, p26, p32, p37
        return make_patch(paths[6], paths[10], paths[2], paths[11])
    elif face == 'E':
        # p56, p26, p12, p15
        return make_patch(paths[5], paths[10], paths[1], paths[9])
    elif face == 'S':
        # p45, p15, p01, p04
        return make_patch(paths[4], paths[9], paths[0], paths[8])        
    elif face == 'W':
        # p47, p37, p03, p04
        return make_patch(paths[7], paths[11], paths[3], paths[8])
    elif face == 'T':
        # p76, p56, p45, p47
        return make_patch(paths[6], paths[5], paths[4], paths[7])
    elif face == 'B':
        # p32, p12, p01, p03
        return make_patch(paths[2], paths[1], paths[0], paths[3])        
    print 'get_face failed'; return None

def findSide(zCoord, xCoord):
    """ This function finds if we're currently on the leeward or windward side 
    of the max thickness line. I think?
    """
    angle = np.arctan((halfspan - r_fb)/(c - c_t))
    adj = zCoord - zStartInlet - c_t
    opp = xCoord - r_fb
    side = np.tan(angle)*adj
    
    if xCoord < r_fb:
        return 'F'
    elif opp > side:
        return 'F'
    elif opp < side:
        return 'A'
    else:
        print 'findSide failed, sorry'
        return None
    
    return



#-------------------------------- Forebody -----------------------------------#

def forebodyLongitudinalCamber(z, C_a=C_a_long, C_p_fraction=0.0):
    """ This function returns the corrected y position given a camber amount
    Equation for the camber line of NACA 4 digit aerofoils
    
    C_a = maximum camber amount
    C_p = distance to max camber location specified as fraction of chord
    0. < C_p < 
    Assume small cambers, avoiding the need to redefine the coordinates normal to 
    the camber line.
    """
    c = abs(r_nose - z)
    chord = abs(L_ramp + r_nose)
    C_p = C_p_fraction*chord
    
    if c < C_p:
        y = C_a*(c/C_p**2.) * (2.*C_p - c/chord)    
    else:
        y = C_a*((chord - c)/(1. - C_p)**2.) * (1. + c/chord - 2.*C_p)
    
    return y

def forebodyTransverseCamber(x, y, z, zStart, C_a=C_a_trans):
    """ This function tadds some transvers camber to the forebody.
    
    """

    localRadius = AR*r_fb
    
    y = C_a*abs(x/localRadius)
    
    return y

def nose(r, s, t, rInner, dr, phis, AR, thetas, face='TOP', block=''):
    """ Function to create central square block positioned on Nose
    Inputs:
        r,s,t   - parametric positons within the volume
        r1      - Nose radius
        dr      - Thickness of the block in the readial direction
        phi1    - angle phi corresponding to N0, N1, N2, N3
        angle_N0, ...   - angular position of differnt points.
    """
    side = -1.
    face = -1. if face == 'BOTTOM' else 1.
    R = rInner + t*dr # t adjusts radius at which mesh is generated
    
    r_NE = R * np.sin(phis[0])
    r_SE = R * np.sin(phis[1])
    r_SW = R * np.sin(phis[2])
    r_NW = R * np.sin(phis[3])
    
    angle_NE = thetas[0]
    angle_SE = thetas[1]
    angle_SW = thetas[2]
    angle_NW = thetas[3]
    
    # Generate 2D cartesian grid defined by corner points
    NE = Node(side*r_NE * np.cos(angle_NE), face*r_NE * np.sin(angle_NE), 0., label="NE") 
    SE = Node(side*r_SE * np.cos(angle_SE), face*r_SE * np.sin(angle_SE), 0., label="SE") 
    SW = Node(side*r_SW * np.cos(angle_SW), face*r_SW * np.sin(angle_SW), 0., label="SW") 
    NW = Node(side*r_NW * np.cos(angle_NW), face*r_NW * np.sin(angle_NW), 0., label="NW")
    	# create Lines
    if side == 1:
        if face == -1:
            SWSE = Line(NW, NE) # N
            SENE = Line(NE, SE) # E
            NWNE = Line(SW, SE) # S
            SWNW = Line(NW,SW) # W
        else:
            NWNE = Line(NW, NE) # N
            SENE = Line(SE, NE) # E
            SWSE = Line(SW, SE) # S
            SWNW = Line(SW, NW) # W
    else:
        if face == -1:
            SWSE = Line(NE, NW) # N
            SWNW = Line(NE, SE) # E
            NWNE = Line(SE, SW) # S
            SENE = Line(NW, SW) # W
        else:
            NWNE = Line(NE, NW) # N
            SWNW = Line(SE, NE) # E
            SWSE = Line(SE, SW) # S
            SENE = Line(SW, NW) # W
    
    # create pach
    Face = make_patch(NWNE, SENE, SWSE, SWNW)
    
    # now use r and s to find locations in x,y coordinates of cartesian cell corners
    # Coords still has .x and .y attributes, it is a 'list' of vector objects
    Coords = Face.eval(r,s)
    
    # convert to 2-D polar coordinates
    # We now have a list of radii angles for the first cartesian list.
    r_2D, theta_2D = cart_to_2D_polar(Coords.x, Coords.y)
    
    # calculate phi
    phi = np.arcsin(r_2D / R)
    
    # convert from spherical polar coordinates to cartesian
    x,y,z = sphericalpolar2car(R, theta_2D, phi)
    y = y + forebodyLongitudinalCamber(z)
    y = y + forebodyTransverseCamber(x, y, z, r_nose)
    x = x*AR
    return x, y, z

def noseJoin(r, s, t, rInner, dr, phis, AR, thetas, camber=1, face='TOP', block=''):
    """ Function to create central square block positioned on Nose
    Inputs:
        r,s,t   - parametric positons within the volume
        r1      - Nose radius
        dr      - Thickness of the block in the readial direction
        phi1    - angle phi corresponding to N0, N1, N2, N3
        angle_N0, ...   - angular position of differnt points.
    """
    side = -1.
    face = -1. if face == 'BOTTOM' else 1.
    R = rInner + t*dr # t adjusts radius at which mesh is generated
    
    r_NE = R * np.sin(phis[0])
    r_SE = R * np.sin(phis[1])
    r_SW = R * np.sin(phis[2])
    r_NW = R * np.sin(phis[3])
    
    #print r_NE, r_SE, r_SW, r_NW    
    
    angle_NE = thetas[0]
    angle_SE = thetas[1]
    angle_SW = thetas[2]
    angle_NW = thetas[3]
    
    # Generate 2D cartesian grid defined by corner points
    NE = Node(side*r_NE * np.cos(angle_NE), face*r_NE * np.sin(angle_NE), 0., label="NE") 
    SE = Node(side*r_SE * np.cos(angle_SE), face*r_SE * np.sin(angle_SE), 0., label="SE") 
    SW = Node(side*r_SW * np.cos(angle_SW), face*r_SW * np.sin(angle_SW), 0., label="SW") 
    NW = Node(side*r_NW * np.cos(angle_NW), face*r_NW * np.sin(angle_NW), 0., label="NW")
    	# create Lines

    #print NE, SE, SW, NW
    #print block
    centre = Node(0., 0., 0.)
    if face == -1:   
        if block == '5': SWSE = Arc(NE, NW, centre)
        elif block == '6': SWSE = Arc(NE, NW, centre)
        else: SWSE = Line(NE, NW) # N
        if block == '3': SWNW = Arc(NE, SE, centre) 
        elif block == '4': SWNW = Arc(NE, SE, centre) 
        else: SWNW = Line(NE, SE) # E
        NWNE = Line(SE, SW) # S
        SENE = Line(NW, SW) # W    
    
    else:
        if block == '0': NWNE = Arc(NE, NW, centre)  
        elif block == '1': NWNE = Arc(NE, NW, centre)
        else: NWNE = Line(NE, NW)# N
        SWNW = Arc(SE, NE, centre) if block == '2' else Line(SE, NE) # E
        SWSE = Line(SE, SW) # S
        SENE = Line(SW, NW) # W
        
        
    
    # create pach
    #print NWNE, SENE, SWSE, SWNW
    Face = make_patch(NWNE, SENE, SWSE, SWNW)
    
    # now use r and s to find locations in x,y coordinates of cartesian cell corners
    # Coords still has .x and .y attributes, it is a 'list' of vector objects
    Coords = Face.eval(r,s)
    
    # convert to 2-D polar coordinates
    # We now have a list of radii angles for the first cartesian list.
    r_2D, theta_2D = cart_to_2D_polar(Coords.x, Coords.y)
    
    # calculate phi
    phi = np.arcsin(r_2D / R)
    
    # convert from spherical polar coordinates to cartesian
    x,y,z = sphericalpolar2car(R, theta_2D, phi)
    if camber != 0:
            y = y + forebodyLongitudinalCamber(z)
            y = y + forebodyTransverseCamber(x, y, z, r_nose)
    else:
            y = y + forebodyTransverseCamber(x, y, z, r_nose)
    x = x*AR
    return x, y, z

def ramp(r, s, t, bottomFace, face, rInnerStart, t_ramp, dr_ramp, phi_mesh, L_ramp, angle_In1, angle_In0, angle_Out1, angle_Out0, phi_join, block=''):
    """Function to create forebody ramp blocks.
    Inputs:
    		r,s,t      - parametric positons within the volume
    		r_inner    - Inner radius
           t_ramp     - Radial thickness of the ramp blocks
    		dr_ramp    - change in radius of ramp (join to base)
    		phi_mesh   - angle phi corresponding to In0, In1, In2, In3; inner nodes
    		L_ramp     - length of ramp
    		angle_In0, ...   - angular position of different points.
    """
    zStart = bottomFace.eval(0, 0).z
    zEnd_Inner = zStart - L_ramp
    zEnd_Outer = zEnd_Inner
    
    # Wall
    rInner = rInnerStart + dr_ramp
    # Mesh    
    rOuter = rInner + t_ramp + L_ramp*np.tan(phi_mesh)
    
    # Create 2D cartesian mesh for the top face
    topIn0 = Node(rInner * np.cos(angle_In0), rInner * np.sin(angle_In0), zEnd_Inner, label="bottomIn0") 
    topIn1 = Node(rInner * np.cos(angle_In1), rInner * np.sin(angle_In1), zEnd_Inner, label="bottomIn1") 
    topOut0 = Node(rOuter * np.cos(angle_Out0), rOuter * np.sin(angle_Out0), zEnd_Outer, label="bottomOut0") 
    topOut1 = Node(rOuter * np.cos(angle_Out1), rOuter * np.sin(angle_Out1), zEnd_Outer, label="bottomOut1")   

    if block == '2':
        topIn1 = Node(-1.*rInner, wingThickness(), zEnd_Inner, label="bottomIn0") 
        topOut1 = Node(-1.*rOuter, wingThickness(), zEnd_Outer, label="bottomOut1")   
    if block == '3':
        topIn1 = Node(-1.*rInner, -1.*wingThickness(), zEnd_Inner, label="bottomIn1") 
        topIn0 = Node(-1.*rInner, wingThickness(), zEnd_Inner, label="bottomIn0")
        topOut0 = Node(-1.*rOuter, wingThickness(), zEnd_Outer, label="bottomOut0") 
        topOut1 = Node(-1.*rOuter, -1.*wingThickness(), zEnd_Outer, label="bottomOut1") 
    if block == '4':
        topIn0 = Node(-1.*rInner, -1.*wingThickness(), zEnd_Inner, label="bottomIn1") 
        topOut0 = Node(-1.*rOuter, -1.*wingThickness(), zEnd_Outer, label="bottomOut0") 
        
    # create paths
    topOut0_Out1 = Arc(topOut1, topOut0, Node(0.,0.,topOut0.z)) if block == '' else Line(topOut1, topOut0)# S# N    
    topIn1_Out1 = Line(topIn0, topOut0) # E
    topIn0_In1 = Arc(topIn1, topIn0, Node(0.,0.,topIn0.z)) if block == '' else Line(topIn1, topIn0)# S
    topIn0_Out0 = Line(topIn1, topOut1) # w

    # create patch (pN, pE, pS, pW)
    topFace = make_patch(topOut0_Out1, topIn1_Out1,
                         topIn0_In1, topIn0_Out0)

    AR=1
    # This function is to correct for some of the mating faces being S not N 
    # (this results in negative volume)
    if face == 'N':
            bottomCoords = topFace.eval(r,s)
            topCoords = bottomFace.eval(r,s)
            x = (1 - t)*topCoords.x + t*AR*bottomCoords.x
            z = (1 - t)*topCoords.z + t*bottomCoords.z
            y = (1 - t)*topCoords.y + t*bottomCoords.y + forebodyLongitudinalCamber(z) 
            y = y + (t)*forebodyTransverseCamber(x, y, z, zStart)
    else:
            topCoords = topFace.eval(r,s)
            bottomCoords = bottomFace.eval(r,s)
            x = (1 - t)*AR*topCoords.x + t*bottomCoords.x
            z = (1 - t)*topCoords.z + t*bottomCoords.z
            y = (1 - t)*topCoords.y + t*bottomCoords.y + forebodyLongitudinalCamber(z) 
            y = y + (1-t)*forebodyTransverseCamber(x, y, z, zStart)
    
    return x,y,z
    
#------------------------------- Centrebody ----------------------------------#
class rootFairing:
    """ This function builds the block right coincident with the inlet, or 
    where the wing is faired into tho body.
    
    BLOCK 0 according to logbook sketch.
    
    The aft node is halfway between the mid nodes and gives block boundaries 
    of 120 deg.
    """
    def __init__(self, length, zStart, xStart, angleWall, angleWing, angle, side='RHS', face='TOP'):
        if side == 'LHS':
            side = 1.
        elif side == 'RHS':
            side = -1.
            self.side = side
        else:
            print 'Incorrect wing side choice'; return None
        
        # Create block
        if face == 'TOP':
            self.topsideBlock(length, zStart, xStart, angleWall, angleWing, angle, side=side)
        elif face == 'BOTTOM':
            self.undersideBlock(length, zStart, xStart, angleWall, angleWing, angle, side=side)
        
    def topsideBlock(self, length, zStart, xStart, angleWall, angleWing, angle, side='RHS'):
        fore_flat = Node(side*xStart, 0, zStart, label='fore')
        midOut_flat = Node(fore_flat.x+length/np.tan(angleWing)*side, 0, fore_flat.z-length, label='midOut')
        midIn_flat = Node(fore_flat.x+length*np.tan(angleWall)*side, 0, fore_flat.z-length, label='midIn')
        aft_flat = Node(midIn_flat.x + (midOut_flat.x-midIn_flat.x)/2, 0, fore_flat.z-length-(abs(midOut_flat.x-midIn_flat.x)/2/np.tan(60.*np.pi/180)), label='aft')    
        
        pN_flat = Line(fore_flat, midOut_flat)     # p32
        pE_flat = Line(aft_flat, midOut_flat)      # p12
        pS_flat = Line(midIn_flat, aft_flat)       # p01
        pW_flat = Line(midIn_flat, fore_flat)  # p03
        
        fore = Node(side*xStart, wingThickness(), zStart, label='fore')
        midOut = Node(fore.x+length/np.tan(angleWing)*side, wingThickness(), fore.z-length, label='midOut')
        midIn = Node(fore.x+length*np.tan(angleWall)*side, wingThickness(), fore.z-length, label='midIn')
        aft = Node(midIn.x + (midOut.x-midIn.x)/2, wingThickness(), fore.z-length-(abs(midOut.x-midIn.x)/2/np.tan(60.*np.pi/180)), label='aft')            
        
        pN = Line(fore, midOut)     # p32
        pE = Line(aft, midOut)      # p12
        pS = Line(midIn, aft)       # p01
        pW = Line(midIn, fore)  # p03        
        
        pNtop = pN_flat.clone().rotate_about_zaxis(side*angle)
        pEtop = pE_flat.clone().rotate_about_zaxis(side*angle)
        pStop = pS_flat.clone().rotate_about_zaxis(side*angle)
        pWtop = pW_flat.clone().rotate_about_zaxis(side*angle)
    
        #foreArc = Arc(fore, pNtop.eval(0.), Node(0.,0.,fore.z))
        foreArc = Line(fore, pNtop.eval(0.))
        #midOutArc = Arc(midOut, pNtop.eval(1.), Node(0.,0.,midOut.z))
        midOutLine = Line(midOut, pNtop.eval(1.))
        #midInArc = Arc(midIn, pStop.eval(0.), Node(0.,0.,midIn.z))
        midInArc = Line(midIn, pStop.eval(0.))
        aftArc = Line(aft, pStop.eval(1.))
        
        #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37    
        paths = [pStop, pEtop, pNtop, pWtop, 
                 pS, pE, pN, pW, 
                 midInArc, aftArc, midOutLine, foreArc]
        
        volume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                 paths[4], paths[5], paths[6], paths[7],
                                 paths[8], paths[9], paths[10], paths[11])   
        
        southFace = make_patch(pStop, aftArc, pS, midInArc) # N E S W
        eastFace = make_patch(pEtop, midOutLine, pE, aftArc)
        
        #self.side = side
        self.volume = volume
        self.southFace = southFace
        self.eastFace = eastFace
        self.surfacePath0 = pN
        self.surfacePath_outer0 = pNtop
        self.TopPaths =  [pNtop, pEtop, pStop, pWtop]
        self.flatPaths =  [pN_flat, pE_flat, pS_flat, pW_flat]
        #self.rOuterEnd = fore_flat.x+length/np.tan(angleWing)*self.side
        

    def undersideBlock(self, length, zStart, xStart, angleWall, angleWing, angle, side='RHS', block='AO'):
        fore_flat = Node(side*xStart, 0, zStart, label='fore')
        midOut_flat = Node(fore_flat.x+length/np.tan(angleWing)*side, 0, fore_flat.z-length, label='midOut')
        midIn_flat = Node(fore_flat.x+length*np.tan(angleWall)*side, 0, fore_flat.z-length, label='midIn')
        aft_flat = Node(midIn_flat.x + (midOut_flat.x-midIn_flat.x)/2, 0, fore_flat.z-length-(abs(midOut_flat.x-midIn_flat.x)/2/np.tan(60.*np.pi/180)), label='aft')    
        
        pN_flat = Line(fore_flat, midOut_flat)     # p32
        pE_flat = Line(aft_flat, midOut_flat)      # p12
        pS_flat = Line(midIn_flat, aft_flat)       # p01
        pW_flat = Line(midIn_flat, fore_flat)  # p03
        
        fore = Node(side*xStart, -1*wingThickness(), zStart, label='fore')
        midOut = Node(fore.x+length/np.tan(angleWing)*side, -1*wingThickness(), fore.z-length, label='midOut')
        midIn = Node(fore.x+length*np.tan(angleWall)*side, -1*wingThickness(), fore.z-length, label='midIn')
        aft = Node(midIn.x + (midOut.x-midIn.x)/2, -1*wingThickness(), fore.z-length-(abs(midOut.x-midIn.x)/2/np.tan(60.*np.pi/180)), label='aft')            
        
        pN = Line(fore, midOut)     # p32
        pE = Line(aft, midOut)      # p12
        pS = Line(midIn, aft)       # p01
        pW = Line(midIn, fore)  # p03        
        
        pNtop = pN_flat.clone().rotate_about_zaxis(side*angle)
        pEtop = pE_flat.clone().rotate_about_zaxis(side*angle)
        pStop = pS_flat.clone().rotate_about_zaxis(side*angle)
        pWtop = pW_flat.clone().rotate_about_zaxis(side*angle)
    
#        foreArc = Line(fore, pNtop.eval(0.))
#        midOutLine = Line(midOut, pNtop.eval(1.))
#        midInArc = Line(midIn, pStop.eval(0.))
#        aftArc = Line(aft, pStop.eval(1.))
        
        foreArc = Line(pNtop.eval(0.), fore)
        midOutLine = Line(pNtop.eval(1.), midOut)
        midInArc = Line(pStop.eval(0.), midIn)
        aftArc = Line(pStop.eval(1.), aft)        
                 
        paths = [pStop, pEtop, pNtop, pWtop, 
                 pS, pE, pN, pW, 
                 midInArc, aftArc, midOutLine, foreArc]
        
#        volume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
#                                 paths[4], paths[5], paths[6], paths[7],
#                                 paths[8], paths[9], paths[10], paths[11])   
        
        
        southFace = make_patch(pStop, aftArc.reverse(), pS, midInArc.reverse()) # N E S W
        eastFace = make_patch(pEtop, midOutLine.reverse(), pE, aftArc)

        # AOpatch(pS, pN , pW , pE)
        self.nx = 5; self.ny = 5
        Nsurf = AOPatch(paths[2], paths[6], paths[11], paths[10].clone().reverse(), 'N', self.nx, self.ny)
        Esurf = AOPatch(paths[1], paths[5], paths[9].clone().reverse(), paths[10].clone().reverse(), 'E', self.nx, self.ny)
        Ssurf = AOPatch(paths[0], paths[4], paths[8].clone().reverse(), paths[9].clone().reverse(), 'S', self.nx, self.ny)
        Wsurf = AOPatch(paths[3], paths[7], paths[8].clone().reverse(), paths[11], 'W', self.nx, self.ny)
        Tsurf = AOPatch(paths[4], paths[6], paths[7], paths[5], 'T', self.nx, self.ny)
        Bsurf = AOPatch(paths[0], paths[2], paths[3], paths[1], 'B', self.nx, self.ny)
        
        volume = ParametricVolume(Nsurf, Esurf, Ssurf, Wsurf, Tsurf, Bsurf)
        
        self.volume = volume
        self.southFace = southFace
        self.eastFace = eastFace
        self.surfacePath0 = pN
        self.surfacePath_outer0 = pNtop
        self.TopPaths =  [pNtop, pEtop, pStop, pWtop]
        self.flatPaths =  [pN_flat, pE_flat, pS_flat, pW_flat]
        self.rOuterEnd = midOut_flat.x*self.side
        #self.side = side
    
    def rotate(self, angleStart, angleEnd, flatPaths, face='TOP'):
        rotation = angleEnd - angleStart
        
        pNflat = flatPaths[0]; pEflat = flatPaths[1];
        pSflat = flatPaths[2]; pWflat = flatPaths[3]
        
        pN = pNflat.clone().rotate_about_zaxis(angleStart)
        pE = pEflat.clone().rotate_about_zaxis(angleStart)
        pS = pSflat.clone().rotate_about_zaxis(angleStart)
        pW = pWflat.clone().rotate_about_zaxis(angleStart)
        
        pNtop = pN.clone().rotate_about_zaxis(rotation)
        pEtop = pE.clone().rotate_about_zaxis(rotation)
        pStop = pS.clone().rotate_about_zaxis(rotation)
        pWtop = pW.clone().rotate_about_zaxis(rotation)

        foreArc = Arc(pN.eval(0), pNtop.eval(0.), Node(0.,0.,pN.eval(0).z))
        midOutArc = Arc(pN.eval(1), pNtop.eval(1.), Node(0.,0.,pN.eval(1).z))
        midInArc = Arc(pS.eval(0), pStop.eval(0.), Node(0.,0.,pS.eval(0).z))
        aftArc = Arc(pS.eval(1), pStop.eval(1.), Node(0.,0.,pS.eval(1).z))
        

        paths = [pS, pE, pN, pW, 
                 pStop, pEtop, pNtop, pWtop, 
                 midInArc, aftArc, midOutArc, foreArc]
        
        if face == 'TOP':
            #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
            rotatedVolume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                            paths[4], paths[5], paths[6], paths[7],
                                            paths[8], paths[9], paths[10], paths[11])  
            return rotatedVolume
            
        else: # BOTTOM
            #pS, pN, pW, pE
            Nsurf = AOPatch(paths[2], paths[6], paths[11], paths[10], 'N', self.nx, self.ny)
            Esurf = AOPatch(paths[1], paths[5], paths[9], paths[10], 'E', self.nx, self.ny)
            Ssurf = AOPatch(paths[0], paths[4], paths[8], paths[9], 'S', self.nx, self.ny)
            Wsurf = AOPatch(paths[3], paths[7], paths[8], paths[11], 'W', self.nx, self.ny)
            Tsurf = AOPatch(paths[4], paths[6], paths[7], paths[5], 'T', self.nx, self.ny)
            Bsurf = AOPatch(paths[0], paths[2], paths[3], paths[1], 'B', self.nx, self.ny)
            
            rotatedVolume = ParametricVolume(Nsurf, Esurf, Ssurf, Wsurf, Tsurf, Bsurf)

            return rotatedVolume, Ssurf, Esurf
         

class foreDiamonds:
    def __init__(self, zEnd, block, endAngle, side, startAngle=0, face='TOP'):
        """ This function builds the rectangular blocks fore of the diamond blocks 
        which sit centrally on the combustor.
        
        Blocks extrude horizontally, fore to aft.
        """
        self.endAngle = endAngle
        self.startAngle = startAngle
        self.zEnd = zEnd
        if side == 'LHS':
            self.side = 1.
        elif side == 'RHS':
            self.side = -1.
        else:
            print 'Incorrect wing side choice'; return None
        
        # This angle defines the middle between the forebody wall and the wing sweep.
        self.phi_mean = ((90.*np.pi/180-phi_sweep) + (phi_fb))/2
       
    def foreDiamonds_1(self, matingFace):

        foreIn = matingFace.eval(0,0)
        foreOut = matingFace.eval(1,0)
        aftOut = Node(foreOut.x+self.side*abs((self.zEnd-foreOut.z)*np.tan(self.phi_mean)), foreOut.y, self.zEnd)
        aftIn = Node(self.side*r_fb, foreOut.y, self.zEnd)

        pN = Line(foreIn, foreOut)
        pE = Line(aftOut, foreOut)
        pS = Line(aftIn, aftOut)
        pW = Line(aftIn, foreIn)

        foreIn_flat = Node(foreIn.x, 0., foreIn.z)
        foreOut_flat = Node(foreOut.x, 0., foreOut.z)
        aftOut_flat = Node(foreOut.x+self.side*abs((self.zEnd-foreOut.z)*np.tan(self.phi_mean)), 0., self.zEnd)
        aftIn_flat = Node(self.side*r_fb, 0., self.zEnd)        

        pN_flat = Line(foreIn_flat, foreOut_flat)
        pE_flat = Line(aftOut_flat, foreOut_flat)
        pS_flat = Line(aftIn_flat, aftOut_flat)
        pW_flat = Line(aftIn_flat, foreIn_flat)
        #print 'pN_flat', pN_flat, pE_flat, 'pE_flat',  pS_flat, 'pS_flat', pW_flat, 'pW_flat'
        
        pNtop = pN_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pEtop = pE_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pStop = pS_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pWtop = pW_flat.clone().rotate_about_zaxis(self.side*self.endAngle)

#        foreInArc = Arc(pNtop.eval(0.), foreIn, Node(0.,0.,foreIn.z))
#        foreOutArc = Arc(pNtop.eval(1.), foreOut, Node(0.,0.,foreOut.z))
#        aftInArc = Arc(pStop.eval(0.), aftIn, Node(0.,0.,aftIn.z))
#        aftOutArc = Arc(pStop.eval(1.), aftOut, Node(0.,0.,aftOut.z))        
#    
#        foreInArc = Arc(foreIn, pNtop.eval(0.), Node(0.,0.,foreIn.z))
#        foreOutArc = Arc(foreOut, pNtop.eval(1.), Node(0.,0.,foreOut.z))
#        aftInArc = Arc(aftIn, pStop.eval(0.), Node(0.,0.,aftIn.z))
#        aftOutArc = Arc(aftOut, pStop.eval(1.), Node(0.,0.,aftOut.z))
        
        foreInArc = Line(foreIn, pNtop.eval(0.))
        foreOutArc = Line(foreOut, pNtop.eval(1.))
        aftInArc = Line(aftIn, pStop.eval(0.))
        aftOutArc = Line(aftOut, pStop.eval(1.))
        
        #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        paths = [pS, pE, pN, pW, pStop, pEtop, pNtop, pWtop, aftInArc, aftOutArc, foreOutArc, foreInArc]
        volume1 = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                 paths[4], paths[5], paths[6], paths[7],
                                 paths[8], paths[9], paths[10], paths[11])          
        #print 'LHS VOLUME 1', volume1.eval(0,0,0), volume1.eval(0,1,0), volume1.eval(1,0,0), volume1.eval(1,1,0)
#       print pS, pS.reverse(), pS.clone().reverse()
        self.combustorTop3 = make_patch(pStop, aftOutArc, pS, aftInArc)
        self.TopPaths1 =  [pNtop, pEtop, pStop, pWtop]
        self.flatPaths1 =  [pN_flat, pE_flat, pS_flat, pW_flat]
        return volume1

    def foreDiamonds_2(self, matingFace):      
#        foreIn = matingFace.eval(0,0)
#        foreOut = matingFace.eval(1,0)       
#        aftIn = Node(foreIn.x+self.side*abs((self.zEnd-foreIn.z)*np.tan(self.phi_mean)), foreIn.y, self.zEnd)
#        aftOut = Node(foreOut.x+self.side*abs((self.zEnd-foreOut.z)/np.tan(phi_sweep)), foreOut.y, self.zEnd)
        #widthLocal = halfspan - abs(zStartboatTail - aftIn.z)/np.tan(phi_sweep) - r_fb
        
        foreIn = matingFace.eval(0,0)
        foreOut = matingFace.eval(1,0)       
        aftIn = Node(foreIn.x+self.side*abs((self.zEnd-foreIn.z)*np.tan(self.phi_mean)), foreIn.y, self.zEnd)
        #aftOut = Node(foreOut.x+self.side*abs((self.zEnd-foreOut.z)/np.tan(phi_sweep)), foreOut.y, self.zEnd)    
        widthLocal = halfspan - abs(zStartboatTail - aftIn.z)/np.tan(phi_sweep)
        aftOut = Node(self.side*widthLocal, foreOut.y, self.zEnd)    
        
        pN = Line(foreIn, foreOut)
        pE = Line(aftOut, foreOut)
        pS = Line(aftIn, aftOut)
        pW = Line(aftIn, foreIn)            
        
        foreIn_flat = Node(foreIn.x, 0., foreIn.z)
        foreOut_flat = Node(foreOut.x, 0., foreOut.z)
        aftIn_flat = Node(foreIn.x+self.side*abs((self.zEnd-foreIn.z)*np.tan(self.phi_mean)), 0., self.zEnd)
        aftOut_flat = Node(self.side*widthLocal, 0., self.zEnd)    
        
        pN_flat = Line(foreIn_flat, foreOut_flat)
        pE_flat = Line(aftOut_flat, foreOut_flat)
        pS_flat = Line(aftIn_flat, aftOut_flat)
        pW_flat = Line(aftIn_flat, foreIn_flat)            
        
        #print 'pN_flat', pN_flat, pE_flat, 'pE_flat',  pS_flat, 'pS_flat', pW_flat, 'pW_flat'
        
        pNtop = pN_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pEtop = pE_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pStop = pS_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pWtop = pW_flat.clone().rotate_about_zaxis(self.side*self.endAngle)

#        foreInArc = Arc(foreIn, pNtop.eval(0.), Node(0.,0.,foreIn.z))
#        foreOutArc = Line(foreOut, pNtop.eval(1.))
#        aftInArc = Arc(aftIn, pStop.eval(0.), Node(0.,0.,aftIn.z))
#        aftOutArc = Line(aftOut, pStop.eval(1.))
        
        foreInArc = Line(foreIn, pNtop.eval(0.))
        foreOutArc = Line(foreOut, pNtop.eval(1.))
        aftInArc = Line(aftIn, pStop.eval(0.))
        aftOutArc = Line(aftOut, pStop.eval(1.))
        
        #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        volume2 = WireFrameVolume(pS, pE, pN, pW, 
                                  pStop, pEtop, pNtop, pWtop, 
                                  aftInArc, aftOutArc, foreOutArc, foreInArc)
        
        self.surfacePath2 = pE
        self.surfacePath_outer2 = pEtop
        self.combustorTop4 = make_patch(pStop, aftOutArc, pS, aftInArc)
        self.TopPaths2 =  [pNtop, pEtop, pStop, pWtop]
        self.flatPaths2 =  [pN_flat, pE_flat, pS_flat, pW_flat]
        self.rOuterEnd = aftOut_flat.x*self.side
        return volume2

    def rotate(self, angleStart, angleEnd, flatPaths):
        rotation = angleEnd - angleStart
        
        pNflat = flatPaths[0]; pEflat = flatPaths[1];
        pSflat = flatPaths[2]; pWflat = flatPaths[3]
        
        pN = pNflat.clone().rotate_about_zaxis(angleStart)
        pE = pEflat.clone().rotate_about_zaxis(angleStart)
        pS = pSflat.clone().rotate_about_zaxis(angleStart)
        pW = pWflat.clone().rotate_about_zaxis(angleStart)
        
        pNtop = pN.clone().rotate_about_zaxis(rotation)
        pEtop = pE.clone().rotate_about_zaxis(rotation)
        pStop = pS.clone().rotate_about_zaxis(rotation)
        pWtop = pW.clone().rotate_about_zaxis(rotation)

        foreInArc = Arc(pN.eval(0), pNtop.eval(0.), Node(0.,0.,pN.eval(0).z))
        foreOutArc = Arc(pN.eval(1), pNtop.eval(1.), Node(0.,0.,pN.eval(1).z))
        aftInArc = Arc(pS.eval(0), pStop.eval(0.), Node(0.,0.,pS.eval(0).z))
        aftOutArc = Arc(pS.eval(1), pStop.eval(1.), Node(0.,0.,pS.eval(1).z))
        

        paths = [pS, pE, pN, pW, 
                 pStop, pEtop, pNtop, pWtop, 
                 aftInArc, aftOutArc, foreOutArc, foreInArc]
        
        #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        rotatedVolume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                        paths[4], paths[5], paths[6], paths[7],
                                        paths[8], paths[9], paths[10], paths[11])  
        
        #self.surfacePath2 = pE
        #self.surfacePath_outer2 = pEtop
        self.combustorTop = make_patch(pStop, aftOutArc, pS, aftInArc)
        return rotatedVolume, self.combustorTop

class diamonds:
    """ This function builds the diamond blocks which sit centrally on the 
    wing.
    
    If we rotate the diamond a little bit, the quality of block 4 is improved 
    (to the detriment of block 3 but still). Denominator in diamondRotation 
    changes angle relatively easily. 1.6 = horizontal, 1.75 seems to work 
    quite well. 2 gives the mean angle between wing and forebody
    
    The block on the LE is the worst, changing diamondRotation changes the 
    quality of predominantly blocks 3 and 4 (also pretty bad).
    """

    def __init__(self, endAngle, diamondRotation, side='RHS', startAngle=0, face='TOP'):
        if side == 'LHS':
            self.side = 1.
        elif side == 'RHS':
            self.side = -1.
        else:
            print 'Incorrect wing side choice'; return None
        
        self.face = -1 if face == 'BOTTOM' else 1
        self.endAngle = endAngle
        self.startAngle = startAngle
        

        if self.face == -1:
            # Underside
            self.t0 = Node(self.side*(rInner_combustor + 0.5*L_combustor*np.tan(wallAngle_combustor)), self.face*wingThickness(), (zStartCombustor - 0.5*L_combustor))
            self.t0_flat = Node(self.side*(rInner_combustor + 0.5*L_combustor*np.tan(wallAngle_combustor)), 0., (zStartCombustor - 0.5*L_combustor))            
            widthLocal = halfspan - abs(zStartboatTail - self.t0.z)/np.tan(phi_sweep) - rInner_combustor + 0.5*L_combustor*np.tan(wallAngle_combustor)
            #diamondRotation = 5.*np.pi/180.
        
        else:            
            # topside
            self.t0 = Node(self.side*r_fb, self.face*wingThickness(), (zStartCombustor - 0.5*L_combustor))
            self.t0_flat = Node(self.side*r_fb, 0., (zStartCombustor - 0.5*L_combustor))
            widthLocal = halfspan - abs(zStartboatTail - self.t0.z)/np.tan(phi_sweep) - r_fb
            #diamondRotation = (90.*np.pi/180-phi_sweep + phi_fb)/2.0
            #diamondRotation = 5.*np.pi/180
        
        #print 'diamondRotation', diamondRotation

        #dz = np.sin(90.*np.pi/180. - diamondRotation)*np.cos(90.*np.pi/180. - diamondRotation)*widthLocal
        tempLength = np.sin(diamondRotation)/np.sin(np.pi - diamondRotation - phi_sweep)*widthLocal
        dz = abs(np.cos(90.*np.pi/180 - phi_sweep)*tempLength)
        dx = (self.side*(widthLocal - np.cos(phi_sweep)*tempLength))
        
        self.t2 = Node(self.t0.x + dx, self.face*wingThickness(), self.t0.z + dz)
        self.t2_flat = Node(self.t0.x + dx, 0., self.t0.z + dz)
        
        
        self.centre = Line(self.t0, self.t2).eval(0.45 if self.face == -1 else 0.5)
        self.centre_flat = Line(self.t0_flat, self.t2_flat).eval(0.45 if self.face == -1 else 0.5)
        centreLength = distance_Node2Node(self.t0, self.centre)
        
        t3_temp = Node(self.t0.x + self.side*centreLength, 
                       self.face*wingThickness(), 
                       self.t0.z + centreLength*np.tan(30.*np.pi/180.))
        
        t3x = np.cos(-1.*(diamondRotation))*(t3_temp.x-self.t0.x) - np.sin(-1.*(diamondRotation))*(t3_temp.z-self.t0.z) + self.t0.x
        t3z = np.sin(-1.*(diamondRotation))*(t3_temp.x-self.t0.x) + np.cos(-1.*(diamondRotation))*(t3_temp.z-self.t0.z) + self.t0.z
        
        self.t3 = Node(t3x, self.face*wingThickness(), t3z)
        self.t3_flat = Node(t3x, 0., t3z)
        
        t1_temp = Node(self.t0.x + self.side*centreLength, 
                       self.face*wingThickness(), 
                       self.t0.z - centreLength*np.tan(30.*np.pi/180.))
        
        t1x = np.cos(-1.*(diamondRotation))*(t1_temp.x-self.t0.x) - np.sin(-1.*(diamondRotation))*(t1_temp.z-self.t0.z) + self.t0.x
        t1z = np.sin(-1.*(diamondRotation))*(t1_temp.x-self.t0.x) + np.cos(-1.*(diamondRotation))*(t1_temp.z-self.t0.z) + self.t0.z
        
        self.t1 = Node(t1x, self.face*wingThickness(), t1z)   
        self.t1_flat = Node(t1x, 0., t1z)   
        
        #print 't0 =', self.t0, 't2 =', self.t2, 't3_temp =', t3_temp, 't3 =', self.t3, 't1_temp', t1_temp, 't1 =', self.t1
        #print 'widthLocal', widthLocal, 'centre', self.centre
        
        self.pN = Line(self.t3, self.t2);  self.pE = Line(self.t1, self.t2)
        self.pS = Line(self.t0, self.t1);  self.pW = Line(self.t0, self.t3)

        self.pN_flat = Line(self.t3_flat, self.t2_flat);  self.pE_flat = Line(self.t1_flat, self.t2_flat)
        self.pS_flat = Line(self.t0_flat, self.t1_flat);  self.pW_flat = Line(self.t0_flat, self.t3_flat)
        
        self.rOuter = self.t2_flat.x*self.side
        self.rInner = self.t0_flat.x*self.side
        self.Nodetwoz = self.t2_flat.z
    
    def diamond5(self):        
        pN_5_flat = Line(self.t3_flat, self.centre_flat)
        pE_5_flat = Line(self.pS_flat.eval(0.5), self.centre_flat)
        pS_5_flat = Line(self.t0_flat, self.pS_flat.eval(0.5))
        pW_5_flat = self.pW_flat           
        
        pN_5 = Line(self.t3, self.centre)
        pE_5 = Line(self.pS.eval(0.5), self.centre)
        pS_5 = Line(self.t0, self.pS.eval(0.5))
        pW_5 = self.pW    
        #print 'pN', pN_5, 'pE', pE_5, 'pS', pS_5, 'pW', pW_5
    
        pNtop = pN_5_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pEtop = pE_5_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pStop = pS_5_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pWtop = pW_5_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        #print 'pNtop', pNtop, 'pEtop', pEtop, 'pStop', pStop, 'pWtop', pWtop
    
#        foreArc = Arc(self.t3, pNtop.eval(0.), Node(0.,0.,self.t3.z))
#        midInArc = Arc(self.t0, pStop.eval(0.), Node(0.,0.,self.t0.z))
#        midOutArc = Arc(self.centre, pNtop.eval(1.), Node(0.,0.,self.centre.z))
#        aftArc = Arc(self.pS.eval(0.5), pStop.eval(1.), Node(0.,0.,self.pS.eval(0.5).z))
        
        if self.face == 1:
            foreArc = Line(self.t3, pNtop.eval(0.))
            midInArc = Line(self.t0, pStop.eval(0.))
            midOutArc = Line(self.centre, pNtop.eval(1.))
            aftArc = Line(self.pS.eval(0.5), pStop.eval(1.))
        else:
            foreArc = Line(pNtop.eval(0.), self.t3)
            midInArc = Line(pStop.eval(0.), self.t0)
            midOutArc = Line(pNtop.eval(1.), self.centre)
            aftArc = Line(pStop.eval(1.), self.pS.eval(0.5))
            
        #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        paths = [pS_5, pE_5, pN_5, pW_5, pStop, pEtop, pNtop, pWtop, midInArc, aftArc, midOutArc, foreArc]
        volume5 = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                  paths[4], paths[5], paths[6], paths[7],
                                  paths[8], paths[9], paths[10], paths[11])   

        #print pW_5, foreArc, pWtop, midInArc         
        self.westFace5 = get_Face(paths, 'W') if self.face == 1 else make_patch(pWtop, foreArc.reverse(), pW_5, midInArc.reverse()) #make_patch(pWtop, foreArc, pW_5, midInArc)
        self.southFace5 = get_Face(paths, 'S') if self.face == 1 else make_patch(pStop, aftArc.reverse(), pS_5, midInArc)
        self.flatPaths5 =  [pN_5_flat, pE_5_flat, pS_5_flat, pW_5_flat]
        return volume5

    def diamond6(self):
        pN_6_flat = self.pN_flat
        pE_6_flat = Line(self.pE_flat.eval(0.5), self.t2_flat)
        pS_6_flat = Line(self.centre_flat, self.pE_flat.eval(0.5))
        pW_6_flat = Line(self.centre_flat, self.t3_flat)

        pN_6 = self.pN
        pE_6 = Line(self.pE.eval(0.5), self.t2)
        pS_6 = Line(self.centre, self.pE.eval(0.5))
        pW_6 = Line(self.centre, self.t3)  
        #print 'pN', pN, 'pE', pE, 'pS', pS, 'pW', pW    
        
        pNtop = pN_6_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pEtop = pE_6_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pStop = pS_6_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pWtop = pW_6_flat.clone().rotate_about_zaxis(self.side*self.endAngle)    
        #print 'pNtop', pNtop, 'pEtop', pEtop, 'pStop', pStop, 'pWtop', pWtop    
        
#        foreArc = Arc(self.t3, pNtop.eval(0.), Node(0.,0.,self.t3.z))
#        midInArc = Arc(self.centre, pWtop.eval(0.), Node(0.,0.,self.centre.z))
#        midOutArc = Line(self.t2, pNtop.eval(1.))
#        aftArc = Arc(pE_6.eval(0.), pStop.eval(1.), Node(0.,0.,pE_6.eval(0.).z))

        if self.face == 1:
            foreArc = Line(self.t3, pNtop.eval(0.))
            midInArc = Line(self.centre, pWtop.eval(0.))
            midOutArc = Line(self.t2, pNtop.eval(1.))
            aftArc = Line(pE_6.eval(0.), pStop.eval(1.))
        else:
            foreArc = Line(pNtop.eval(0.), self.t3)
            midInArc = Line(pWtop.eval(0.), self.centre)
            midOutArc = Line(pNtop.eval(1.), self.t2)
            aftArc = Line(pStop.eval(1.), pE_6.eval(0.))
            
        #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        paths = [pS_6, pE_6, pN_6, pW_6, pStop, pEtop, pNtop, pWtop, midInArc, aftArc, midOutArc, foreArc]
        volume6 = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                  paths[4], paths[5], paths[6], paths[7],
                                  paths[8], paths[9], paths[10], paths[11])

        self.northFace6 = get_Face(paths, 'N') if self.face == 1 else make_patch(pNtop, midOutArc.reverse(), pN_6, foreArc.reverse())# make_patch(pNtop, midOutArc.reverse(), pN_6, foreArc.reverse()) # N E S W
        self.eastFace6 = get_Face(paths, 'E') if self.face == 1 else make_patch(pEtop, midOutArc, pE_6, aftArc.reverse())
        self.flatPaths6 = [pN_6_flat, pE_6_flat, pS_6_flat, pW_6_flat]
        return volume6
        
    def diamond7(self):
        pN_7_flat = Line(self.centre_flat, self.pE_flat.eval(0.5))
        pE_7_flat = Line(self.t1_flat, self.pE_flat.eval(0.5))
        pS_7_flat = Line(self.pS_flat.eval(0.5), self.t1_flat)
        pW_7_flat = Line(self.pS_flat.eval(0.5), self.centre_flat)
        
        pN_7 = Line(self.centre, self.pE.eval(0.5))
        pE_7 = Line(self.t1, self.pE.eval(0.5))
        pS_7 = Line(self.pS.eval(0.5), self.t1)
        pW_7 = Line(self.pS.eval(0.5), self.centre)
        #print 'pN', pN, 'pE', pE, 'pS', pS, 'pW', pW    
        
        pNtop = pN_7_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pEtop = pE_7_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pStop = pS_7_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pWtop = pW_7_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        #print 'pNtop', pNtop, 'pEtop', pEtop, 'pStop', pStop, 'pWtop', pWtop    
        
#        foreArc = Arc(self.centre, pNtop.eval(0.), Node(0.,0.,self.centre.z))
#        midInArc = Arc(pS_7.eval(0.), pStop.eval(0.), Node(0.,0.,pS_7.eval(0.).z))
#        midOutArc = Arc(pE_7.eval(1.), pEtop.eval(1.), Node(0.,0.,pE_7.eval(1.).z))
#        aftArc = Arc(self.t1, pStop.eval(1.), Node(0.,0.,self.t1.z))

        if self.face == 1:
            foreArc = Line(self.centre, pNtop.eval(0.))
            midInArc = Line(pS_7.eval(0.), pStop.eval(0.))
            midOutArc = Line(pE_7.eval(1.), pEtop.eval(1.))
            aftArc = Line(self.t1, pStop.eval(1.)) 
        else:
            foreArc = Line(pNtop.eval(0.),self.centre)
            midInArc = Line(pStop.eval(0.), pS_7.eval(0.))
            midOutArc = Line(pEtop.eval(1.), pE_7.eval(1.))
            aftArc = Line(pStop.eval(1.), self.t1)
            
        #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        paths = [pS_7, pE_7, pN_7, pW_7, pStop, pEtop, pNtop, pWtop, midInArc, aftArc, midOutArc, foreArc]
        volume7 = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                  paths[4], paths[5], paths[6], paths[7],
                                  paths[8], paths[9], paths[10], paths[11])
        
        self.southFace7 = get_Face(paths, 'S') if self.face == 1 else make_patch(pStop, aftArc.reverse(), pS_7, midInArc.reverse()) #make_patch(pStop, aftArc, pS_7, midInArc)
        self.eastFace7 = get_Face(paths, 'E') if self.face == 1 else make_patch(pEtop, midOutArc.reverse(), pE_7, aftArc) #make_patch(pEtop, midOutArc, pE_7, aftArc)
        self.flatPaths7 =  [pN_7_flat, pE_7_flat, pS_7_flat, pW_7_flat]        
        return volume7

    def diamondTEMP(self):
        pNtop = self.pN.clone().rotate_about_zaxis(self.side*self.angle)
        pEtop = self.pE.clone().rotate_about_zaxis(self.side*self.angle)
        pStop = self.pS.clone().rotate_about_zaxis(self.side*self.angle)
        pWtop = self.pW.clone().rotate_about_zaxis(self.side*self.angle)
        
        foreInArc = Arc(foreIn, pNtop.eval(0.), Node(0.,0.,foreIn.z))
        foreOutArc = Arc(foreOut, pNtop.eval(1.), Node(0.,0.,foreOut.z))
        aftInArc = Arc(aftIn, pStop.eval(0.), Node(0.,0.,aftIn.z))
        aftOutArc = Arc(aftOut, pStop.eval(1.), Node(0.,0.,aftOut.z))
            
            
        #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        volumeTEMP = WireFrameVolume(pSt, pEt, pNt, pWt, 
                                 pStop, pEtop, pNtop, pWtop, 
                                 aftInArc, aftOutArc, foreOutArc, foreInArc)
        
        return volumeTEMP
    
    def rotate(self, angleStart, angleEnd, flatPaths):
        rotation = angleEnd - angleStart
        
        pNflat = flatPaths[0]; pEflat = flatPaths[1];
        pSflat = flatPaths[2]; pWflat = flatPaths[3]
        
        pN = pNflat.clone().rotate_about_zaxis(angleStart)
        pE = pEflat.clone().rotate_about_zaxis(angleStart)
        pS = pSflat.clone().rotate_about_zaxis(angleStart)
        pW = pWflat.clone().rotate_about_zaxis(angleStart)
        
        pNtop = pN.clone().rotate_about_zaxis(rotation)
        pEtop = pE.clone().rotate_about_zaxis(rotation)
        pStop = pS.clone().rotate_about_zaxis(rotation)
        pWtop = pW.clone().rotate_about_zaxis(rotation)

        foreArc = Arc(pN.eval(0), pNtop.eval(0.), Node(0.,0.,pN.eval(0).z))
        midOutArc = Arc(pN.eval(1), pNtop.eval(1.), Node(0.,0.,pN.eval(1).z))
        midInArc = Arc(pS.eval(0), pStop.eval(0.), Node(0.,0.,pS.eval(0).z))
        aftArc = Arc(pS.eval(1), pStop.eval(1.), Node(0.,0.,pS.eval(1).z))
        

        paths = [pS, pE, pN, pW, 
                 pStop, pEtop, pNtop, pWtop, 
                 midInArc, aftArc, midOutArc, foreArc]
        
        #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        rotatedVolume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                        paths[4], paths[5], paths[6], paths[7],
                                        paths[8], paths[9], paths[10], paths[11])  
        
        self.southFace7 = get_Face(paths, 'S')
        self.eastFace7 = get_Face(paths, 'E')
        self.northFace6 = get_Face(paths, 'N')
        self.eastFace6 = get_Face(paths, 'E')
        self.westFace5 = get_Face(paths, 'W')
        self.southFace5 = get_Face(paths, 'S')
        return rotatedVolume    

class rearBlocks:
    def __init__(self, startAngle, endAngle, side, face='TOP'):
        """ Function to create the blocks along the leading edge. Firstly, 
        the surfaces are created, surf2surf is then used to return the coordinates.
        Inputs:
            r,s,t   		  - parametric positons within the volume
            r_leadingEdge  	  - radius of leading edge [m]
            theta_leadingEdge - angle of leading edge blocks
            halfspan	      - halfspan of wings (from centreline) [m]
            sweep			  - sweep angle of wings (from centreline) [rads]
            zStart			  - z position of wing leading edge at root [m]
            
        Returns
        Coordinates defining a block.
        """
        if side == 'LHS':
            self.side = 1.
        elif side == 'RHS':
            self.side = -1.
        else:
            print 'Incorrect wing side choice'; return None
        if face == 'TOP':
            self.face = 1.
        else:
            self.face -1
        self.endAngle = endAngle
        self.startAngle = startAngle
    
    def block(self, length, width, zStart, xStart, block, rotation=0, TEflag=0):
        dxFore = abs(zStartboatTail - zStart)/np.tan(phi_sweep)
        dxAft = abs(zStartboatTail - zStart + length)/np.tan(phi_sweep)
        self.TEflag = TEflag
        #print 'this', xFore, xAft, side

        if block =='EDGE':
            foreIn_flat = Node(self.side*xStart, 0., zStart, label='foreIn_flat')
            foreOut_flat = Node(self.side*(halfspan-dxFore), 0., foreIn_flat.z, label='foreOut_flat')
            aftIn_flat = Node(foreIn_flat.x, 0., foreIn_flat.z - length, label='aftOut_flat')
            aftOut_flat = Node(self.side*(halfspan-dxAft), 0., aftIn_flat.z, label='aftOut_flat')

            foreIn = Node(self.side*xStart, wingThickness(), zStart, label='foreIn')
            foreOut = Node(self.side*(halfspan-dxFore), wingThickness(), foreIn.z, label='foreOut')
            aftIn = Node(foreIn.x, wingThickness(), foreIn.z - length, label='aftOut')
            aftOut = Node(self.side*(halfspan-dxAft), wingThickness(), aftIn.z, label='aftOut')
        
        else:
            foreIn_flat = Node(self.side*xStart, 0., zStart, label='foreIn_flat')
            foreOut_flat = Node(foreIn_flat.x+self.side*width, 0., foreIn_flat.z, label='foreOut_flat')
            aftIn_flat = Node(foreIn_flat.x, 0., foreIn_flat.z - length, label='aftOut_flat')
            aftOut_flat = Node(foreOut_flat.x, 0., aftIn_flat.z, label='aftOut_flat')
            
            foreIn = Node(self.side*xStart, wingThickness(), zStart, label='foreIn')
            foreOut = Node(foreIn.x+self.side*width, wingThickness(), foreIn.z, label='foreOut')
            aftIn = Node(foreIn.x, wingThickness(), foreIn.z - length, label='aftOut')
            aftOut = Node(foreOut.x, wingThickness(), aftIn.z, label='aftOut')
        
        
        pN_flat = Line(foreIn_flat, foreOut_flat)
        self.pE_flat = Line(aftOut_flat, foreOut_flat)
        self.pS_flat = Line(aftIn_flat, aftOut_flat)
        pW_flat = Line(aftIn_flat, foreIn_flat)        
        
        pN = Line(foreIn, foreOut)
        self.pE = Line(aftOut, foreOut)
        self.pS = Line(aftIn, aftOut)
        pW = Line(aftIn, foreIn)
       
        pNtop = pN_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        self.pEtop = self.pE_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        self.pStop = self.pS_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pWtop = pW_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
    
    
        if block == 'EDGE':
            foreInArc = Line(foreIn, pNtop.eval(0.))
            foreOutArc = Line(foreOut, pNtop.eval(1.))
            aftInArc = Line(aftIn, self.pStop.eval(0.))
            aftOutArc = Line(aftOut, self.pStop.eval(1.))

        else:
            foreInArc = Line(foreIn, pNtop.eval(0.))
            foreOutArc = Line(foreOut, pNtop.eval(1.))
            aftInArc = Line(aftIn, self.pStop.eval(0.))
            aftOutArc = Line(aftOut, self.pStop.eval(1.))
            
            
        paths = [self.pS, self.pE, pN, pW, 
                 self.pStop, self.pEtop, pNtop, pWtop, 
                 aftInArc, aftOutArc, foreOutArc, foreInArc]        
        self.volume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                      paths[4], paths[5], paths[6], paths[7],
                                      paths[8], paths[9], paths[10], paths[11])
                                     
        #return volume, pE, pEtop, pS, pStop
        self.surfacePath = self.pE
        self.surfacePath_outer = self.pEtop
        self.flatPaths =  [pN_flat, self.pE_flat, self.pS_flat, pW_flat]                
        self.northFace = get_Face(paths, 'N')
        self.rOuterStart = foreOut_flat.x*self.side
        self.rOuterEnd = aftOut_flat.x*self.side
        return self.volume

    def undersideBlock(self, length, width, zStart, xStart, block, wallAngle=0, rotation=0, TEflag=0):
        print length, width, zStart, xStart, block, wallAngle, rotation
        xFore = abs(zStartboatTail - zStart)/np.tan(phi_sweep)
        xAft = abs(zStartboatTail - zStart + length)/np.tan(phi_sweep)
        self.TEflag = TEflag
        #print 'this', xFore, xAft, side

        if block == 'EDGE':
            foreIn_flat = Node(self.side*xStart, 0., zStart, label='foreIn_flat')
            foreOut_flat = Node(self.side*(halfspan-xFore), 0., foreIn_flat.z, label='foreOut_flat')
            aftIn_flat = Node(foreIn_flat.x, 0., foreIn_flat.z - length, label='aftOut_flat')
            aftOut_flat = Node(self.side*(halfspan-xAft), 0., aftIn_flat.z, label='aftOut_flat')

            foreIn = Node(self.side*xStart, -1*wingThickness(), zStart, label='foreIn')
            foreOut = Node(self.side*(halfspan-xFore), -1*wingThickness(), foreIn.z, label='foreOut')
            aftIn = Node(foreIn.x, -1*wingThickness(), foreIn.z - length, label='aftOut')
            aftOut = Node(self.side*(halfspan-xAft), -1*wingThickness(), aftIn.z, label='aftOut')
        
        elif block =='WALL':
            foreIn_flat = Node(self.side*xStart, 0., zStart, label='foreIn_flat')
            foreOut_flat = Node(foreIn_flat.x+self.side*width, 0., foreIn_flat.z, label='foreOut_flat')
            aftIn_flat = Node(foreIn_flat.x+self.side*length*np.tan(wallAngle), 0., foreIn_flat.z - length, label='aftOut_flat')
            aftOut_flat = Node(foreOut_flat.x, 0., aftIn_flat.z, label='aftOut_flat')
            
            foreIn = Node(self.side*xStart, -1*wingThickness(), zStart, label='foreIn')
            foreOut = Node(foreIn.x+self.side*width, -1*wingThickness(), foreIn.z, label='foreOut')
            aftIn = Node(foreIn.x+self.side*length*np.tan(wallAngle), -1*wingThickness(), foreIn.z - length, label='aftOut')
            aftOut = Node(foreOut.x, -1*wingThickness(), aftIn.z, label='aftOut')                     
        
        else:
            foreIn_flat = Node(self.side*xStart, 0., zStart, label='foreIn_flat')
            foreOut_flat = Node(foreIn_flat.x+self.side*width, 0., foreIn_flat.z, label='foreOut_flat')
            aftIn_flat = Node(foreIn_flat.x, 0., foreIn_flat.z - length, label='aftOut_flat')
            aftOut_flat = Node(foreOut_flat.x, 0., aftIn_flat.z, label='aftOut_flat')
            
            foreIn = Node(self.side*xStart, -1*wingThickness(), zStart, label='foreIn')
            foreOut = Node(foreIn.x+self.side*width, -1*wingThickness(), foreIn.z, label='foreOut')
            aftIn = Node(foreIn.x, -1*wingThickness(), foreIn.z - length, label='aftOut')
            aftOut = Node(foreOut.x, -1*wingThickness(), aftIn.z, label='aftOut')
        
        pN_flat = Line(foreIn_flat, foreOut_flat)
        self.pE_flat = Line(aftOut_flat, foreOut_flat)
        self.pS_flat = Line(aftIn_flat, aftOut_flat)
        pW_flat = Line(aftIn_flat, foreIn_flat)        
        
        pN = Line(foreIn, foreOut)
        self.pE = Line(aftOut, foreOut)
        self.pS = Line(aftIn, aftOut)
        pW = Line(aftIn, foreIn)
       
        pNtop = pN_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        self.pEtop = self.pE_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        self.pStop = self.pS_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pWtop = pW_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
    
    
        if block == 'EDGE':
            foreInArc = Line(pNtop.eval(0.), foreIn)
            foreOutArc = Line(pNtop.eval(1.), foreOut)
            aftInArc = Line(self.pStop.eval(0.), aftIn)
            aftOutArc = Line(self.pStop.eval(1.), aftOut)

        else:
            foreInArc = Line(pNtop.eval(0.), foreIn)
            foreOutArc = Line(pNtop.eval(1.), foreOut)
            aftInArc = Line(self.pStop.eval(0.), aftIn)
            aftOutArc = Line(self.pStop.eval(1.), aftOut)
            
        paths = [self.pS, self.pE, pN, pW, 
                 self.pStop, self.pEtop, pNtop, pWtop, 
                 aftInArc, aftOutArc, foreOutArc, foreInArc]        
        self.volume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                      paths[4], paths[5], paths[6], paths[7],
                                      paths[8], paths[9], paths[10], paths[11])
                                     
        #return volume, pE, pEtop, pS, pStop
        self.surfacePath = self.pE
        self.surfacePath_outer = self.pEtop
        self.flatPaths =  [pN_flat, self.pE_flat, self.pS_flat, pW_flat]                
        self.northFace = make_patch(pNtop, foreOutArc.reverse(), pN, foreInArc.reverse())
        self.rOuterStart = foreOut_flat.x*self.side
        self.rOuterEnd = aftOut_flat.x*self.side
        return self.volume

    def rotate(self, angleStart, angleEnd, flatPaths):
        rotation = angleEnd - angleStart
        
        pNflat = flatPaths[0]; pEflat = flatPaths[1];
        pSflat = flatPaths[2]; pWflat = flatPaths[3]
        
        pN = pNflat.clone().rotate_about_zaxis(angleStart)
        pE = pEflat.clone().rotate_about_zaxis(angleStart)
        pS = pSflat.clone().rotate_about_zaxis(angleStart)
        pW = pWflat.clone().rotate_about_zaxis(angleStart)
        
        pNtop = pN.clone().rotate_about_zaxis(rotation)
        pEtop = pE.clone().rotate_about_zaxis(rotation)
        pStop = pS.clone().rotate_about_zaxis(rotation)
        pWtop = pW.clone().rotate_about_zaxis(rotation)

        foreInArc = Arc(pN.eval(0), pNtop.eval(0.), Node(0.,0.,pN.eval(0).z))
        foreOutArc = Arc(pN.eval(1), pNtop.eval(1.), Node(0.,0.,pN.eval(1).z))
        aftInArc = Arc(pS.eval(0), pStop.eval(0.), Node(0.,0.,pS.eval(0).z)) if self.TEflag==0 and self.face==1 else Line(pS.eval(0), pStop.eval(0.))
        aftOutArc = Arc(pS.eval(1), pStop.eval(1.), Node(0.,0.,pS.eval(1).z))
        

        paths = [pS, pE, pN, pW, 
                 pStop, pEtop, pNtop, pWtop, 
                 aftInArc, aftOutArc, foreOutArc, foreInArc]
        #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        rotatedVolume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                        paths[4], paths[5], paths[6], paths[7],
                                        paths[8], paths[9], paths[10], paths[11])  
        self.surfacePath = self.pE
        self.surfacePath_outer = self.pEtop
        self.northFace = get_Face(paths, 'N')
        return rotatedVolume    

class outerBlocks:
    def __init__(self, endAngle, zStart, zEnd, rInnerStart, rOuterStart, rInnerEnd, side, face='TOP'):
        """ Function to create the blocks along the leading edge. Firstly, 
        the surfaces are created, surf2surf is then used to return the coordinates.
        Inputs:
            r,s,t   		  - parametric positons within the volume
            r_leadingEdge  	  - radius of leading edge [m]
            theta_leadingEdge - angle of leading edge blocks
            halfspan	      - halfspan of wings (from centreline) [m]
            sweep			  - sweep angle of wings (from centreline) [rads]
            zStart			  - z position of wing leading edge at root [m]
            
        Returns
        Coordinates defining a block.
        """
        if side == 'LHS':
            self.side = 1.
        elif side == 'RHS':
            self.side = -1.
        else:
            print 'Incorrect wing side choice'; return None
        if face == 'TOP':
            self.face = 1.
        else:
            self.face -1
        self.endAngle = endAngle
        
        
        
        self.zStart = zStart
        self.zEnd = zEnd
        self.rInnerStart = rInnerStart
        self.rOuterStart = rOuterStart
        self.rInnerEnd = rInnerEnd
        self.rOuterEnd = rOuterStart + abs(self.zEnd - self.zStart)*np.tan((90.*np.pi/180. - phi_sweep) + phi_meshWing)
        
        
    def undersideBlock(self):
        foreIn = Node(self.side*(self.rInnerStart), -1.*wingThickness(), self.zStart)
        foreOut = Node(self.side*(self.rOuterStart), -1.*wingThickness(), self.zStart)
        aftIn = Node(self.side*(self.rInnerEnd), -1.*wingThickness(), self.zEnd)
        aftOut = Node(self.side*(self.rOuterEnd), -1.*wingThickness(), self.zEnd)

        foreIn_flat = Node(self.side*self.rInnerStart, 0., self.zStart, label='foreIn_flat')
        foreOut_flat = Node(self.side*(self.rOuterStart), 0., self.zStart, label='foreOut_flat')
        aftIn_flat = Node(self.side*self.rInnerEnd, 0., self.zEnd, label='aftOut_flat')
        aftOut_flat = Node(self.side*(self.rOuterEnd), 0., self.zEnd, label='aftOut_flat')

        pN_flat = Line(foreIn_flat, foreOut_flat)
        self.pE_flat = Line(aftOut_flat, foreOut_flat)
        self.pS_flat = Line(aftIn_flat, aftOut_flat)
        pW_flat = Line(aftIn_flat, foreIn_flat)        
        
        pN = Line(foreIn, foreOut)
        self.pE = Line(aftOut, foreOut)
        self.pS = Line(aftIn, aftOut)
        pW = Line(aftIn, foreIn)
       
        pNtop = pN_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        self.pEtop = self.pE_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        self.pStop = self.pS_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pWtop = pW_flat.clone().rotate_about_zaxis(self.side*self.endAngle)

        foreInArc = Line(pNtop.eval(0.), foreIn)
        foreOutArc = Line(pNtop.eval(1.), foreOut)
        aftInArc = Line(self.pStop.eval(0.), aftIn)
        aftOutArc = Line(self.pStop.eval(1.), aftOut)
            
        paths = [self.pS, self.pE, pN, pW, 
                 self.pStop, self.pEtop, pNtop, pWtop, 
                 aftInArc, aftOutArc, foreOutArc, foreInArc]        
        self.volume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                      paths[4], paths[5], paths[6], paths[7],
                                      paths[8], paths[9], paths[10], paths[11])
                                     
        self.flatPaths =  [pN_flat, self.pE_flat, self.pS_flat, pW_flat]                
        return self.volume
        
    def topsideBlock(self):
        foreIn = Node(self.side*(self.rInnerStart), wingThickness(), self.zStart)
        foreOut = Node(self.side*(self.rOuterStart), wingThickness(), self.zStart)
        aftIn = Node(self.side*(self.rInnerEnd), wingThickness(), self.zEnd)
        aftOut = Node(self.side*(self.rOuterEnd), wingThickness(), self.zEnd)

        foreIn_flat = Node(self.side*self.rInnerStart, 0., self.zStart, label='foreIn_flat')
        foreOut_flat = Node(self.side*(self.rOuterStart), 0., self.zStart, label='foreOut_flat')
        aftIn_flat = Node(self.side*self.rInnerEnd, 0., self.zEnd, label='aftOut_flat')
        aftOut_flat = Node(self.side*(self.rOuterEnd), 0., self.zEnd, label='aftOut_flat')

        pN_flat = Line(foreIn_flat, foreOut_flat)
        self.pE_flat = Line(aftOut_flat, foreOut_flat)
        self.pS_flat = Line(aftIn_flat, aftOut_flat)
        pW_flat = Line(aftIn_flat, foreIn_flat)        
        
        pN = Line(foreIn, foreOut)
        self.pE = Line(aftOut, foreOut)
        self.pS = Line(aftIn, aftOut)
        pW = Line(aftIn, foreIn)
       
        pNtop = pN_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        self.pEtop = self.pE_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        self.pStop = self.pS_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pWtop = pW_flat.clone().rotate_about_zaxis(self.side*self.endAngle)

        foreInArc = Line(foreIn, pNtop.eval(0.))
        foreOutArc = Line(foreOut, pNtop.eval(1.))
        aftInArc = Line(aftIn, self.pStop.eval(0.))
        aftOutArc = Line(aftOut, self.pStop.eval(1.))
            
        paths = [self.pS, self.pE, pN, pW, 
                 self.pStop, self.pEtop, pNtop, pWtop, 
                 aftInArc, aftOutArc, foreOutArc, foreInArc]        
        self.volume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                      paths[4], paths[5], paths[6], paths[7],
                                      paths[8], paths[9], paths[10], paths[11])
                                     
        self.flatPaths =  [pN_flat, self.pE_flat, self.pS_flat, pW_flat]                
        return self.volume        
        
    def rotate(self, angleStart, angleEnd, flatPaths):
        rotation = angleEnd - angleStart
        
        pNflat = flatPaths[0]; pEflat = flatPaths[1];
        pSflat = flatPaths[2]; pWflat = flatPaths[3]
        
        pN = pNflat.clone().rotate_about_zaxis(angleStart)
        pE = pEflat.clone().rotate_about_zaxis(angleStart)
        pS = pSflat.clone().rotate_about_zaxis(angleStart)
        pW = pWflat.clone().rotate_about_zaxis(angleStart)
        
        pNtop = pN.clone().rotate_about_zaxis(rotation)
        pEtop = pE.clone().rotate_about_zaxis(rotation)
        pStop = pS.clone().rotate_about_zaxis(rotation)
        pWtop = pW.clone().rotate_about_zaxis(rotation)

        foreInArc = Arc(pN.eval(0), pNtop.eval(0.), Node(0, 0, pN.eval(0).z))
        foreOutArc = Arc(pN.eval(1), pNtop.eval(1.), Node(0.,0.,pN.eval(1).z))
        aftInArc = Arc(pS.eval(0), pStop.eval(0.), Node(0, 0, pS.eval(0).z))
        aftOutArc = Arc(pS.eval(1), pStop.eval(1.), Node(0.,0.,pS.eval(1).z))
        

        paths = [pS, pE, pN, pW, 
                 pStop, pEtop, pNtop, pWtop, 
                 aftInArc, aftOutArc, foreOutArc, foreInArc]
        #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        rotatedVolume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                        paths[4], paths[5], paths[6], paths[7],
                                        paths[8], paths[9], paths[10], paths[11])  
        #self.surfacePath = self.pE
        #self.surfacePath_outer = self.pEtop
        #self.northFace = get_Face(paths, 'N')
        return rotatedVolume 

#---------------------------------- Wing -------------------------------------#
def wingThickness(xCoord=0, zCoord=0, side=0, face='TOP'):
    """ This function returns the local wing thickness given the (x, y) 
    coordinate.
    
    """

#    tMaxLine = np.arctan((halfspan-r_fb)/(c_t))
#    leadAngleRoot = np.arctan((r_leading*np.sin(theta_leadingEdge) - t_root)/c_t)
#    trailAngleRoot = np.arctan((r_leading*np.sin(theta_leadingEdge) - t_root)/c_t)
#    leadAngle = lambda x : (x - r_fb)/(halfspan - r_fb)*leadAngleRoot
#    trailAngle = lambda x :(x - r_fb)/(halfspan - r_fb)*trailAngleRoot 
#
#    t_root

    if face == 'TOP':
        return r_leadingEdge
    else:
        return -1.*r_leadingEdge
        
class edgeRounding:
    """ Function to create the blocks along the leading edge. Firstly, 
    the surfaces are created, surf2surf is then used to return the coordinates.
    Inputs:
        r,s,t   		  - parametric positons within the volume
        r_leadingEdge  	  - radius of leading edge [m]
        theta_leadingEdge - angle of leading edge blocks(angle of airfoil)
        r_root      - radius of root LE rounding (in x-z) [m]
        chordRoot       - Chord length of wings at root [m]
        halfspan	      - halfspan of wings (from centreline) [m]
        sweep			  - sweep angle of wings (from centreline) [rads]
        zStart			  - z position of wing leading edge at root [m]
        side        - choose LHS or RHS to build wing?        
        t_leadingEdge - Thickness of LE blocks [m]
    
    NOTE: 
    """
    def __init__(self, side, face='TOP'):
        if side == 'LHS':
                self.side = 1.
        elif side == 'RHS':
            self.side = -1.
        else:
            print 'Incorrect wing side choice'; return None
        
        
    def block_LEfairing(self, r_leadingEdge, surfacePath, surfacePath_outer, angleWing, theta_leadingEdge=90.*np.pi/180, face='TOP'):
    
        tBlock = abs(surfacePath.eval(0).y - surfacePath_outer.eval(0).y)
        rOuter = r_leadingEdge + tBlock
        
        angleBottomFace = wallAngle_inlet
    
        # Meets Forebody here
        centreStart = Node(surfacePath.eval(0).x - self.side*r_leadingEdge*np.sin(phi_fb), 0., surfacePath.eval(0).z + r_leadingEdge*np.cos(phi_fb))
        centreStart_outer = Node(surfacePath.eval(0).x - self.side*rOuter*np.sin(phi_fb), 0., surfacePath.eval(0).z + rOuter*np.cos(phi_fb))
        
        # Rest of the wing down this end
        self.centreEnd = Node(surfacePath.eval(1).x + self.side*r_leadingEdge*np.cos(angleBottomFace), 0., surfacePath.eval(1).z + r_leadingEdge*np.sin(angleBottomFace))
        self.centreEnd_outer = Node(surfacePath.eval(1).x + self.side*rOuter*np.cos(angleBottomFace), 0., surfacePath.eval(1).z + rOuter*np.sin(angleBottomFace))
        
        centrePath = Line(centreStart, self.centreEnd)
        centrePath_outer = Line(centreStart_outer, self.centreEnd_outer)
    
        pN = Line(centrePath_outer.eval(0), surfacePath_outer.eval(0))
        pE = Line(surfacePath.eval(0), surfacePath_outer.eval(0))
        pS = Line(centreStart, surfacePath.eval(0))
        pW = Line(centreStart, centrePath_outer.eval(0))
        #print 'pN', pN, 'pE', pE, 'pS', pS, 'pW', pW    
    
        pNtop = Line(centrePath_outer.eval(1), surfacePath_outer.eval(1))
        pEtop = Line(surfacePath.eval(1), surfacePath_outer.eval(1))
        pStop = Line(self.centreEnd, surfacePath.eval(1))
        pWtop = Line(self.centreEnd, centrePath_outer.eval(1))
        #print 'pN', pN, 'pE', pE, 'pS', pS, 'pW', pW    
    
        # p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        # bottom, top, between
        self.volume = WireFrameVolume(pS, pE, pN, pW, 
                                 pStop, pEtop, pNtop, pWtop, 
                                 centrePath_outer, surfacePath_outer, surfacePath, centrePath)
        return self.volume
        
    def block_LE(self, r_leadingEdge, surfacePath, surfacePath_outer, centreStart,
                 centreStart_outer, angleWing, theta_leadingEdge=90.*np.pi/180, endFlag=0, TEFlag=0, face='TOP'):
                     
        tBlock = abs(surfacePath.eval(0).y - surfacePath_outer.eval(0).y)
        rOuter = r_leadingEdge + tBlock
    
        surfacePathEnd = surfacePath.eval(1); surfacePathStart = surfacePath.eval(0)
    
        if endFlag == 0 and TEFlag == 0: 
            self.centreEnd = Node(surfacePathStart.x + self.side*r_leadingEdge*np.sin(phi_sweep), 0., surfacePathStart.z + r_leadingEdge*np.cos(phi_sweep))
            self.centreEnd_outer = Node(surfacePathStart.x + self.side*rOuter*np.sin(phi_sweep), 0., surfacePathStart.z + rOuter*np.cos(phi_sweep))
        elif TEFlag == 1:
            self.centreEnd = Node(surfacePathStart.x, 0., surfacePathStart.z - r_leadingEdge)
            self.centreEnd_outer = Node(surfacePathStart.x, 0., surfacePathStart.z - rOuter)
        else:
            angleMitre = (360.*np.pi/180. - phi_sweep)/2. - 90.*np.pi/180.
            self.centreEnd = Node(surfacePathStart.x + self.side*r_leadingEdge*np.tan(angleMitre), 0., surfacePathStart.z - r_leadingEdge)
            self.centreEnd_outer = Node(surfacePathStart.x + self.side*rOuter*np.tan(angleMitre), 0., surfacePathStart.z - rOuter)
        
        centrePath = Line(self.centreEnd, centreStart)
        centrePath_outer = Line(self.centreEnd_outer, centreStart_outer)
    
        pN = Line(surfacePath_outer.eval(1), centrePath_outer.eval(1))
        pE = Line(self.centreEnd, centrePath_outer.eval(1))
        #pS = Arc(surfacePathEnd, centreEnd, Node(surfacePathEnd.x, 0., centreEnd.z))
        pS = Line(surfacePathEnd, self.centreEnd)
        pW = Line(surfacePathEnd, surfacePath_outer.eval(1))
        #print 'pN', pN, 'pE', pE, 'pS', pS, 'pW', pW    
    
        pNtop = Line(surfacePath_outer.eval(0), centrePath_outer.eval(0))
        pEtop = Line(centreStart, centrePath_outer.eval(0))
        #pStop = Arc(surfacePathStart, centreStart, Node(surfacePathStart.x, 0., centreStart.z))
        pStop = Line(surfacePathStart, centreStart)
        pStop = Line(surfacePathStart, centreStart)
        pWtop = Line(surfacePathStart, surfacePath_outer.eval(0))
        #print 'pN', pN, 'pE', pE, 'pS', pS, 'pW', pW    
    
        # p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        self.volume = WireFrameVolume(pS, pE, pN, pW, 
                                      pStop, pEtop, pNtop, pWtop, 
                                      centrePath, surfacePath, surfacePath_outer, centrePath_outer)   
        
        
        surfacePathEnd = Node(surfacePath.eval(1).x, -1*surfacePath.eval(1).y, surfacePath.eval(1).z); 
        surfacePathStart = Node(surfacePath.eval(0).x, -1*surfacePath.eval(0).y, surfacePath.eval(0).z);
        
        surfacePath_outerEnd = Node(surfacePath_outer.eval(1).x, -1*surfacePath_outer.eval(1).y, surfacePath_outer.eval(1).z); 
        surfacePath_outerStart = Node(surfacePath_outer.eval(0).x, -1*surfacePath_outer.eval(0).y, surfacePath_outer.eval(0).z);
        
        pNb = Line(surfacePath_outerEnd, centrePath_outer.eval(1))
        pEb = Line(self.centreEnd, centrePath_outer.eval(1))
        pSb = Line(surfacePathEnd, self.centreEnd)
        pWb = Line(surfacePathEnd, surfacePath_outerEnd)
        #print 'pN', pN, 'pE', pE, 'pS', pS, 'pW', pW    
    
        pNtopb = Line(surfacePath_outerStart, centrePath_outer.eval(0))
        pEtopb = Line(centreStart, centrePath_outer.eval(0))
        pStopb = Line(surfacePathStart, centreStart)
        pWtopb = Line(surfacePathStart, surfacePath_outerStart)
        
        self.underVolume = WireFrameVolume(pSb, pEb, pNb, pWb, 
                                           pStopb, pEtopb, pNtopb, pWtopb, 
                                           centrePath, surfacePath, surfacePath_outer, centrePath_outer)   
        return self.volume

def wing(flatPaths, block=''):
    flatSurf = make_patch(flatPaths[0], flatPaths[1], flatPaths[2], flatPaths[3])

    topPaths = [flatPaths[0].clone().translate(0., wingThickness(), 0.), 
                flatPaths[1].clone().translate(0., wingThickness(), 0.),
                flatPaths[2].clone().translate(0., wingThickness(), 0.),
                flatPaths[3].clone().translate(0., wingThickness(), 0.) ]
    
    bottomPaths = [flatPaths[0].clone().translate(0., -1*wingThickness(), 0.),
                   flatPaths[1].clone().translate(0., -1*wingThickness(), 0.),
                   flatPaths[2].clone().translate(0., -1*wingThickness(), 0.),
                   flatPaths[3].clone().translate(0., -1*wingThickness(), 0.)]
    
    foreInArc = Line(bottomPaths[0].eval(0), topPaths[0].eval(0.))
    foreOutArc = Line(bottomPaths[0].eval(1), topPaths[0].eval(1.))
    aftInArc = Line(bottomPaths[2].eval(0), topPaths[2].eval(0.))
    aftOutArc = Line(bottomPaths[2].eval(1), topPaths[2].eval(1.))

    paths = [bottomPaths[0] , bottomPaths[1], bottomPaths[2], bottomPaths[3],
             topPaths[0], topPaths[1], topPaths[2], topPaths[3], 
             aftInArc, aftOutArc, foreOutArc, foreInArc]         
    
    volume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                  paths[4], paths[5], paths[6], paths[7],
                                  paths[8], paths[9], paths[10], paths[11])
    return volume
        
#------------------------- Boat Tail and Wake --------------------------------#        
class BTandWake:
    def __init__(self, startAngle, endAngle, side, face='TOP'):
        """
        """
        if side == 'LHS':
            self.side = 1.
        elif side == 'RHS':
            self.side = -1.
        else:
            print 'Incorrect wing side choice'; return None
        if face == 'TOP':
            self.face = 1.
        else:
            self.face = -1
        self.endAngle = endAngle
        self.startAngle = startAngle
    
    def block(self, length, width, zStart, xStart, block, wallAngle=0, meshAngle=0):
        self.meshAngle = meshAngle
        xFore = abs(zStartboatTail - zStart)*np.sin(meshAngle)
        xAft = length*np.sin(meshAngle)
        #print 'this', xFore, xAft, side

        if block =='EDGE':
            foreIn_flat = Node(self.side*xStart, 0., zStart, label='foreIn_flat')
            foreOut_flat = Node(self.side*(halfspan+xFore), 0., foreIn_flat.z, label='foreOut_flat')
            aftIn_flat = Node(foreIn_flat.x, 0., foreIn_flat.z - length, label='aftOut_flat')
            aftOut_flat = Node(self.side*(halfspan+xFore+xAft), 0., aftIn_flat.z, label='aftOut_flat')

            foreIn = Node(self.side*xStart, wingThickness(), zStart, label='foreIn')
            foreOut = Node(self.side*(halfspan+xFore), wingThickness(), foreIn.z, label='foreOut')
            aftIn = Node(foreIn.x, wingThickness(), foreIn.z - length, label='aftOut')
            aftOut = Node(self.side*(halfspan+xFore+xAft), wingThickness(), aftIn.z, label='aftOut')
        
        else:
            foreIn_flat = Node(self.side*xStart, 0., zStart, label='foreIn_flat')
            foreOut_flat = Node(foreIn_flat.x+self.side*width, 0., foreIn_flat.z, label='foreOut_flat')
            aftIn_flat = Node(foreIn_flat.x+self.side*length*np.tan(wallAngle), 0., foreIn_flat.z - length, label='aftOut_flat')
            aftOut_flat = Node(foreOut_flat.x, 0., aftIn_flat.z, label='aftOut_flat')
            
            foreIn = Node(self.side*xStart, wingThickness(), zStart, label='foreIn')
            foreOut = Node(foreIn.x+self.side*width, wingThickness(), foreIn.z, label='foreOut')
            aftIn = Node(foreIn.x+self.side*length*np.tan(wallAngle), wingThickness(), foreIn.z - length, label='aftOut')
            aftOut = Node(foreOut.x, wingThickness(), aftIn.z, label='aftOut')
        
        
        pN_flat = Line(foreIn_flat, foreOut_flat)
        self.pE_flat = Line(aftOut_flat, foreOut_flat)
        self.pS_flat = Line(aftIn_flat, aftOut_flat)
        pW_flat = Line(aftIn_flat, foreIn_flat)        
        
        pN = Line(foreIn, foreOut)
        self.pE = Line(aftOut, foreOut)
        self.pS = Line(aftIn, aftOut)
        pW = Line(aftIn, foreIn)
       
        pNtop = pN_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        self.pEtop = self.pE_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        self.pStop = self.pS_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pWtop = pW_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
    
    
        if block == 'EDGE':
            foreInArc = Line(foreIn, pNtop.eval(0.))
            foreOutArc = Line(foreOut, pNtop.eval(1.))
            aftInArc = Line(aftIn, self.pStop.eval(0.))
            aftOutArc = Line(aftOut, self.pStop.eval(1.))

        else:
            foreInArc = Line(foreIn, pNtop.eval(0.))
            foreOutArc = Line(foreOut, pNtop.eval(1.))
            aftInArc = Line(aftIn, self.pStop.eval(0.))
            aftOutArc = Line(aftOut, self.pStop.eval(1.))
            
            
        paths = [self.pS, self.pE, pN, pW, 
                 self.pStop, self.pEtop, pNtop, pWtop, 
                 aftInArc, aftOutArc, foreOutArc, foreInArc]        
        self.volume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                      paths[4], paths[5], paths[6], paths[7],
                                      paths[8], paths[9], paths[10], paths[11])
                                     
        self.surfacePath = self.pE
        self.surfacePath_outer = self.pEtop
        self.flatPaths =  [pN_flat, self.pE_flat, self.pS_flat, pW_flat]                
        self.northFace = get_Face(paths, 'N')
        self.rOuterStart = foreOut_flat.x*self.side
        self.rOuterEnd = aftOut_flat.x*self.side
        return self.volume

    def undersideBlock(self, length, width, zStart, xStart, block, wallAngle=0, meshAngle=0):
        self.meshAngle = meshAngle
        xFore = abs(zStartboatTail - zStart)*np.sin(meshAngle)
        xAft = length*np.sin(meshAngle)
        #print 'this', xFore, xAft, side

        if block == 'EDGE':
            foreIn_flat = Node(self.side*xStart, 0., zStart, label='foreIn_flat')
            foreOut_flat = Node(self.side*(halfspan+xFore), 0., foreIn_flat.z, label='foreOut_flat')
            aftIn_flat = Node(foreIn_flat.x, 0., foreIn_flat.z - length, label='aftOut_flat')
            aftOut_flat = Node(self.side*(halfspan+xFore+xAft), 0., aftIn_flat.z, label='aftOut_flat')

            foreIn = Node(self.side*xStart, -1*wingThickness(), zStart, label='foreIn')
            foreOut = Node(self.side*(halfspan+xFore), -1*wingThickness(), foreIn.z, label='foreOut')
            aftIn = Node(foreIn.x, -1*wingThickness(), foreIn.z - length, label='aftOut')
            aftOut = Node(self.side*(halfspan+xFore+xAft), -1*wingThickness(), aftIn.z, label='aftOut')
        
        elif block =='WALL':
            self.block = 'WALL'
            foreIn_flat = Node(self.side*xStart, 0., zStart, label='foreIn_flat')
            foreOut_flat = Node(foreIn_flat.x+self.side*width, 0., foreIn_flat.z, label='foreOut_flat')
            aftIn_flat = Node(foreIn_flat.x+self.side*length*np.tan(wallAngle), 0., foreIn_flat.z - length, label='aftOut_flat')
            aftOut_flat = Node(foreOut_flat.x, 0., aftIn_flat.z, label='aftOut_flat')
            
            foreIn = Node(self.side*xStart, -1*wingThickness(), zStart, label='foreIn')
            foreOut = Node(foreIn.x+self.side*width, -1*wingThickness(), foreIn.z, label='foreOut')
            aftIn = Node(foreIn.x+self.side*length*np.tan(wallAngle), -1*wingThickness(), foreIn.z - length, label='aftOut')
            aftOut = Node(foreOut.x, -1*wingThickness(), aftIn.z, label='aftOut')                     
        
        else:
            foreIn_flat = Node(self.side*xStart, 0., zStart, label='foreIn_flat')
            foreOut_flat = Node(foreIn_flat.x+self.side*width, 0., foreIn_flat.z, label='foreOut_flat')
            aftIn_flat = Node(foreIn_flat.x, 0., foreIn_flat.z - length, label='aftOut_flat')
            aftOut_flat = Node(foreOut_flat.x, 0., aftIn_flat.z, label='aftOut_flat')   
            
            foreIn = Node(self.side*xStart, -1*wingThickness(), zStart, label='foreIn')
            foreOut = Node(foreIn.x+self.side*width, -1*wingThickness(), foreIn.z, label='foreOut')
            aftIn = Node(foreIn.x, -1*wingThickness(), foreIn.z - length, label='aftOut')
            aftOut = Node(foreOut.x, -1*wingThickness(), aftIn.z, label='aftOut')
        
        pN_flat = Line(foreIn_flat, foreOut_flat)
        self.pE_flat = Line(aftOut_flat, foreOut_flat)
        self.pS_flat = Line(aftIn_flat, aftOut_flat)
        pW_flat = Line(aftIn_flat, foreIn_flat)        
        
        pN = Line(foreIn, foreOut)
        self.pE = Line(aftOut, foreOut)
        self.pS = Line(aftIn, aftOut)
        pW = Line(aftIn, foreIn)
       
        pNtop = pN_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        self.pEtop = self.pE_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        self.pStop = self.pS_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
        pWtop = pW_flat.clone().rotate_about_zaxis(self.side*self.endAngle)
    
    
        if block == 'EDGE':
            foreInArc = Line(pNtop.eval(0.), foreIn)
            foreOutArc = Line(pNtop.eval(1.), foreOut)
            aftInArc = Line(self.pStop.eval(0.), aftIn)
            aftOutArc = Line(self.pStop.eval(1.), aftOut)

        else:
            foreInArc = Line(pNtop.eval(0.), foreIn)
            foreOutArc = Line(pNtop.eval(1.), foreOut)
            aftInArc = Line(self.pStop.eval(0.), aftIn)
            aftOutArc = Line(self.pStop.eval(1.), aftOut)
            
        paths = [self.pS, self.pE, pN, pW, 
                 self.pStop, self.pEtop, pNtop, pWtop, 
                 aftInArc, aftOutArc, foreOutArc, foreInArc]        
        self.volume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                      paths[4], paths[5], paths[6], paths[7],
                                      paths[8], paths[9], paths[10], paths[11])
                                     
        #return volume, pE, pEtop, pS, pStop
        self.surfacePath = self.pE
        self.surfacePath_outer = self.pEtop
        self.flatPaths =  [pN_flat, self.pE_flat, self.pS_flat, pW_flat]                
        self.northFace = make_patch(pNtop, foreOutArc.reverse(), pN, foreInArc.reverse())
        self.rOuterStart = foreOut_flat.x*self.side
        self.rOuterEnd = aftOut_flat.x*self.side
        return self.volume

    def rotate(self, angleStart, angleEnd, flatPaths, block=''):
        self.block = block
        rotation = angleEnd - angleStart
        
        pNflat = flatPaths[0]; pEflat = flatPaths[1];
        pSflat = flatPaths[2]; pWflat = flatPaths[3]
        
        pN = pNflat.clone().rotate_about_zaxis(angleStart)
        pE = pEflat.clone().rotate_about_zaxis(angleStart)
        pS = pSflat.clone().rotate_about_zaxis(angleStart)
        pW = pWflat.clone().rotate_about_zaxis(angleStart)
        
        pNtop = pN.clone().rotate_about_zaxis(rotation)
        pEtop = pE.clone().rotate_about_zaxis(rotation)
        pStop = pS.clone().rotate_about_zaxis(rotation)
        pWtop = pW.clone().rotate_about_zaxis(rotation)

        foreInArc = Line(pN.eval(0), pNtop.eval(0.)) if self.block == 'WALL' and self.face==-1 else Arc(pN.eval(0), pNtop.eval(0.), Node(0, 0, pN.eval(0).z))
        foreOutArc = Arc(pN.eval(1), pNtop.eval(1.), Node(0.,0.,pN.eval(1).z))
        aftInArc = Line(pS.eval(0), pStop.eval(0.))  if self.block == 'WALL' and self.face==-1 else Arc(pS.eval(0), pStop.eval(0.), Node(0, 0, pS.eval(0).z))
        aftOutArc = Arc(pS.eval(1), pStop.eval(1.), Node(0.,0.,pS.eval(1).z))
        

        paths = [pS, pE, pN, pW, 
                 pStop, pEtop, pNtop, pWtop, 
                 aftInArc, aftOutArc, foreOutArc, foreInArc]
        #p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37
        rotatedVolume = WireFrameVolume(paths[0], paths[1], paths[2], paths[3],
                                        paths[4], paths[5], paths[6], paths[7],
                                        paths[8], paths[9], paths[10], paths[11])  
        self.surfacePath = self.pE
        self.surfacePath_outer = self.pEtop
        self.northFace = get_Face(paths, 'N')
        return rotatedVolume    

    def centreBlocks(self):
        dx = r_bt/3.
        dy = wingThickness()
        centre = Node(0., 0., zStartWake)
        
        # This function quickly and easily creates the outer nodes.
        emsTop = lambda angle: Node(self.side*r_bt*np.cos(angle), r_bt*np.sin(angle), zStartWake)
        emsBottom = lambda angle: Node(self.side*r_bt*np.cos(angle), r_bt*np.sin(angle), zStartWake)
        
        M0 = emsTop(angle_90)
        M1 = emsTop(angle_Tail)
        M2 = emsTop(angle_WingTop)
        M3 = Node(self.side*r_bt, wingThickness(), zStartWake)    
        M4 = Node(self.side*r_bt, -1*wingThickness(), zStartWake)    
        M5 = emsBottom(-1*angle_WingBottom)
        M6 = emsBottom(engineSmile2[1])
        M7 = emsBottom(engineSmile2[0])
        #print engineSmile2
        #print 'M0', M0, 'M1', M1, 'M2', M2, 'M3', M3, 'M4', M4, 'M5', M5, 'M6', M6
        
        N0 = Node(self.side*0., 3.25*dy, zStartWake)
        N1 = Node(self.side*dx, 3.25*dy, zStartWake)
        N2 = Node(self.side*2.*dx, dy*1.75, zStartWake)
        N3 = Node(self.side*2.*dx, dy, zStartWake)
        
        N4 = Node(self.side*2.*dx, -1.*dy, zStartWake)
        N5 = Node(self.side*2.*dx, -1.75*dy, zStartWake)
        N6 = Node(self.side*dx, -3.25*dy, zStartWake)
        N7 = Node(self.side*0., -3.25*dy, zStartWake)
        
        P0 = Node(0., dy, zStartWake)
        P1 = Node(self.side*dx, dy, zStartWake)
        P2 = Node(self.side*dx, -1.*dy, zStartWake)
        P3 = Node(0., -1.*dy, zStartWake)
        
        
        # M circumferential paths
        M0_M1 = Arc(M0, M1, centre)
        M1_M2 = Arc(M1, M2, centre)
        M3_M2 = Line(M3, M2)
        M4_M3 = Line(M4, M3)
        M5_M4 = Line(M5, M4)
        M6_M5 = Line(M6, M5)
        M7_M6 = Line(M7, M6)
        
        
        # N circumferential paths
        N0_N1 = Line(N0, N1)
        N1_N2 = Line(N1, N2)  
        N3_N2 = Line(N3, N2)  
        N4_N3 = Line(N4, N3)  
        N5_N4 = Line(N5, N4)  
        N6_N5 = Line(N6, N5) 
        N7_N6 = Line(N7, N6)
        
        # P circumferential paths
        P0_P1 = Line(P0, P1)
        P2_P1 = Line(P2, P1)
        P3_P2 = Line(P3, P2)
        P3_P0 = Line(P3, P0)
        P1_N3 = Line(P1, N3)
        P2_N4 = Line(P2, N4)
        # Inner radial paths
        P0_N0 = Line(P0, N0)
        P1_N1 = Line(P1, N1)
        N6_P2 = Line(N6, P2)
        N7_P3 = Line(N7, P3)
        
        # Outer radial paths        
        N0_M0 = Line(N0, M0)
        N1_M1 = Line(N1, M1)
        N2_M2 = Line(N2, M2)
        N3_M3 = Line(N3, M3)
        N4_M4 = Line(N4, M4)
        N5_M5 = Line(N5, M5)
        M6_N6 = Line(M6, N6); N6_M6 = Line(N6, M6)
        M7_N7 = Line(M7, N7)
        
        # This function helps extrude the centre wake blocks the req. length.
        paths = lambda node: Line(node, Node(node.x, node.y, zStartWake-L_aftFarField))
        
        # pN, pE, pS, pW, grid_type='TFI' OR 'AO'
        surf0 = make_patch(M0_M1, N1_M1, N0_N1, N0_M0)
        surf1 = make_patch(M1_M2, N2_M2, N1_N2, N1_M1)
        surf2 = make_patch(N2_M2, M3_M2, N3_M3, N3_N2)
        surf3 = make_patch(N3_M3, M4_M3, N4_M4, N4_N3)
        surf4 = make_patch(N4_M4, M5_M4, N5_M5, N5_N4)
        surf5 = make_patch(N5_M5, M6_M5, N6_M6, N6_N5)
        surf6 = make_patch(N7_N6, M6_N6, M7_M6, M7_N7) 
        surf7 = make_patch(P3_P2, N6_P2, N7_N6, N7_P3)
        surf8 = make_patch(P2_N4, N5_N4, N6_N5, N6_P2)
        surf9 = make_patch(P0_P1, P2_P1, P3_P2, P3_P0)
        surf10 = make_patch(P1_N3, N4_N3, P2_N4, P2_P1)
        surf11 = make_patch(N1_N2, N3_N2, P1_N3, P1_N1)
        surf12 = make_patch(N0_N1, P1_N1, P0_P1, P0_N0)

#        surf0 = make_patch(M0_M1, N1_M1, N0_N1, N0_M0, grid_type='AO')
#        surf1 = make_patch(M1_M2, N2_M2, N1_N2, N1_M1, grid_type='AO')
#        surf2 = make_patch(N2_M2, M3_M2, N3_M3, N3_N2, grid_type='AO')
#        surf3 = make_patch(N3_M3, M4_M3, N4_M4, N4_N3, grid_type='AO')
#        surf4 = make_patch(N4_M4, M5_M4, N5_M5, N5_N4, grid_type='AO')
#        surf5 = make_patch(N6_N5, M5_N5, M6_M5, M6_N6, grid_type='AO')
#        surf6 = make_patch(P0_P1, N5_P1, N6_N5, N6_P0, grid_type='AO') 
#        surf7 = make_patch(P1_N3, N4_N3, N5_N4, N5_P1, grid_type='AO')
#        surf8 = make_patch(N1_N2, N3_N2, P1_N3, P1_N1, grid_type='AO')
#        surf9 = make_patch(N0_N1, P1_N1, P0_P1, P0_N0, grid_type='AO')

        self.centre0 = WireFrameVolume(surf0, paths(N0))
        self.centre1 = WireFrameVolume(surf1, paths(N1))
        self.centre2 = WireFrameVolume(surf2, paths(N3))
        self.centre3 = WireFrameVolume(surf3, paths(N4))
        self.centre4 = WireFrameVolume(surf4, paths(N5))
        self.centre5 = WireFrameVolume(surf5, paths(N6))
        self.centre6 = WireFrameVolume(surf6, paths(M7))
        self.centre7 = WireFrameVolume(surf7, paths(N7))
        self.centre8 = WireFrameVolume(surf8, paths(N6))
        self.centre9 = WireFrameVolume(surf9, paths(P3))
        self.centre10 = WireFrameVolume(surf10, paths(P2))
        self.centre11 = WireFrameVolume(surf11, paths(P1))
        self.centre12 = WireFrameVolume(surf12, paths(P0))


#------------------------------  End of File ---------------------------------#
