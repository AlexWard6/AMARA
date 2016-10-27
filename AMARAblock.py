# AMARAbin.py

"""
This file contains all the block definitions.

"""

#--------------------------------- FOREBODY ----------------------------------#
### NOSE TIP
#""    NOSE - TIP
#This block is the very central square block on the nosetip of the vehicle.
#""
#  (NORTH, EAST, SOUTH, WEST, TOP, BOTTOM)
tipBCs = [ExtrapolateOutBC(), ExtrapolateOutBC(), ExtrapolateOutBC(), ExtrapolateOutBC(), SupInBC(inflow), FixedTBC(Twall=wall_temp)]
#SlipWallBC(label='symmetry')
tipClust = [None,]*4 + \
           [None,]*4 + \
           [cf_wallNormal,]*4


# NE, SE, SW, NW

angle = np.pi/7.
phis0 = [phiTip, phiTip/2., phiTip/2.22, phiTip]
thetas0 = [angle_90-angle, angle_90-angle, angle_90, angle_90]

phis1 = [phiTip, phiTip/1.75, phis0[1], phiTip]
thetas1 = [angle_90-2.0*angle, angle_90-2.*angle, thetas0[1], thetas0[0]]

phis2 = [phiTip, phiTip, phiTip/2., phis1[1]]
thetas2 = [thetas1[0], angle/2., angle/2.0, thetas1[1]]

phis3 = [phiTip, phiTip, -1.*phis2[2], phis2[2]]
thetas3 = [thetas2[1], -1.*thetas2[1], angle_180-thetas2[2], thetas2[2]]

phis9 = [phiTip/4.3, phiTip/4.3, phiTip/5.0, phiTip/5.0]
thetas9 = [thetas1[0], thetas1[0]-3.*angle, -1*angle_90, angle_90]

phis10 = [phis3[3], phis3[2], phis9[1], phis9[0]]
thetas10 = [thetas3[3], thetas3[2], thetas9[1], thetas9[0]]

phis11 = [phis1[1], phis10[0], phis10[3], phis1[2]]
thetas11 = [thetas1[1], thetas10[0], thetas10[3], thetas1[2]]

phis12 = [phis1[2], phis11[2], phis9[3], phis0[2]]
thetas12 = [thetas1[2], thetas11[2], thetas9[3], thetas0[2]]

# NOSE 0
pyfunction_Nose_0 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis0, AR, thetas0)
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_0)
BLOCK_0t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (0)", 
                  nni=nCircumTail, nnj=nNormalWingTop, nnk=nRadialWall,
                  #bc_list=[ExtrapolateOutBC(), ExtrapolateOutBC(), ExtrapolateOutBC(), SlipWallBC(label='symmetry'), SupInBC(inflow), FixedTBC(Twall=wall_temp)],
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_0t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_0t.bc_list[TOP] = SupInBC(inflow)
BLOCK_0t.bc_list[EAST] = SlipWallBC(label='symmetry')


# NOSE 1
pyfunction_Nose_1 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis1, AR, thetas1)
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_1)
BLOCK_1t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (1)", 
                  nni=nCircumFuselage, nnj=nNormalWingTop, nnk=nRadialWall,
                  #bc_list=tipBCs,
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_1t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_1t.bc_list[TOP] = SupInBC(inflow)


# NOSE 2
pyfunction_Nose_2 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis2, AR, thetas2)
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_2)
BLOCK_2t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (2)", 
                  nni=nNormalWingTop, nnj=nNormalWingTop, nnk=nRadialWall,
                  bc_list=tipBCs,
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_2t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_2t.bc_list[TOP] = SupInBC(inflow)


# NOSE 3
pyfunction_Nose_3 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis3, AR, thetas3)
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_3)
BLOCK_3t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (3)", 
                  nni=nNormalWingTop, nnj=nWing, nnk=nRadialWall,
                  #bc_list=tipBCs,
                  cf_list = tipClust,
                 fill_condition = initial)
                 #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_3t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_3t.bc_list[TOP] = SupInBC(inflow)

# NOSE 4
pyfunction_Nose_4 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis2, AR, thetas2, 'BOTTOM')
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_4)
BLOCK_4t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (4)", 
                  nni=nNormalWingBottom, nnj=nNormalWingBottom, nnk=nRadialWall,
                  #bc_list=tipBCs,
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_4t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_4t.bc_list[TOP] = SupInBC(inflow)

# NOSE 5
pyfunction_Nose_5 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis1, AR, thetas1, 'BOTTOM')
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_5)
BLOCK_5t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (5)", 
                  nni=nCircumEngine, nnj=nNormalWingBottom, nnk=nRadialWall,
                  #bc_list=[ExtrapolateOutBC(), ExtrapolateOutBC(), ExtrapolateOutBC(), SlipWallBC(label='symmetry'), SupInBC(inflow), FixedTBC(Twall=wall_temp)],
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_5t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_5t.bc_list[TOP] = SupInBC(inflow)


# NOSE 6
pyfunction_Nose_6 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis0, AR, thetas0, 'BOTTOM')
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_6)
BLOCK_6t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (6)", 
                  nni=nCircumEngine, nnj=nNormalWingBottom, nnk=nRadialWall,
                  #bc_list=[ExtrapolateOutBC(), ExtrapolateOutBC(), ExtrapolateOutBC(), SlipWallBC(label='symmetry'), SupInBC(inflow), FixedTBC(Twall=wall_temp)],
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_6t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_6t.bc_list[TOP] = SupInBC(inflow)
BLOCK_6t.bc_list[EAST] = SlipWallBC(label='symmetry')

# NOSE 7
pyfunction_Nose_7 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis12, AR, thetas12, 'BOTTOM')
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_7)
BLOCK_7t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (7)", 
                  nni=nCircumEngine, nnj=nNormalWingTop, nnk=nRadialWall,
                  #bc_list=tipBCs,
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_7t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_7t.bc_list[TOP] = SupInBC(inflow)
BLOCK_7t.bc_list[EAST] = SlipWallBC(label='symmetry')

# NOSE 8
pyfunction_Nose_8 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis11, AR, thetas11, 'BOTTOM')
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_8)
BLOCK_8t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (8)", 
                  nni=nCircumFuselage, nnj=nNormalWingTop, nnk=nRadialWall,
                  #bc_list=tipBCs,
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_8t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_8t.bc_list[TOP] = SupInBC(inflow)              

# NOSE 9
pyfunction_Nose_9 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis9, AR, thetas9)
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_9)
BLOCK_9t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (9)", 
                  nni=nCircumTail, nnj=nWing, nnk=nRadialWall,
                  #bc_list=[ExtrapolateOutBC(), ExtrapolateOutBC(), ExtrapolateOutBC(), SlipWallBC(label='symmetry'), SupInBC(inflow), FixedTBC(Twall=wall_temp)],
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_9t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_9t.bc_list[TOP] = SupInBC(inflow)
BLOCK_9t.bc_list[EAST] = SlipWallBC(label='symmetry')

# NOSE 10
pyfunction_Nose_10 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis10, AR, thetas10)
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_10)
BLOCK_10t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (10)", 
                  nni=nCircumEngine, nnj=nWing, nnk=nRadialWall,
                  #bc_list=tipBCs,
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_10t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_10t.bc_list[TOP] = SupInBC(inflow)

# NOSE 11
pyfunction_Nose_11 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis11, AR, thetas11)
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_11)
BLOCK_11t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (11)", 
                  nni=nCircumFuselage, nnj=nNormalWingTop, nnk=nRadialWall,
                  #bc_list=tipBCs,
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])        
BLOCK_11t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_11t.bc_list[TOP] = SupInBC(inflow)

# NOSE 12
pyfunction_Nose_12 = lambda r,s,t: nose(r, s, t, r_nose, t_nose, phis12, AR, thetas12)
BLOCK_NOSE_VOL = PyFunctionVolume(pyfunction_Nose_12)
BLOCK_12t = Block3D(BLOCK_NOSE_VOL, label="NOSE-TIP (12)", 
                  nni=nCircumTail, nnj=nNormalWingTop, nnk=nRadialWall,
                  #bc_list=[ExtrapolateOutBC(), ExtrapolateOutBC(), ExtrapolateOutBC(), SlipWallBC(label='symmetry'), SupInBC(inflow), FixedTBC(Twall=wall_temp)],
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_12t.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_12t.bc_list[TOP] = SupInBC(inflow)
BLOCK_12t.bc_list[EAST] = SlipWallBC(label='symmetry')

### NOSE JOIN
#""    NOSE-JOIN
#These blocks are the final ring sitting on the spherical nose. They connect 
#to the central square blocks - RING_x. 
#"
joinBCs = bc_list=[ExtrapolateOutBC(), ExtrapolateOutBC(), ExtrapolateOutBC(),  ExtrapolateOutBC(), SupInBC(inflow), FixedTBC(Twall=wall_temp)]
joinClust = [None,]*4 + \
            [None,]*4 + \
            [cf_wallNormal,]*4
RingClust = [None,]*4 + \
            [None,]*4 + \
            [cf_wallNormal,]*4


# JOIN 0
# NE, SE, SW, NW
phis0 = [phiJoin, phis0[3], phis0[3], phiJoin]
thetas0 = [thetas0[0], thetas0[0], thetas0[3], thetas0[3]]
pyfunction_Join_0 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis0, AR, thetas0, block='0')
BLOCK_JOIN0_VOL = PyFunctionVolume(pyfunction_Join_0)
BLOCK_0j = Block3D(BLOCK_JOIN0_VOL, label="NOSE-JOIN (0)", 
                  nni=nCircumTail, nnj=nLongJoin, nnk=nRadialWall,
                  #bc_list=[ExtrapolateOutBC(), ExtrapolateOutBC(), ExtrapolateOutBC(), SlipWallBC(label='symmetry'), SupInBC(inflow), FixedTBC(Twall=wall_temp)],
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_0j.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_0j.bc_list[TOP] = SupInBC(inflow)
BLOCK_0j.bc_list[EAST] = SlipWallBC(label='symmetry')
pyfunction_Join_0_nocamber = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis0, AR, thetas0, camber=0, block='0')
BLOCK_JOIN0_VOL = PyFunctionVolume(pyfunction_Join_0_nocamber)


# JOIN 1
phis1 = [phiJoin, phis1[3], phis1[3], phiJoin]
thetas1 = [thetas1[0]-thetas2[0]/2., thetas1[0], thetas1[3], thetas1[3]]
pyfunction_Join_1 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis1, AR, thetas1, block='1')
BLOCK_JOIN1_VOL = PyFunctionVolume(pyfunction_Join_1)
BLOCK_1j = Block3D(BLOCK_JOIN1_VOL, label="NOSE-JOIN (1)", 
                  nni=nCircumFuselage, nnj=nLongJoin, nnk=nRadialWall,
                  #bc_list=tipBCs,
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_1j.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_1j.bc_list[TOP] = SupInBC(inflow)
pyfunction_Join_1 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis1, AR, thetas1, camber=0, block='1')
BLOCK_JOIN1_VOL = PyFunctionVolume(pyfunction_Join_1)

# JOIN 2
phis2 = [phiJoin, phiJoin, phis2[1], phis2[0]]
thetas2 = [thetas2[0]/2., thetas2[1]/2., thetas2[1], thetas2[0]]
pyfunction_Join_2 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis2, AR, thetas2, block='2')
BLOCK_JOIN2_VOL = PyFunctionVolume(pyfunction_Join_2)
BLOCK_2j = Block3D(BLOCK_JOIN2_VOL, label="NOSE-JOIN (2)", 
                  nni=nLongJoin, nnj=nNormalWingTop, nnk=nRadialWall,
                  #bc_list=tipBCs,
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_2j.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_2j.bc_list[TOP] = SupInBC(inflow)       
pyfunction_Join_2 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis2, AR, thetas2, camber=0, block='2')
BLOCK_JOIN2_VOL = PyFunctionVolume(pyfunction_Join_2)

# JOIN 3
phis3 = [phiJoin, phiJoin, phis3[1], phis3[0]]
thetas3 = [thetas3[0]/2., thetas3[1]/2., thetas3[1], thetas3[0]]
pyfunction_Join_3 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis3, AR, thetas3, face='BOTTOM', block='3')
BLOCK_JOIN3_VOL = PyFunctionVolume(pyfunction_Join_3)
BLOCK_3j = Block3D(BLOCK_JOIN3_VOL, label="NOSE-JOIN (3)", 
                  nni=nLongJoin, nnj=nWing, nnk=nRadialWall,
                  #bc_list=tipBCs,
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_3j.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_3j.bc_list[TOP] = SupInBC(inflow)
pyfunction_Join_3 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis3, AR, thetas3, camber=0, face='BOTTOM', block='3')
BLOCK_JOIN3_VOL = PyFunctionVolume(pyfunction_Join_3)

# JOIN 4
pyfunction_Join_4 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis2, AR, thetas2, face='BOTTOM', block='4')
BLOCK_JOIN4_VOL = PyFunctionVolume(pyfunction_Join_4)
BLOCK_4j = Block3D(BLOCK_JOIN4_VOL, label="NOSE-JOIN (4)", 
                  nni=nLongJoin, nnj=nNormalWingBottom, nnk=nRadialWall,
                  #bc_list=tipBCs,
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_4j.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_4j.bc_list[TOP] = SupInBC(inflow)
pyfunction_Join_4 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis2, AR, thetas2, camber=0, face='BOTTOM', block='4')
BLOCK_JOIN4_VOL = PyFunctionVolume(pyfunction_Join_4)

# JOIN 5
pyfunction_Join_5 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis1, AR, thetas1, face='BOTTOM', block='5')
BLOCK_JOIN5_VOL = PyFunctionVolume(pyfunction_Join_5)
BLOCK_5j = Block3D(BLOCK_JOIN5_VOL, label="NOSE-JOIN (5)", 
                  nni=nCircumEngine, nnj=nLongJoin, nnk=nRadialWall,
                  #bc_list=[ExtrapolateOutBC(), ExtrapolateOutBC(), ExtrapolateOutBC(), SlipWallBC(label='symmetry'), SupInBC(inflow), FixedTBC(Twall=wall_temp)],
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_5j.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_5j.bc_list[TOP] = SupInBC(inflow)
pyfunction_Join_5 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis1, AR, thetas1, camber=0, face='BOTTOM', block='5')
BLOCK_JOIN5_VOL = PyFunctionVolume(pyfunction_Join_5)

# JOIN 6
pyfunction_Join_6 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis0, AR, thetas0, face='BOTTOM', block='6')
BLOCK_JOIN6_VOL = PyFunctionVolume(pyfunction_Join_6)
BLOCK_6j = Block3D(BLOCK_JOIN6_VOL, label="NOSE-JOIN (5)", 
                  nni=nCircumEngine, nnj=nLongJoin, nnk=nRadialWall,
                  #bc_list=[ExtrapolateOutBC(), ExtrapolateOutBC(), ExtrapolateOutBC(), SlipWallBC(label='symmetry'), SupInBC(inflow), FixedTBC(Twall=wall_temp)],
                  cf_list = tipClust,
                  fill_condition = initial)
                  #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_6j.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_6j.bc_list[TOP] = SupInBC(inflow)
BLOCK_6j.bc_list[EAST] = SlipWallBC(label='symmetry')
pyfunction_Join_6 = lambda r,s,t: noseJoin(r, s, t, r_nose, t_nose, phis0, AR, thetas0, camber=0, face='BOTTOM', block='6')
BLOCK_JOIN6_VOL = PyFunctionVolume(pyfunction_Join_6)

### FOREBODY RAMP
#" FOREBODY-RAMP
#These blocks forming the ramp of the forebody. At the nose they are defined 
#by the mating face of the corresponding Join blocks. Along the length of the 
#forebody ramp, the blocks polar angle change to correctly form the blocks 
#wrapping on the wings, engines and tail.
#""
rampClustN = [None, cf_wallNormal,]*2 + \
             [None, cf_wallNormal,]*2 + \
             [None,]*4

rampClustS = [None, cf_wallNormal,]*2 + \
             [None, cf_wallNormal,]*2 + \
             [None,]*4

# Tail, right side
nRampSubdivide = 1
print "RAMP_0"
matingFace = evalFace(BLOCK_JOIN0_VOL, 'N')
pyfunction_Ramp_4 = lambda r,s,t: ramp(r, s, t, matingFace, 'S', r_join, t_ramp, dr_ramp, phi_mesh, L_ramp, angle_180-angle_Tail, angle_90, angle_180-angle_Tail, angle_90, phiJoin)
BLOCK_25_VOL = PyFunctionVolume(pyfunction_Ramp_4)
BLOCK_0r = Block3D(BLOCK_25_VOL, label="RAMP-4 (25)", 
                   nni=nCircumTail, nnj=nRadialWall, nnk=nLongRamp,
                   #nbi=1, nbj=1, nbk=nRampSubdivide,
                   #bc_list=rampBCs,
                   cf_list = rampClustS,
                   fill_condition = initial)
                   #xforce_list = [0, 0, 1, 0, 0, 0])
#for blk in BLOCK_0r.blks[0][0]:
#    blk.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)
#    blk.bc_list[NORTH] = SupInBC(inflow)
#    blk.bc_list[EAST] = SlipWallBC(label='symmetry')
#BLOCK_0r.bc_list[BOTTOM] = ExtrapolateOutBC()
BLOCK_0r.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)
BLOCK_0r.bc_list[NORTH] = SupInBC(inflow)
BLOCK_0r.bc_list[EAST] = SlipWallBC(label='symmetry')


# Fuselage, right side
print "RAMP_1"
matingFace = evalFace(BLOCK_JOIN1_VOL, 'N')  
pyfunction_Ramp_5 = lambda r,s,t: ramp(r, s, t, matingFace, '', r_join, t_ramp, dr_ramp, phi_mesh, L_ramp, angle_180-angle_WingTop, angle_180-angle_Tail, angle_180-angle_WingTop, angle_180-angle_Tail, phiJoin)
BLOCK_26_VOL = PyFunctionVolume(pyfunction_Ramp_5)
BLOCK_1r = Block3D(BLOCK_26_VOL, label="RAMP-5 (26)", 
                   nni=nCircumEngine, nnj=nRadialWall, nnk=nLongRamp,
                   #nbi=1, nbj=1, nbk=nRampSubdivide,
                   #bc_list=rampBCs,
                   cf_list = rampClustN,
                   fill_condition = initial)
                   #xforce_list = [0, 0, 1, 0, 0, 0])
#for blk in BLOCK_0r.blks[0][0]:
#    blk.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)
#    blk.bc_list[NORTH] = SupInBC(inflow)
#BLOCK_1r.bc_list[BOTTOM] = ExtrapolateOutBC()
BLOCK_1r.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)
BLOCK_1r.bc_list[NORTH] = SupInBC(inflow)

# Right wing topside
print "RAMP_2"
matingFace = evalFace(BLOCK_JOIN2_VOL, 'W')
pyfunction_Ramp_6 = lambda r,s,t: ramp(r, s, t, matingFace, '', r_join, t_ramp, dr_ramp, phi_mesh, L_ramp, angle_180, angle_180-angle_WingTop, angle_180, angle_180-angle_WingTop, phiJoin, block='2')
BLOCK_27_VOL = PyFunctionVolume(pyfunction_Ramp_6)
BLOCK_2r = Block3D(BLOCK_27_VOL, label="RAMP-6 (27)", 
                   nni=nNormalWingTop, nnj=nRadialWall, nnk=nLongRamp,
                   #nbi=1, nbj=1, nbk=nRampSubdivide,
                   #bc_list=rampBCs,
                   cf_list = rampClustN,
                   fill_condition = initial)
                   #xforce_list = [0, 0, 1, 0, 0, 0])
#for blk in BLOCK_0r.blks[0][0]:
#    blk.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)
#    blk.bc_list[NORTH] = SupInBC(inflow)
#BLOCK_2r.bc_list[BOTTOM] = ExtrapolateOutBC()
BLOCK_2r.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)
BLOCK_2r.bc_list[NORTH] = SupInBC(inflow)

# Wing
print "RAMP_3"
matingFace = evalFace(BLOCK_JOIN3_VOL, 'W')
pyfunction_Ramp_7= lambda r,s,t: ramp(r, s, t, matingFace, 'S', r_join, t_ramp, dr_ramp, phi_mesh, L_ramp, angle_180+angle_WingBottom, angle_180, angle_180+angle_WingBottom, angle_180, phiJoin, block='3')
BLOCK_28_VOL = PyFunctionVolume(pyfunction_Ramp_7)
BLOCK_3r = Block3D(BLOCK_28_VOL, label="RAMP-7 (28)", 
                   nni=nWing, nnj=nRadialWall, nnk=nLongRamp,
                   #nbi=1, nbj=1, nbk=nRampSubdivide,
                   #bc_list=rampBCs,
                   cf_list = rampClustS,
                   fill_condition = initial)
                   #xforce_list = [0, 0, 1, 0, 0, 0])
#for blk in BLOCK_0r.blks[0][0]:
#    blk.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)
#    blk.bc_list[NORTH] = SupInBC(inflow)
#BLOCK_3r.bc_list[BOTTOM] = ExtrapolateOutBC()
BLOCK_3r.bc_list[NORTH] = SupInBC(inflow)
BLOCK_3r.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)

# Right wing underside
print "RAMP_4"
matingFace = evalFace(BLOCK_JOIN4_VOL, 'W')
pyfunction_Ramp_8 = lambda r,s,t: ramp(r, s, t, matingFace, 'S', r_join, t_ramp, dr_ramp, phi_mesh, L_ramp, angle_180+angle_WingBottom, angle_Engine4+inletSmile/2, angle_180+angle_WingBottom, angle_Engine4+inletSmile/2, phiJoin, block='4')
BLOCK_29_VOL = PyFunctionVolume(pyfunction_Ramp_8)
BLOCK_4r = Block3D(BLOCK_29_VOL, label="RAMP-8 (29)", 
                   nni=nNormalWingBottom, nnj=nRadialWall, nnk=nLongRamp,
                   #nbi=1, nbj=1, nbk=nRampSubdivide,
                   #bc_list=rampBCs,
                   cf_list = rampClustN,
                   fill_condition = initial)
                   #xforce_list = [0, 0, 1, 0, 0, 0])
#for blk in BLOCK_0r.blks[0][0]:
#    blk.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)
#    blk.bc_list[NORTH] = SupInBC(inflow)
#BLOCK_4r.bc_list[BOTTOM] = ExtrapolateOutBC()
BLOCK_4r.bc_list[NORTH] = SupInBC(inflow)
BLOCK_4r.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)

# Engine 4
print "RAMP_5"
matingFace = evalFace(BLOCK_JOIN5_VOL, 'S')  
pyfunction_Ramp_9 = lambda r,s,t: ramp(r, s, t, matingFace, 'N', r_join, t_ramp, dr_ramp, phi_mesh, L_ramp, angle_180+angle_WingBottom, angle_Engine4+inletSmile/2, angle_180+angle_WingBottom, angle_Engine4+inletSmile/2, phiJoin)
BLOCK_30_VOL = PyFunctionVolume(pyfunction_Ramp_9)
BLOCK_5r = Block3D(BLOCK_30_VOL, label="RAMP-9 (30)", 
                   nni=nCircumEngine, nnj=nRadialWall, nnk=nLongRamp,
                   #nbi=1, nbj=1, nbk=nRampSubdivide,
                   #bc_list=[SupInBC(inflow), ExtrapolateOutBC(), ExtrapolateOutBC(), FixedTBC(Twall=wall_temp), None, ExtrapolateOutBC()],
                   cf_list = rampClustN,
                   fill_condition = initial)
                   #xforce_list = [0, 0, 1, 0, 0, 0])
#for blk in BLOCK_0r.blks[0][0]:
#    blk.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)
#    blk.bc_list[NORTH] = SupInBC(inflow)
#BLOCK_5r.bc_list[TOP] = ExtrapolateOutBC()                   
BLOCK_5r.bc_list[NORTH] = SupInBC(inflow)
BLOCK_5r.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)

# Engine 3
print "RAMP_6"
matingFace = evalFace(BLOCK_JOIN6_VOL, 'S')  
pyfunction_Ramp_9 = lambda r,s,t: ramp(r, s, t, matingFace, 'N', r_join, t_ramp, dr_ramp, phi_mesh, L_ramp, angle_Engine4+inletSmile/2, angle_Engine3+inletSmile/2, angle_Engine4+inletSmile/2, angle_Engine3+inletSmile/2, phiJoin)
BLOCK_31_VOL = PyFunctionVolume(pyfunction_Ramp_9)
BLOCK_6r = Block3D(BLOCK_31_VOL, label="RAMP-9 (30)", 
                   nni=nCircumEngine, nnj=nRadialWall, nnk=nLongRamp,
                   #nbi=1, nbj=1, nbk=nRampSubdivide,
                   #bc_list=[SupInBC(inflow), ExtrapolateOutBC(), ExtrapolateOutBC(), FixedTBC(Twall=wall_temp), None, ExtrapolateOutBC()],
                   cf_list = rampClustN,
                   fill_condition = initial)
                   #xforce_list = [0, 0, 1, 0, 0, 0])
#for blk in BLOCK_0r.blks[0][0]:
#    blk.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)
#    blk.bc_list[NORTH] = SupInBC(inflow)
#    blk.bc_list[EAST] = SlipWallBC(label='symmetry')
#BLOCK_6r.bc_list[TOP] = ExtrapolateOutBC()    
BLOCK_6r.bc_list[NORTH] = SupInBC(inflow)
BLOCK_6r.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)
BLOCK_6r.bc_list[EAST] = SlipWallBC(label='symmetry')


#------------------------------ TOP SURFACE --------------------------------#
#bc_list=[NORTH, EAST, SOUTH, WEST, TOP, BOTTOM]
rOuterStart = rInner_inlet+t_ramp + L_ramp*np.tan(phi_mesh)
outer_inletBlocks = outerBlocks(angle_WingTop, zStartInlet, zStartInletRamp, rInner_inlet, rOuterStart, rInner_inletRamp, side='RHS')


### ROOT FAIRING (0)
# Right wing, topside
rootFairing_RHwingTop = rootFairing(L_inlet, zStartInlet, rInner_inlet, phi_fb, (90.*np.pi/180-wallAngle_inlet), angle_WingTop, side='RHS')
rootFairing_flatPaths = rootFairing_RHwingTop.flatPaths
BLOCK_0 = Block3D(rootFairing_RHwingTop.volume, label="BLOCK-0", 
	    	  nni=nLongInlet, nnj=nLongInlet, nnk=nNormalWingTop,
	    	  fill_condition=initial)
	          #xforce_list = [0, 0, 0, 1, 0, 1])
BLOCK_0.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_0.bc_list[WEST] = FixedTBC(Twall=wall_temp)

outerinlet_wingTop = outer_inletBlocks.topsideBlock()
wingTop_flatpaths = outer_inletBlocks.flatPaths
BLOCK_M3 = Block3D(outerinlet_wingTop, label="BLOCK-M3", 
                   nni=nRadialWall , nnj=nLongInlet, nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

# Fuselage, right side
rootFairing_RHfuselage = rootFairing_RHwingTop.rotate(-1*angle_WingTop, -1*angle_Tail, rootFairing_flatPaths)
BLOCK_0 = Block3D(rootFairing_RHfuselage, label="BLOCK-0", 
	    	  nni=nLongInlet, nnj=nLongInlet, nnk=nCircumFuselage,
	    	  fill_condition=initial)
                  #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_0.bc_list[WEST] = FixedTBC(Twall=wall_temp)

outerinlet_fuselage = outer_inletBlocks.rotate(-1*angle_WingTop, -1*angle_Tail, wingTop_flatpaths)
fuselage_flatpaths = outer_inletBlocks.flatPaths
BLOCK_M4 = Block3D(outerinlet_fuselage, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongInlet, nnk=nCircumFuselage,
                   fill_condition=initial)        
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)

# Tail, right side
rootFairing_RHtail = rootFairing_RHwingTop.rotate(-1*angle_Tail, -1*angle_90, rootFairing_flatPaths)
BLOCK_0 = Block3D(rootFairing_RHtail, label="BLOCK-0", 
	    	  nni=nLongInlet, nnj=nLongInlet, nnk=nCircumTail,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_0.bc_list[TOP] = SlipWallBC(label='symmetry')
BLOCK_0.bc_list[WEST] = FixedTBC(Twall=wall_temp)

outerinlet_tail = outer_inletBlocks.rotate(-1*angle_Tail, -1*angle_90, fuselage_flatpaths)
BLOCK_M5 = Block3D(outerinlet_tail, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongInlet, nnk=nCircumFuselage,
                   fill_condition=initial)     
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[TOP] = SlipWallBC(label='symmetry')
# LE Rounding
#edgeRoundingTop = edgeRounding('RHS')
#edgeRounding20 = edgeRoundingTop.block_LEfairing(r_leadingEdge, rootFairing_RHwingTop.surfacePath0, rootFairing_RHwingTop.surfacePath_outer0, (90.*np.pi/180-wallAngle_inlet))
#centreEnd20 = edgeRoundingTop.centreEnd
#centreEnd_outer20 = edgeRoundingTop.centreEnd_outer
#BLOCK_20 = Block3D(edgeRounding20, label="BLOCK-20 - LE FAIRING ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nLongInlet,
#	    	  fill_condition=initial)

#edgeRoundingTop = edgeRounding('RHS')
#
#edgeRounding20 = edgeRoundingTop.block_LEfairing(r_leadingEdge, rootFairing_RHwingTop.surfacePath0, rootFairing_RHwingTop.surfacePath_outer0, (90.*np.pi/180-wallAngle_inlet))
#centreEnd20 = edgeRoundingTop.centreEnd
#centreEnd_outer20 = edgeRoundingTop.centreEnd_outer



### FORE DIAMONDS - INLET RAMP (1-2)
foreDiamonds_RHwingTop = foreDiamonds(zStartCombustor, 1, angle_WingTop, 'RHS')

# Right wing, topside
matingFace1 = rootFairing_RHwingTop.southFace
BLOCK_1 = Block3D(foreDiamonds_RHwingTop.foreDiamonds_1(matingFace1), label="BLOCK-1", 
	    	  nni=nLongInlet, nnj=nLongInletRamp, nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 1])
BLOCK_1.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_1.bc_list[WEST] = FixedTBC(Twall=wall_temp)

matingFace2 = rootFairing_RHwingTop.eastFace
BLOCK_2 = Block3D(foreDiamonds_RHwingTop.foreDiamonds_2(matingFace2), label="BLOCK-2", 
	    	  nni=nLongInlet, nnj=nLongInletRamp, nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_2.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

rOuterStart = outer_inletBlocks.rOuterEnd
rInnerEnd = foreDiamonds_RHwingTop.rOuterEnd
outer_inletRampBlocks = outerBlocks(angle_WingTop, zStartInletRamp, zStartCombustor, rInner_inletRamp, rOuterStart, rInnerEnd, side='RHS')
outerInletRamp_wingTop = outer_inletRampBlocks.topsideBlock()
wingTop_flatpaths = outer_inletRampBlocks.flatPaths
BLOCK_M3 = Block3D(outerInletRamp_wingTop, label="BLOCK-M3", 
                   nni=nRadialWall, nnj=nLongInletRamp, nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

# Fuselage, right side
foreDiamonds_RHfuselage_1, combustorTop_fuselage1 = foreDiamonds_RHwingTop.rotate(-1*angle_WingTop, -1*angle_Tail, foreDiamonds_RHwingTop.flatPaths1)
BLOCK_0 = Block3D(foreDiamonds_RHfuselage_1, label="BLOCK-0", 
	    	  nni=nLongInlet, nnj=nLongInletRamp, nnk=nCircumFuselage,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_0.bc_list[WEST] = FixedTBC(Twall=wall_temp)

foreDiamonds_RHfuselage_2, combustorTop_fuselage2 = foreDiamonds_RHwingTop.rotate(-1*angle_WingTop, -1*angle_Tail, foreDiamonds_RHwingTop.flatPaths2)
BLOCK_0 = Block3D(foreDiamonds_RHfuselage_2, label="BLOCK-0", 
	    	  nni=nLongInlet, nnj=nLongInletRamp, nnk=nCircumFuselage,
	    	  fill_condition=initial)

outerInletRamp_fuselage = outer_inletRampBlocks.rotate(-1*angle_WingTop, -1*angle_Tail, wingTop_flatpaths)
fuselage_flatpaths = outer_inletRampBlocks.flatPaths
BLOCK_M4 = Block3D(outerInletRamp_fuselage, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongInletRamp, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)

# Tail, right side
foreDiamonds_RHtail_1, combustorTop_tail1 = foreDiamonds_RHwingTop.rotate(-1*angle_Tail, -1*angle_90, foreDiamonds_RHwingTop.flatPaths1)
BLOCK_0 = Block3D(foreDiamonds_RHtail_1, label="BLOCK-0", 
	    	  nni=nLongInlet, nnj=nLongInletRamp, nnk=nCircumTail,
	    	  fill_condition=initial)
BLOCK_0.bc_list[WEST] = FixedTBC(Twall=wall_temp)
BLOCK_0.bc_list[TOP] = SlipWallBC(label='symmetry')

foreDiamonds_RHtail_2, combustorTop_tail2 = foreDiamonds_RHwingTop.rotate(-1*angle_Tail, -1*angle_90, foreDiamonds_RHwingTop.flatPaths2)
BLOCK_0 = Block3D(foreDiamonds_RHtail_2, label="BLOCK-0", 
	    	  nni=nLongInlet, nnj=nLongInletRamp, nnk=nCircumTail,
	    	  fill_condition=initial)
BLOCK_0.bc_list[TOP] = SlipWallBC(label='symmetry')

outerInletRamp_tail = outer_inletRampBlocks.rotate(-1*angle_Tail, -1*angle_90, fuselage_flatpaths)
BLOCK_M5 = Block3D(outerInletRamp_tail, label="BLOCK-M5", 
                   nni=nRadialWall , nnj=nLongInletRamp, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[TOP] = SlipWallBC(label='symmetry')
# LE Rounding
#print 'surfacePath_outer', foreDiamonds_RHwingTop.surfacePath_outer2
#print 'surfacePath', foreDiamonds_RHwingTop.surfacePath2
#print 'centreEnd', centreEnd
#print 'centreEnd_outer', centreEnd_outer
#edgeRounding21 = edgeRoundingTop.block_LE(r_leadingEdge, foreDiamonds_RHwingTop.surfacePath2, foreDiamonds_RHwingTop.surfacePath_outer2, centreEnd20, centreEnd_outer20, angle_WingTop, theta_leadingEdge=90.*np.pi/180.)
#centreEnd21 = edgeRoundingTop.centreEnd
#centreEnd_outer21 = edgeRoundingTop.centreEnd_outer
#BLOCK_21 = Block3D(edgeRounding21, label="BLOCK-21 - LE FORE DIAMONDS ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nLongInletRamp,
#	    	  fill_condition=initial)

### DIAMONDS (5-7)
diamonds_RHwingTop = diamonds(angle_WingTop, diamondRotation, side='RHS')

# Right wing, topside
BLOCK_5 = Block3D(diamonds_RHwingTop.diamond5(), label="BLOCK-5", 
	    	  nni=nElevonNeighbour , nnj=nLongInlet , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 0, 1])
westFace5_wingTop = diamonds_RHwingTop.westFace5
southFace5_wingTop = diamonds_RHwingTop.southFace5
BLOCK_5.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

BLOCK_6 = Block3D(diamonds_RHwingTop.diamond6(), label="BLOCK-6", 
	    	  nni=nLongInlet , nnj=nElevonNeighbour , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 0, 1])
northFace6_wingTop = diamonds_RHwingTop.northFace6
eastFace6_wingTop = diamonds_RHwingTop.eastFace6
BLOCK_6.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

BLOCK_7 = Block3D(diamonds_RHwingTop.diamond7(), label="BLOCK-7", 
	    	  nni=nLongInlet , nnj=nLongInlet , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 0, 1])
southFace7_wingTop = diamonds_RHwingTop.southFace7
eastFace7_wingTop = diamonds_RHwingTop.eastFace7
BLOCK_7.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

# Fuselage, right side
BLOCK_5 = Block3D(diamonds_RHwingTop.rotate(-1*angle_WingTop, -1*angle_Tail, diamonds_RHwingTop.flatPaths5), label="BLOCK-5", 
	    	  nni=nElevonNeighbour , nnj=nLongInlet , nnk=nCircumFuselage,
	    	  fill_condition=initial)
westFace5_fuselage = diamonds_RHwingTop.westFace5
southFace5_fuselage = diamonds_RHwingTop.southFace5

BLOCK_6 = Block3D(diamonds_RHwingTop.rotate(-1*angle_WingTop, -1*angle_Tail, diamonds_RHwingTop.flatPaths6), label="BLOCK-6", 
	    	  nni=nLongInlet , nnj=nElevonNeighbour , nnk=nCircumFuselage,
	    	  fill_condition=initial)
northFace6_fuselage = diamonds_RHwingTop.northFace6
eastFace6_fuselage = diamonds_RHwingTop.eastFace6


BLOCK_7 = Block3D(diamonds_RHwingTop.rotate(-1*angle_WingTop, -1*angle_Tail, diamonds_RHwingTop.flatPaths7), label="BLOCK-7", 
	    	  nni=nLongInlet , nnj=nLongInlet , nnk=nCircumFuselage,
	    	  fill_condition=initial)
southFace7_fuselage = diamonds_RHwingTop.southFace7
eastFace7_fuselage = diamonds_RHwingTop.eastFace7

     
# Tail, right side
BLOCK_5 = Block3D(diamonds_RHwingTop.rotate(-1*angle_Tail, -1*angle_90, diamonds_RHwingTop.flatPaths5), label="BLOCK-5", 
	    	  nni=nElevonNeighbour, nnj=nLongInlet, nnk=nCircumTail,
	    	  fill_condition=initial)
westFace5_tail = diamonds_RHwingTop.westFace5
southFace5_tail = diamonds_RHwingTop.southFace5
BLOCK_5.bc_list[TOP] = SlipWallBC(label='symmetry')

BLOCK_6 = Block3D(diamonds_RHwingTop.rotate(-1*angle_Tail, -1*angle_90, diamonds_RHwingTop.flatPaths6), label="BLOCK-6", 
	    	  nni=nLongInlet, nnj=nElevonNeighbour, nnk=nCircumTail,
	    	  fill_condition=initial)
northFace6_tail = diamonds_RHwingTop.northFace6
eastFace6_tail = diamonds_RHwingTop.eastFace6
BLOCK_6.bc_list[TOP] = SlipWallBC(label='symmetry')

BLOCK_7 = Block3D(diamonds_RHwingTop.rotate(-1*angle_Tail, -1*angle_90, diamonds_RHwingTop.flatPaths7), label="BLOCK-7", 
	    	  nni=nLongInlet, nnj=nLongInlet, nnk=nCircumTail,
	    	  fill_condition=initial)
southFace7_tail = diamonds_RHwingTop.southFace7
eastFace7_tail = diamonds_RHwingTop.eastFace7
BLOCK_7.bc_list[TOP] = SlipWallBC(label='symmetry')

### FORE DIAMONDS - COMBUSTOR (3-4)
# Right wing, topside
topFace = foreDiamonds_RHwingTop.combustorTop3 #evalFace(BLOCK_1_VOL, 'S', rsReversal=1)
bottomFace = westFace5_wingTop #diamonds_RHwingTop.westFace5 #evalFace(BLOCK_5_VOL, 'W')
pyfunction_fore_3 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_3_VOL = PyFunctionVolume(pyfunction_fore_3)
BLOCK_3 = Block3D(BLOCK_3_VOL, label="BLOCK-3", 
	    	  nni=nLongInlet , nnj=nNormalWingTop , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 1, 1, 0, 0])
BLOCK_3.bc_list[WEST] = FixedTBC(Twall=wall_temp)
BLOCK_3.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)

topFace = foreDiamonds_RHwingTop.combustorTop4 #= evalFace(BLOCK_6_VOL, 'N')
bottomFace = northFace6_wingTop #diamonds_RHwingTop.northFace6 #evalFace(BLOCK_2_VOL, 'S')
pyfunction_fore_4 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_4_VOL = PyFunctionVolume(pyfunction_fore_4)
BLOCK_4 = Block3D(BLOCK_4_VOL, label="BLOCK-4", 
	    	  nni=nLongInlet , nnj=nNormalWingTop , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [1, 0, 0, 0, 0, 1])  
BLOCK_4.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)

rOuterStart = outer_inletRampBlocks.rOuterEnd
zDiamonds = diamonds_RHwingTop.Nodetwoz
outer_foreDiamonds2Blocks = outerBlocks(angle_WingTop, zStartCombustor, zDiamonds, rInnerEnd, rOuterStart, diamonds_RHwingTop.rOuter, side='RHS')
outerforeDiamonds2_wingTop = outer_foreDiamonds2Blocks.topsideBlock()
wingTop_flatpaths = outer_foreDiamonds2Blocks.flatPaths
BLOCK_M3 = Block3D(outerforeDiamonds2_wingTop, label="BLOCK-M3", 
                   nni=nRadialWall, nnj=int(nLongCombust/2), nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

# Fuselage, right side
topFace = combustorTop_fuselage1  #evalFace(BLOCK_1_VOL, 'S', rsReversal=1)
bottomFace = westFace5_fuselage #diamonds_RHwingTop.westFace5       #evalFace(BLOCK_5_VOL, 'W')
pyfunction_fore_3 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_3_VOL = PyFunctionVolume(pyfunction_fore_3)
BLOCK_3 = Block3D(BLOCK_3_VOL, label="BLOCK-3", 
	    	  nni=nLongInlet , nnj=nCircumFuselage , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_3.bc_list[WEST] = FixedTBC(Twall=wall_temp)

topFace = combustorTop_fuselage2 #= evalFace(BLOCK_6_VOL, 'N')
bottomFace = northFace6_fuselage #diamonds_RHwingTop.northFace6 #evalFace(BLOCK_2_VOL, 'S')
pyfunction_fore_4 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_4_VOL = PyFunctionVolume(pyfunction_fore_4)
BLOCK_4 = Block3D(BLOCK_4_VOL, label="BLOCK-4", 
	    	  nni=nLongInlet , nnj=nCircumFuselage, nnk=int(nLongCombust/2),
	    	  fill_condition=initial)  

foreDiamonds2_fuselage = outer_foreDiamonds2Blocks.rotate(-1*angle_WingTop, -1*angle_Tail, wingTop_flatpaths)
fuselage_flatpaths = outer_foreDiamonds2Blocks.flatPaths
BLOCK_M4 = Block3D(foreDiamonds2_fuselage, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=int(nLongCombust/2), nnk=nCircumFuselage,
                   fill_condition=initial)   
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)

# Tail, right side
topFace = combustorTop_tail1  #evalFace(BLOCK_1_VOL, 'S', rsReversal=1)
bottomFace = westFace5_tail #diamonds_RHwingTop.westFace5       #evalFace(BLOCK_5_VOL, 'W')
pyfunction_fore_3 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_3_VOL = PyFunctionVolume(pyfunction_fore_3)
BLOCK_3 = Block3D(BLOCK_3_VOL, label="BLOCK-3", 
	    	  nni=nLongInlet , nnj=nCircumTail, nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_3.bc_list[WEST] = FixedTBC(Twall=wall_temp)
BLOCK_3.bc_list[NORTH] = SlipWallBC(label='symmetry')

topFace = combustorTop_tail2 #= evalFace(BLOCK_6_VOL, 'N')
bottomFace = northFace6_tail #diamonds_RHwingTop.northFace6 #evalFace(BLOCK_2_VOL, 'S')
pyfunction_fore_4 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_4_VOL = PyFunctionVolume(pyfunction_fore_4)
BLOCK_4 = Block3D(BLOCK_4_VOL, label="BLOCK-4", 
	    	  nni=nLongInlet , nnj=nCircumTail, nnk=int(nLongCombust/2),
	    	  fill_condition=initial)  
BLOCK_4.bc_list[NORTH] = SlipWallBC(label='symmetry')

foreDiamonds2_tail = outer_foreDiamonds2Blocks.rotate(-1*angle_Tail, -1*angle_90, wingTop_flatpaths)
BLOCK_M5 = Block3D(foreDiamonds2_tail, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=int(nLongCombust/2), nnk=nCircumTail,
                   fill_condition=initial)  
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[TOP] = SlipWallBC(label='symmetry')
# LE Rounding
#surfacePath_outer4 = Line(northFace6_wingTop.eval(1, 1), foreDiamonds_RHwingTop.combustorTop4.eval(1, 1))
#surfacePath4 = Line(northFace6_wingTop.eval(1, 0), foreDiamonds_RHwingTop.combustorTop4.eval(1, 0))

#print 'surfacePath_outer', surfacePath_outer4, 'surfacePath', surfacePath4, 'centreEnd', centreEnd, 'centreEnd_outer', centreEnd_outer

#BLOCK_22_VOL, centreEnd, centreEnd_outer = edgeRounding(r_leadingEdge, surfacePath4, surfacePath_outer4, centreEnd, centreEnd_outer, angle_WingTop, theta_leadingEdge=90.*np.pi/180, side='RHS')
#edgeRounding22 = edgeRoundingTop.block_LE(r_leadingEdge, surfacePath4, surfacePath_outer4, centreEnd21, centreEnd_outer21, angle_WingTop, theta_leadingEdge=90.*np.pi/180,)
#centreEnd22 = edgeRoundingTop.centreEnd
#centreEnd_outer22 = edgeRoundingTop.centreEnd_outer
#BLOCK_22 = Block3D(edgeRounding22, label="BLOCK-22 - LE COMBUSTOR ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=int(nLongCombust/2),
#	    	  fill_condition=initial)  


### OUTLET - FORE (12-15)
rears_wingTop = rearBlocks(0., angle_WingTop, 'RHS')

# Right wing, topside
rears_wingTop_12 = rears_wingTop.block(L_outlet-L_elevon, (xStartElevon - r_fb), zStartOutlet, r_fb, '12')
rears_flatPaths12 = rears_wingTop.flatPaths
BLOCK_12 = Block3D(rears_wingTop_12, label="BLOCK-12", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletFore , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 1])
BLOCK_12.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_12.bc_list[WEST] = FixedTBC(Twall=wall_temp)

rears_wingTop_13 = rears_wingTop.block(L_outlet-L_elevon, W_elevon/2., zStartOutlet, xStartElevon, '13')
rears_flatPaths13 = rears_wingTop.flatPaths
BLOCK_13 = Block3D(rears_wingTop_13, label="BLOCK-13", 
	    	  nni=nLongInlet, nnj=nLongOutletFore , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_13.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

rears_wingTop_14 = rears_wingTop.block(L_outlet-L_elevon, W_elevon/2., zStartOutlet, xStartElevon+W_elevon/2., '14')
rears_flatPaths14 = rears_wingTop.flatPaths
BLOCK_14 = Block3D(rears_wingTop_14, label="BLOCK-14", 
	    	  nni=nLongInlet, nnj=nLongOutletFore , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 0, 1])      
BLOCK_14.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

rears_wingTop_15 = rears_wingTop.block(L_outlet-L_elevon, '', zStartOutlet, xStartElevon+W_elevon, 'EDGE')
rears_flatPaths15 = rears_wingTop.flatPaths
surfacePath15 = rears_wingTop.pE
surfacePath_outer15 = rears_wingTop.pEtop
BLOCK_15 = Block3D(rears_wingTop_15, label="BLOCK-15", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletFore , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 0, 1])        
BLOCK_15.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

rOuterStart = outer_foreDiamonds2Blocks.rOuterEnd
rInnerEnd = rears_wingTop.rOuterStart
zDiamonds = diamonds_RHwingTop.Nodetwoz
rInner_diamonds = diamonds_RHwingTop.rOuter
outer_aftDiamondsBlocks = outerBlocks(angle_WingTop, zDiamonds, zStartOutlet, rInner_diamonds, rOuterStart, rInnerEnd, side='RHS')
outeraftDiamonds_wingTop = outer_aftDiamondsBlocks.topsideBlock()
wingTop_flatpaths = outer_aftDiamondsBlocks.flatPaths
BLOCK_M3 = Block3D(outeraftDiamonds_wingTop, label="BLOCK-M3", 
                   nni=nRadialWall , nnj=int(nLongCombust/2), nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

# Fuselage, right side
rears_fuselage_12 = rears_wingTop.rotate(-1*angle_WingTop, -1*angle_Tail, rears_flatPaths12)
BLOCK_12 = Block3D(rears_fuselage_12, label="BLOCK-12", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletFore , nnk=nCircumFuselage,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_12.bc_list[WEST] = FixedTBC(Twall=wall_temp)

rears_fuselage_13 = rears_wingTop.rotate(-1*angle_WingTop, -1*angle_Tail, rears_flatPaths13)
BLOCK_13 = Block3D(rears_fuselage_13, label="BLOCK-13", 
	    	  nni=nLongInlet, nnj=nLongOutletFore , nnk=nCircumFuselage,
	    	  fill_condition=initial)
   
rears_fuselage_14 = rears_wingTop.rotate(-1*angle_WingTop, -1*angle_Tail, rears_flatPaths14)
BLOCK_14 = Block3D(rears_fuselage_14, label="BLOCK-14", 
	    	  nni=nLongInlet, nnj=nLongOutletFore , nnk=nCircumFuselage,
	    	  fill_condition=initial)      

rears_fuselage_15 = rears_wingTop.rotate(-1*angle_WingTop, -1*angle_Tail, rears_flatPaths15)
BLOCK_15 = Block3D(rears_fuselage_15, label="BLOCK-15", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletFore , nnk=nCircumFuselage,
	    	  fill_condition=initial)  

outeraftDiamonds_fuselage = outer_aftDiamondsBlocks.rotate(-1*angle_WingTop, -1*angle_Tail, wingTop_flatpaths)
BLOCK_M4 = Block3D(outeraftDiamonds_fuselage, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=int(nLongCombust/2), nnk=nCircumEngine,
                   fill_condition=initial) 
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)

# Tail, right side
rears_tail_12 = rears_wingTop.rotate(-1*angle_Tail, -1*angle_90, rears_flatPaths12)
BLOCK_12 = Block3D(rears_tail_12, label="BLOCK-12", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletFore , nnk=nCircumTail,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_12.bc_list[WEST] = FixedTBC(Twall=wall_temp)
BLOCK_12.bc_list[TOP] = SlipWallBC(label='symmetry')

rears_tail_13 = rears_wingTop.rotate(-1*angle_Tail, -1*angle_90, rears_flatPaths13)
BLOCK_13 = Block3D(rears_tail_13, label="BLOCK-13", 
	    	  nni=nLongInlet, nnj=nLongOutletFore , nnk=nCircumTail,
	    	  fill_condition=initial)
BLOCK_13.bc_list[TOP] = SlipWallBC(label='symmetry')
   
rears_tail_14 = rears_wingTop.rotate(-1*angle_Tail, -1*angle_90, rears_flatPaths14)
BLOCK_14 = Block3D(rears_tail_14, label="BLOCK-14", 
	    	  nni=nLongInlet, nnj=nLongOutletFore , nnk=nCircumTail,
	    	  fill_condition=initial)      
BLOCK_14.bc_list[TOP] = SlipWallBC(label='symmetry')

rears_tail_15 = rears_wingTop.rotate(-1*angle_Tail, -1*angle_90, rears_flatPaths15)
BLOCK_15 = Block3D(rears_tail_15, label="BLOCK-15", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletFore , nnk=nCircumTail,
	    	  fill_condition=initial)  
BLOCK_15.bc_list[TOP] = SlipWallBC(label='symmetry')

outeraftDiamonds_fuselage = outer_aftDiamondsBlocks.rotate(-1*angle_Tail, -1*angle_90, wingTop_flatpaths)
BLOCK_M5 = Block3D(outeraftDiamonds_fuselage, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=int(nLongCombust/2), nnk=nCircumTail,
                   fill_condition=initial) 
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[TOP] = SlipWallBC(label='symmetry')
# LE ROUNDING DEFINED BELOW IN AFT OF DIAMONDS ROW


### AFT DIAMONDS - COMBUSTOR (8-11)
# Right wing, topside
topFace = southFace5_wingTop
bottomFace = evalFace(rears_wingTop_12, 'N')
pyfunction_rear_8 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_8_VOL = PyFunctionVolume(pyfunction_rear_8)
BLOCK_8 = Block3D(BLOCK_8_VOL, label="BLOCK-8", 
	    	  nni=nElevonNeighbour , nnj=nNormalWingTop , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 1, 1, 0, 0])
BLOCK_8.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)
BLOCK_8.bc_list[WEST] = FixedTBC(Twall=wall_temp)

topFace = southFace7_wingTop
bottomFace = evalFace(rears_wingTop_13, 'N')
pyfunction_rear_9 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_9_VOL = PyFunctionVolume(pyfunction_rear_9)
BLOCK_9 = Block3D(BLOCK_9_VOL, label="BLOCK-9", 
	    	  nni=nLongInlet, nnj=nNormalWingTop , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 1, 0, 0, 0])
BLOCK_9.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)

topFace = eastFace7_wingTop
bottomFace = evalFace(rears_wingTop_14, 'N')
pyfunction_rear_10 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_10_VOL = PyFunctionVolume(pyfunction_rear_10)
BLOCK_10 = Block3D(BLOCK_10_VOL, label="BLOCK-10", 
	    	  nni=nLongInlet, nnj=nNormalWingTop , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 1, 0, 0, 0])
BLOCK_10.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)

topFace = eastFace6_wingTop
northFace15 = evalFace(rears_wingTop_15, 'N')
pyfunction_rear_11 = lambda r,s,t: surf2surf(r, s, t, topFace, northFace15)
BLOCK_11_VOL = PyFunctionVolume(pyfunction_rear_11)
BLOCK_11 = Block3D(BLOCK_11_VOL, label="BLOCK-11", 
	    	  nni=nElevonNeighbour , nnj=nNormalWingTop , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 1, 0, 0, 0])
BLOCK_11.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)

rOuterStart = outer_aftDiamondsBlocks.rOuterEnd
rInnerEnd = rears_wingTop.rOuterEnd
rInnerStart = rears_wingTop.rOuterStart
outer_elevon1Blocks = outerBlocks(angle_WingTop, zStartOutlet, zStartElevon, rInnerStart, rOuterStart, rInnerEnd, side='RHS')
outer_elevon1_wingBottom = outer_elevon1Blocks.topsideBlock()
wingTop_flatpaths = outer_elevon1Blocks.flatPaths
BLOCK_M3 = Block3D(outer_elevon1_wingBottom, label="BLOCK-M3", 
                   nni=nRadialWall , nnj=nLongOutletFore, nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

# Right fuselage
topFace = southFace5_fuselage
bottomFace = evalFace(rears_fuselage_12, 'N')
pyfunction_rear_8 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_8_VOL = PyFunctionVolume(pyfunction_rear_8)
BLOCK_8 = Block3D(BLOCK_8_VOL, label="BLOCK-8", 
	    	  nni=nElevonNeighbour , nnj=nCircumFuselage , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_8.bc_list[WEST] = FixedTBC(Twall=wall_temp)

topFace = southFace7_fuselage
bottomFace = evalFace(rears_fuselage_13, 'N')
pyfunction_rear_9 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_9_VOL = PyFunctionVolume(pyfunction_rear_9)
BLOCK_9 = Block3D(BLOCK_9_VOL, label="BLOCK-9", 
	    	  nni=nLongInlet, nnj=nCircumFuselage , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)

topFace = eastFace7_fuselage
bottomFace = evalFace(rears_fuselage_14, 'N')
pyfunction_rear_10 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_10_VOL = PyFunctionVolume(pyfunction_rear_10)
BLOCK_10 = Block3D(BLOCK_10_VOL, label="BLOCK-10", 
	    	  nni=nLongInlet, nnj=nCircumFuselage , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)

topFace = eastFace6_fuselage
bottomFace = evalFace(rears_fuselage_15, 'N')
pyfunction_rear_11 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_11_VOL = PyFunctionVolume(pyfunction_rear_11)
BLOCK_11 = Block3D(BLOCK_11_VOL, label="BLOCK-11", 
	    	  nni=nElevonNeighbour , nnj=nCircumFuselage , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
 
outerelevon1_fuselage = outer_elevon1Blocks.rotate(-1*angle_WingTop, -1*angle_Tail,wingTop_flatpaths)
BLOCK_M4 = Block3D(outerelevon1_fuselage, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongOutletFore, nnk=nCircumFuselage,
                   fill_condition=initial) 
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)

# Right tail
topFace = southFace5_tail
bottomFace = evalFace(rears_tail_12, 'N')
pyfunction_rear_8 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_8_VOL = PyFunctionVolume(pyfunction_rear_8)
BLOCK_8 = Block3D(BLOCK_8_VOL, label="BLOCK-8", 
	    	  nni=nElevonNeighbour , nnj=nCircumTail , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_8.bc_list[NORTH] = SlipWallBC(label='symmetry')
BLOCK_8.bc_list[WEST] = FixedTBC(Twall=wall_temp)

topFace = southFace7_tail
bottomFace = evalFace(rears_tail_13, 'N')
pyfunction_rear_9 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_9_VOL = PyFunctionVolume(pyfunction_rear_9)
BLOCK_9 = Block3D(BLOCK_9_VOL, label="BLOCK-9", 
	    	  nni=nLongInlet, nnj=nCircumTail , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
BLOCK_9.bc_list[NORTH] = SlipWallBC(label='symmetry')

topFace = eastFace7_tail
bottomFace = evalFace(rears_tail_14, 'N')
pyfunction_rear_10 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_10_VOL = PyFunctionVolume(pyfunction_rear_10)
BLOCK_10 = Block3D(BLOCK_10_VOL, label="BLOCK-10", 
	    	  nni=nLongInlet, nnj=nCircumTail , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
BLOCK_10.bc_list[NORTH] = SlipWallBC(label='symmetry')

topFace = eastFace6_tail
bottomFace = evalFace(rears_tail_15, 'N')
pyfunction_rear_11 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_11_VOL = PyFunctionVolume(pyfunction_rear_11)
BLOCK_11 = Block3D(BLOCK_11_VOL, label="BLOCK-11", 
	    	  nni=nElevonNeighbour , nnj=nCircumTail , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
BLOCK_11.bc_list[NORTH] = SlipWallBC(label='symmetry')

outerelevon1_tail = outer_elevon1Blocks.rotate(-1*angle_Tail, -1*angle_90, wingTop_flatpaths)
BLOCK_M5 = Block3D(outerelevon1_tail, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongOutletFore, nnk=nCircumTail,
                   fill_condition=initial) 
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[TOP] = SlipWallBC(label='symmetry')

# LE ROUNDING - AFT OF DIAMONDS ROW
surfacePath11 = Line(northFace15.eval(1, 0), eastFace6_wingTop.eval(1, 0))
surfacePath_outer11 = Line(northFace15.eval(1, 1), eastFace6_wingTop.eval(1, 1))

#print 'surfacePath_outer', surfacePath_outer4
#print 'surfacePath', surfacePath4
#print 'centreEnd', centreEnd
#print 'centreEnd_outer', centreEnd_outer

#edgeRounding23 = edgeRoundingTop.block_LE(r_leadingEdge, surfacePath11, surfacePath_outer11, centreEnd22, centreEnd_outer22, angle_WingTop, theta_leadingEdge=90.*np.pi/180.)
#centreEnd23 = edgeRoundingTop.centreEnd
#centreEnd_outer23 = edgeRoundingTop.centreEnd_outer
#BLOCK_23 = Block3D(edgeRounding23, label="BLOCK-23 - AFT DIAMONDS ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=int(nLongCombust/2),
#	    	  fill_condition=initial)


# LE Rounding FORE OUTLET ROW
#BLOCK_24_VOL, centreEnd, centreEnd_outer = edgeRounding(r_leadingEdge, surfacePath15, surfacePath_outer15, centreEnd, centreEnd_outer, angle_WingTop, theta_leadingEdge=90.*np.pi/180, side='RHS')
#edgeRounding24 = edgeRoundingTop.block_LE(r_leadingEdge, surfacePath15, surfacePath_outer15, centreEnd23, centreEnd_outer23, angle_WingTop, theta_leadingEdge=90.*np.pi/180.)
#centreEnd24 = edgeRoundingTop.centreEnd
#centreEnd_outer24 = edgeRoundingTop.centreEnd_outer
#BLOCK_24 = Block3D(edgeRounding24, label="BLOCK-24 - LE COMBUSTOR ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nLongOutletFore,
#	    	  fill_condition=initial) 


### OUTLET - AFT - ELEVONS (16-19)
rears_wingTop_19 = rears_wingTop.block(L_elevon, '', zStartElevon, xStartElevon+W_elevon, 'EDGE', TEflag=0)
surfacePath19 = rears_wingTop.pE; surfacePath_outer19 = rears_wingTop.pEtop; TE_surfacePath19 = rears_wingTop.pS; TE_surfacePath_outer19 = rears_wingTop.pStop
rears_flatPaths19 = rears_wingTop.flatPaths
BLOCK_19 = Block3D(rears_wingTop_19, label="BLOCK-19", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletElevon , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_19.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

# LE Rounding AFT OUTLET ROW
#BLOCK_25_VOL, centreEnd, centreEnd_outer = edgeRounding(r_leadingEdge, surfacePath19, surfacePath_outer19, centreEnd, centreEnd_outer, angle_WingTop, theta_leadingEdge=90.*np.pi/180, endFlag=1, side='RHS')
#edgeRounding25 = edgeRoundingTop.block_LE(r_leadingEdge, surfacePath19, surfacePath_outer19, centreEnd24, centreEnd_outer24, angle_WingTop, theta_leadingEdge=90.*np.pi/180, endFlag=1)
#centreEnd25 = edgeRoundingTop.centreEnd
#centreEnd_outer25 = edgeRoundingTop.centreEnd_outer
#BLOCK_25 = Block3D(edgeRounding25, label="BLOCK-25 - LE ELEVON ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nLongOutletElevon,
#	    	  fill_condition=initial) 

#BLOCK_26_VOL, centreEnd, centreEnd_outer = edgeRounding(r_leadingEdge, TE_surfacePath19, TE_surfacePath_outer19, centreEnd, centreEnd_outer, angle_WingTop, theta_leadingEdge=90.*np.pi/180, endFlag=1, TEFlag=1, side='RHS')
#edgeRounding26 = edgeRoundingTop.block_LE(r_leadingEdge, TE_surfacePath19, TE_surfacePath_outer19, centreEnd25, centreEnd_outer25, angle_WingTop, theta_leadingEdge=90.*np.pi/180, endFlag=1, TEFlag=1)
#centreEnd26 = edgeRoundingTop.centreEnd
#centreEnd_outer26 = edgeRoundingTop.centreEnd_outer
#BLOCK_26 = Block3D(edgeRounding26, label="BLOCK-26 - TE ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nElevonNeighbour,
#	    	  fill_condition=initial) 



rears_wingTop_18 = rears_wingTop.block(L_elevon, W_elevon/2., zStartElevon, xStartElevon+W_elevon/2., '18', TEflag=0)
TE_surfacePath18 = rears_wingTop.pS; TE_surfacePath_outer18 = rears_wingTop.pStop; 
rears_flatPaths18 = rears_wingTop.flatPaths
BLOCK_18 = Block3D(rears_wingTop_18, label="BLOCK-18", 
	    	  nni=nLongInlet, nnj=nLongOutletElevon , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_18.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

#edgeRounding27 = edgeRoundingTop.block_LE(r_leadingEdge, TE_surfacePath18, TE_surfacePath_outer18, centreEnd26, centreEnd_outer26, angle_WingTop, theta_leadingEdge=90.*np.pi/180, endFlag=1, TEFlag=1)
#centreEnd27 = edgeRoundingTop.centreEnd
#centreEnd_outer27 = edgeRoundingTop.centreEnd_outer
#BLOCK_27 = Block3D(edgeRounding27, label="BLOCK-27 - TE ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nLongInlet,
#	    	  fill_condition=initial) 




rears_wingTop_17 = rears_wingTop.block(L_elevon, W_elevon/2., zStartElevon, xStartElevon, '17', TEflag=0)
TE_surfacePath17 = rears_wingTop.pS; TE_surfacePath_outer17 = rears_wingTop.pStop
rears_flatPaths17 = rears_wingTop.flatPaths
BLOCK_17 = Block3D(rears_wingTop_17, label="BLOCK-17", 
	    	  nni=nLongInlet, nnj=nLongOutletElevon , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_17.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
  
#edgeRounding28 = edgeRoundingTop.block_LE(r_leadingEdge, TE_surfacePath17, TE_surfacePath_outer17, centreEnd27, centreEnd_outer27, angle_WingTop, theta_leadingEdge=90.*np.pi/180, endFlag=1, TEFlag=1)
#centreEnd28 = edgeRoundingTop.centreEnd
#centreEnd_outer28 = edgeRoundingTop.centreEnd_outer
#BLOCK_28 = Block3D(edgeRounding28, label="BLOCK-28 - TE ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nLongInlet,
#	    	  fill_condition=initial) 



rears_wingTop_16 = rears_wingTop.block(L_elevon, (xStartElevon - r_fb), zStartElevon, r_fb, '16', TEflag=0)
TE_surfacePath16 = rears_wingTop.pS; TE_surfacePath_outer16 = rears_wingTop.pStop
rears_flatPaths16 = rears_wingTop.flatPaths
BLOCK_16 = Block3D(rears_wingTop_16, label="BLOCK-16", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletElevon , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 1])
BLOCK_16.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_16.bc_list[WEST] = FixedTBC(Twall=wall_temp)


#edgeRounding29 = edgeRoundingTop.block_LE(r_leadingEdge, TE_surfacePath16, TE_surfacePath_outer16, centreEnd28, centreEnd_outer28, angle_WingTop, theta_leadingEdge=90.*np.pi/180, endFlag=1, TEFlag=1)
#centreEnd29 = edgeRoundingTop.centreEnd
#centreEnd_outer29 = edgeRoundingTop.centreEnd_outer
#BLOCK_29 = Block3D(edgeRounding29, label="BLOCK-29 - TE ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nElevonNeighbour,
#	    	  fill_condition=initial)

rears_wingTop.block(L_elevon, '', zStartElevon, xStartElevon+W_elevon, 'EDGE', TEflag=0)
rOuterStart = outer_elevon1Blocks.rOuterEnd
rInnerEnd = rears_wingTop.rOuterEnd
rInnerStart = rears_wingTop.rOuterStart
outer_elevon2Blocks = outerBlocks(angle_WingTop, zStartElevon, zStartboatTail, rInnerStart, rOuterStart, rInnerEnd, side='RHS')
outer_elevon2_wingTop = outer_elevon2Blocks.topsideBlock()
wingTop_flatpaths = outer_elevon2Blocks.flatPaths
BLOCK_M3 = Block3D(outer_elevon2_wingTop, label="BLOCK-M3", 
                   nni=nRadialWall , nnj=nLongOutletElevon, nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

# Fuselage, right side
rears_fuselage_16 = rears_wingTop.rotate(-1*angle_WingTop, -1*angle_Tail, rears_flatPaths16)
BLOCK_16 = Block3D(rears_fuselage_16, label="BLOCK-16", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletElevon , nnk=nCircumFuselage,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_16.bc_list[WEST] = FixedTBC(Twall=wall_temp)

rears_fuselage_17 = rears_wingTop.rotate(-1*angle_WingTop, -1*angle_Tail, rears_flatPaths17)
BLOCK_17 = Block3D(rears_fuselage_17, label="BLOCK-17", 
	    	  nni=nLongInlet, nnj=nLongOutletElevon , nnk=nCircumFuselage,
	    	  fill_condition=initial)
        
        
rears_fuselage_18 = rears_wingTop.rotate(-1*angle_WingTop, -1*angle_Tail, rears_flatPaths18)
BLOCK_18 = Block3D(rears_fuselage_18, label="BLOCK-18", 
	    	  nni=nLongInlet, nnj=nLongOutletElevon , nnk=nCircumFuselage,
	    	  fill_condition=initial)  

rears_fuselage_19 = rears_wingTop.rotate(-1*angle_WingTop, -1*angle_Tail, rears_flatPaths19)
BLOCK_19 = Block3D(rears_fuselage_19, label="BLOCK-19", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletElevon , nnk=nCircumFuselage,
	    	  fill_condition=initial)

outerelevon2_fuselage = outer_elevon2Blocks.rotate(-1*angle_WingTop, -1*angle_Tail,wingTop_flatpaths)
BLOCK_M4 = Block3D(outerelevon2_fuselage, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongOutletElevon, nnk=nCircumFuselage,
                   fill_condition=initial) 
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)

# Tail, right side
rears_tail_16 = rears_wingTop.rotate(-1*angle_Tail, -1*angle_90, rears_flatPaths16)
BLOCK_16 = Block3D(rears_tail_16, label="BLOCK-16", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletElevon , nnk=nCircumTail,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_16.bc_list[WEST] = FixedTBC(Twall=wall_temp)        
BLOCK_16.bc_list[TOP] = SlipWallBC(label='symmetry')

rears_tail_17 = rears_wingTop.rotate(-1*angle_Tail, -1*angle_90, rears_flatPaths17)
BLOCK_17 = Block3D(rears_tail_17, label="BLOCK-17", 
	    	  nni=nLongInlet, nnj=nLongOutletElevon , nnk=nCircumTail,
	    	  fill_condition=initial)
BLOCK_17.bc_list[TOP] = SlipWallBC(label='symmetry')

rears_tail_18 = rears_wingTop.rotate(-1*angle_Tail, -1*angle_90, rears_flatPaths18)
BLOCK_18 = Block3D(rears_tail_18, label="BLOCK-18", 
	    	  nni=nLongInlet, nnj=nLongOutletElevon, nnk=nCircumTail,
	    	  fill_condition=initial)  
BLOCK_18.bc_list[TOP] = SlipWallBC(label='symmetry')

rears_tail_19 = rears_wingTop.rotate(-1*angle_Tail, -1*angle_90, rears_flatPaths19)
BLOCK_19 = Block3D(rears_tail_19, label="BLOCK-19", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletElevon, nnk=nCircumTail,
	    	  fill_condition=initial)
BLOCK_19.bc_list[TOP] = SlipWallBC(label='symmetry')
    
outerelevon2_fuselage = outer_elevon2Blocks.rotate(-1*angle_Tail, -1*angle_90, wingTop_flatpaths)
BLOCK_M5 = Block3D(outerelevon2_fuselage, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongOutletElevon, nnk=nCircumTail,
                   fill_condition=initial) 
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[TOP] = SlipWallBC(label='symmetry')

### BOAT TAIL (20-23)
boatTailBlocks = BTandWake(0., angle_WingTop, 'RHS', face='TOP')
boatTail_wingTop_20 = boatTailBlocks.block(L_boatTail, (xStartElevon - rInner_boatTail), zStartboatTail, rInner_boatTail, wallAngle=wallAngle_boatTail, block='WALL')
boatTail_flatPaths20 = boatTailBlocks.flatPaths
combustorBottom20_wingTop = boatTailBlocks.northFace

BLOCK_20 = Block3D(boatTail_wingTop_20, label="BLOCK-20", 
	    	  nni=nElevonNeighbour , nnj=nLongBoatTail , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_20.bc_list[WEST] = FixedTBC(Twall=wall_temp)

boatTail_wingTop_21 = boatTailBlocks.block(L_boatTail, W_elevon/2., zStartboatTail, xStartElevon, block='21')
boatTail_flatPaths21 = boatTailBlocks.flatPaths
combustorBottom21_wingTop = boatTailBlocks.northFace
BLOCK_21 = Block3D(boatTail_wingTop_21, label="BLOCK-20", 
	    	  nni=nLongInlet , nnj=nLongBoatTail , nnk=nNormalWingTop,
	    	  fill_condition=initial)
  
boatTail_wingTop_22 = boatTailBlocks.block(L_boatTail, W_elevon/2., zStartboatTail, xStartElevon+W_elevon/2., block='22')
boatTail_flatPaths22 = boatTailBlocks.flatPaths
combustorBottom22_wingTop = boatTailBlocks.northFace
BLOCK_22 = Block3D(boatTail_wingTop_22, label="BLOCK-22", 
	    	  nni=nLongInlet , nnj=nLongBoatTail , nnk=nNormalWingTop,
	    	  fill_condition=initial)      

boatTail_wingTop_23 = boatTailBlocks.block(L_boatTail, '', zStartboatTail, xStartElevon+W_elevon, meshAngle=5.*np.pi/180, block='EDGE')
boatTail_flatPaths23 = boatTailBlocks.flatPaths
combustorBottom23_wingTop = boatTailBlocks.northFace
surfacePath23 = boatTailBlocks.pE
surfacePath_outer23 = boatTailBlocks.pEtop
BLOCK_23 = Block3D(boatTail_wingTop_23, label="BLOCK-23", 
                   nni=nElevonNeighbour , nnj=nLongBoatTail , nnk=nNormalWingTop,
                   fill_condition=initial)    

rOuterStart = outer_elevon2Blocks.rOuterEnd
rInnerEnd = boatTailBlocks.rOuterEnd
rInnerStart = boatTailBlocks.rOuterStart
outer_boatTailBlocks = outerBlocks(angle_WingTop, zStartboatTail, zStartWake, rInnerStart, rOuterStart, rInnerEnd, side='RHS')
outer_boatTail_wingTop = outer_boatTailBlocks.topsideBlock()
wingTop_flatpaths = outer_boatTailBlocks.flatPaths
BLOCK_M3 = Block3D(outer_boatTail_wingTop, label="BLOCK-M3", 
                   nni=nRadialWall, nnj=nLongBoatTail, nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

# Fuselage, right side
boatTail_fuselage_20 = boatTailBlocks.rotate(-1*angle_WingTop, -1*angle_Tail, boatTail_flatPaths20, block='WALL')
BLOCK_20 = Block3D(boatTail_fuselage_20, label="BLOCK-20", 
                   nni=nElevonNeighbour, nnj=nLongBoatTail, nnk=nCircumFuselage,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_20.bc_list[WEST] = FixedTBC(Twall=wall_temp)

boatTail_fuselage_21 = boatTailBlocks.rotate(-1*angle_WingTop, -1*angle_Tail, boatTail_flatPaths21)
BLOCK_21 = Block3D(boatTail_fuselage_21, label="BLOCK-21", 
                   nni=nLongInlet, nnj=nLongBoatTail, nnk=nCircumFuselage,
                   fill_condition=initial)
        
        
boatTail_fuselage_22 = boatTailBlocks.rotate(-1*angle_WingTop, -1*angle_Tail, boatTail_flatPaths22)
BLOCK_22 = Block3D(boatTail_fuselage_22, label="BLOCK-2", 
                   nni=nLongInlet, nnj=nLongBoatTail, nnk=nCircumFuselage,
                   fill_condition=initial)  

boatTail_fuselage_23 = boatTailBlocks.rotate(-1*angle_WingTop, -1*angle_Tail, boatTail_flatPaths23)
BLOCK_23= Block3D(boatTail_fuselage_23, label="BLOCK-23", 
                  nni=nElevonNeighbour , nnj=nLongBoatTail , nnk=nCircumFuselage,
                  fill_condition=initial)

outerboatTail_fuselage = outer_boatTailBlocks.rotate(-1*angle_WingTop, -1*angle_Tail,wingTop_flatpaths)
BLOCK_M4 = Block3D(outerboatTail_fuselage, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongBoatTail, nnk=nCircumFuselage,
                   fill_condition=initial) 
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)

# Tail, right side
boatTail_tail_20 = boatTailBlocks.rotate(-1*angle_Tail, -1*angle_90, boatTail_flatPaths20, block='WALL')
BLOCK_20 = Block3D(boatTail_tail_20, label="BLOCK-20", 
                   nni=nElevonNeighbour, nnj=nLongBoatTail, nnk=nCircumTail,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_20.bc_list[TOP] = SlipWallBC(label='symmetry')
BLOCK_20.bc_list[WEST] = FixedTBC(Twall=wall_temp)

boatTail_tail_21 = boatTailBlocks.rotate(-1*angle_Tail, -1*angle_90, boatTail_flatPaths21)
BLOCK_21 = Block3D(boatTail_tail_21, label="BLOCK-21", 
                   nni=nLongInlet, nnj=nLongBoatTail , nnk=nCircumTail,
                   fill_condition=initial)
BLOCK_21.bc_list[TOP] = SlipWallBC(label='symmetry')

boatTail_tail_22 = boatTailBlocks.rotate(-1*angle_Tail, -1*angle_90, boatTail_flatPaths22)
BLOCK_22 = Block3D(boatTail_tail_22, label="BLOCK-22", 
                   nni=nLongInlet, nnj=nLongBoatTail , nnk=nCircumTail,
                   fill_condition=initial)  
BLOCK_22.bc_list[TOP] = SlipWallBC(label='symmetry')

boatTail_tail_23 = boatTailBlocks.rotate(-1*angle_Tail, -1*angle_90, boatTail_flatPaths23)
BLOCK_23 = Block3D(boatTail_tail_23, label="BLOCK-23", 
                   nni=nElevonNeighbour , nnj=nLongBoatTail , nnk=nCircumTail,
                   fill_condition=initial)
BLOCK_23.bc_list[TOP] = SlipWallBC(label='symmetry')


outerboatTail_tail = outer_boatTailBlocks.rotate(-1*angle_Tail, -1*angle_90, wingTop_flatpaths)
BLOCK_M5 = Block3D(outerboatTail_tail, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongBoatTail, nnk=nCircumTail,
                   fill_condition=initial) 
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[TOP] = SlipWallBC(label='symmetry')

### WAKE (24-27)
# Right wing, top side
wakeBlocks = BTandWake(0., angle_WingTop, 'RHS', face='TOP')
wake_wingTop_24 = wakeBlocks.block(L_aftFarField, (xStartElevon - r_bt), zStartWake, r_bt, wallAngle=0., block='WALL')
wake_flatPaths24 = wakeBlocks.flatPaths
combustorBottom24_wingTop = wakeBlocks.northFace
BLOCK_24 = Block3D(wake_wingTop_24, label="BLOCK-24", 
	    	  nni=nElevonNeighbour, nnj=nLongWake, nnk=nNormalWingTop,
	    	  fill_condition=initial)
BLOCK_24.bc_list[SOUTH] = ExtrapolateOutBC()

wake_wingTop_25 = wakeBlocks.block(L_aftFarField, W_elevon/2., zStartWake, xStartElevon, block='25')
wake_flatPaths25 = wakeBlocks.flatPaths
combustorBottom25_wingTop = wakeBlocks.northFace
BLOCK_25 = Block3D(wake_wingTop_25, label="BLOCK-25", 
	    	  nni=nLongInlet, nnj=nLongWake, nnk=nNormalWingTop,
	    	  fill_condition=initial)
BLOCK_25.bc_list[SOUTH] = ExtrapolateOutBC()

wake_wingTop_26 = wakeBlocks.block(L_aftFarField, W_elevon/2., zStartWake, xStartElevon+W_elevon/2., block='26')
wake_flatPaths26 = wakeBlocks.flatPaths
combustorBottom26_wingTop = wakeBlocks.northFace
BLOCK_26 = Block3D(wake_wingTop_26, label="BLOCK-26", 
	    	  nni=nLongInlet, nnj=nLongWake, nnk=nNormalWingTop,
	    	  fill_condition=initial)      
BLOCK_26.bc_list[SOUTH] = ExtrapolateOutBC()

wake_wingTop_27 = wakeBlocks.block(L_aftFarField, '', zStartWake, xStartElevon+W_elevon, meshAngle=5.*np.pi/180, block='EDGE')
wake_flatPaths27 = wakeBlocks.flatPaths
combustorBottom27_wingTop = wakeBlocks.northFace
surfacePath27 = wakeBlocks.pE
surfacePath_outer27 = boatTailBlocks.pEtop
BLOCK_27 = Block3D(wake_wingTop_27, label="BLOCK-27", 
                   nni=nElevonNeighbour, nnj=nLongWake , nnk=nNormalWingTop,
                   fill_condition=initial)   
BLOCK_27.bc_list[SOUTH] = ExtrapolateOutBC()

rOuterStart = outer_boatTailBlocks.rOuterEnd
rInnerEnd = wakeBlocks.rOuterEnd
rInnerStart = wakeBlocks.rOuterStart
outer_wakeBlocks = outerBlocks(angle_WingTop, zStartWake, zStartWake-L_aftFarField, rInnerStart, rOuterStart, rInnerEnd, side='RHS')
outer_wake_wingTop = outer_wakeBlocks.topsideBlock()
wingTop_flatpaths = outer_wakeBlocks.flatPaths
BLOCK_M3 = Block3D(outer_wake_wingTop, label="BLOCK-M3", 
                   nni=nRadialWall, nnj=nLongWake, nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)
BLOCK_M3.bc_list[SOUTH] = ExtrapolateOutBC()

# Fuselage, right side
wake_fuselage_24 = wakeBlocks.rotate(-1*angle_WingTop, -1*angle_Tail, wake_flatPaths24, block='WALL')
BLOCK_24 = Block3D(wake_fuselage_24, label="BLOCK-24", 
                   nni=nElevonNeighbour, nnj=nLongWake, nnk=nCircumFuselage,
                   fill_condition=initial)
BLOCK_24.bc_list[SOUTH] = ExtrapolateOutBC()

wake_fuselage_25 = wakeBlocks.rotate(-1*angle_WingTop, -1*angle_Tail, wake_flatPaths25)
BLOCK_25 = Block3D(wake_fuselage_25, label="BLOCK-25", 
                   nni=nLongInlet, nnj=nLongWake, nnk=nCircumFuselage,
                   fill_condition=initial)
BLOCK_25.bc_list[SOUTH] = ExtrapolateOutBC()
    
wake_fuselage_26 = wakeBlocks.rotate(-1*angle_WingTop, -1*angle_Tail, wake_flatPaths26)
BLOCK_26 = Block3D(wake_fuselage_26, label="BLOCK-26", 
                   nni=nLongInlet, nnj=nLongWake, nnk=nCircumFuselage,
                   fill_condition=initial)  
BLOCK_26.bc_list[SOUTH] = ExtrapolateOutBC()

wake_fuselage_27 = wakeBlocks.rotate(-1*angle_WingTop, -1*angle_Tail, wake_flatPaths27)
BLOCK_27 = Block3D(wake_fuselage_27, label="BLOCK-27",
                  nni=nElevonNeighbour, nnj=nLongWake, nnk=nCircumFuselage,
                  fill_condition=initial)
BLOCK_27.bc_list[SOUTH] = ExtrapolateOutBC()

outer_wake_fuselage = outer_wakeBlocks.rotate(-1*angle_WingTop, -1*angle_Tail, wingTop_flatpaths)
BLOCK_M4 = Block3D(outer_wake_fuselage, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongWake, nnk=nCircumTail,
                   fill_condition=initial) 
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)
BLOCK_M4.bc_list[SOUTH] = ExtrapolateOutBC()

# Tail, right side
wake_tail_24 = wakeBlocks.rotate(-1*angle_Tail, -1*angle_90, wake_flatPaths24, block='WALL')
BLOCK_24 = Block3D(wake_tail_24, label="BLOCK-24", 
                   nni=nElevonNeighbour, nnj=nLongWake, nnk=nCircumTail,
                   fill_condition=initial)
BLOCK_24.bc_list[TOP] = SlipWallBC(label='symmetry')
BLOCK_24.bc_list[SOUTH] = ExtrapolateOutBC()

wake_tail_25 = wakeBlocks.rotate(-1*angle_Tail, -1*angle_90, wake_flatPaths25)
BLOCK_25 = Block3D(wake_tail_25, label="BLOCK-25", 
                   nni=nLongInlet, nnj=nLongWake, nnk=nCircumTail,
                   fill_condition=initial)
BLOCK_25.bc_list[TOP] = SlipWallBC(label='symmetry')
BLOCK_25.bc_list[SOUTH] = ExtrapolateOutBC()

wake_tail_26 = wakeBlocks.rotate(-1*angle_Tail, -1*angle_90, wake_flatPaths26)
BLOCK_26 = Block3D(wake_tail_26, label="BLOCK-26", 
                   nni=nLongInlet, nnj=nLongWake, nnk=nCircumTail,
                   fill_condition=initial)  
BLOCK_26.bc_list[TOP] = SlipWallBC(label='symmetry')
BLOCK_26.bc_list[SOUTH] = ExtrapolateOutBC()

wake_tail_27 = wakeBlocks.rotate(-1*angle_Tail, -1*angle_90, wake_flatPaths27)
BLOCK_27 = Block3D(wake_tail_27, label="BLOCK-27", 
                   nni=nElevonNeighbour , nnj=nLongWake, nnk=nCircumTail,
                   fill_condition=initial)
BLOCK_27.bc_list[TOP] = SlipWallBC(label='symmetry')
BLOCK_27.bc_list[SOUTH] = ExtrapolateOutBC()

outer_wake_tail = outer_wakeBlocks.rotate(-1*angle_Tail, -1*angle_90, wingTop_flatpaths)
BLOCK_M5 = Block3D(outer_wake_tail, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongWake, nnk=nCircumTail,
                   fill_condition=initial) 
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[TOP] = SlipWallBC(label='symmetry')
BLOCK_M5.bc_list[SOUTH] = ExtrapolateOutBC()



#------------------------------ BOTTOM SURFACE --------------------------------#
### FORE DIAMONDS - INLET
rOuterStart = rInner_inlet+t_ramp + L_ramp*np.tan(phi_mesh)
outer_inletBlocks = outerBlocks(-1*angle_WingBottom, zStartInlet, zStartInletRamp, rInner_inlet, rOuterStart, rInner_inletRamp, side='RHS')

# WING UNDERSIDE
mate_wingBottom = outer_inletBlocks.undersideBlock()
mate_wingBottom_flatpaths = outer_inletBlocks.flatPaths
BLOCK_M3 = Block3D(mate_wingBottom, label="BLOCK-M3", 
                   nni=nRadialWall , nnj=nLongInlet, nnk=nNormalWingTop,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)
BLOCK_M3.bc_list[WEST] = FixedTBC(Twall=wall_temp)

# ENGINE 4
mate_engine4 = outer_inletBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], mate_wingBottom_flatpaths)
mate_engine4_flatpaths = outer_inletBlocks.flatPaths
BLOCK_M4 = Block3D(mate_engine4, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongInlet, nnk=nCircumEngine,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)
BLOCK_M4.bc_list[WEST] = FixedTBC(Twall=wall_temp)

# ENGINE 3
mate_engine3 = outer_inletBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], mate_engine4_flatpaths)
BLOCK_M5 = Block3D(mate_engine3, label="BLOCK-M5", 
                   nni=nRadialWall , nnj=nLongInlet, nnk=nCircumEngine,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[BOTTOM] = SlipWallBC(label='symmetry')
BLOCK_M5.bc_list[WEST] = FixedTBC(Twall=wall_temp)

# Wing
mate_wing = wing(mate_engine4_flatpaths)
BLOCK_W = Block3D(mate_wing, label="BLOCK-M5", 
                   nni=nRadialWall , nnj=nLongInlet, nnk=nWing,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_W.bc_list[EAST] = SupInBC(inflow)
BLOCK_W.bc_list[WEST] = FixedTBC(Twall=wall_temp)

### FORE DIAMONDS - INLET RAMP (1,2)
# Right wing, underside
foreDiamonds_RHwingBottom = rootFairing(L_inletRamp, zStartInletRamp, rInner_inletRamp, (wallAngle_inletRamp), phi_sweep, -1*angle_WingBottom, side='RHS', face='BOTTOM')
rootFairing_flatPaths = foreDiamonds_RHwingBottom.flatPaths
combustorTop3_wingBottom = foreDiamonds_RHwingBottom.southFace
combustorTop4_wingBottom = foreDiamonds_RHwingBottom.eastFace
BLOCK_1_2 = Block3D(foreDiamonds_RHwingBottom.volume, label="BLOCK-1_2", 
	    	  nni=nLongInletRamp, nnj=nLongInletRamp, nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 1, 0])
BLOCK_1_2.bc_list[TOP] = FixedTBC(Twall=wall_temp)
BLOCK_1_2.bc_list[WEST] = FixedTBC(Twall=wall_temp)

rOuterStart = outer_inletBlocks.rOuterEnd
rInnerEnd = foreDiamonds_RHwingBottom.rOuterEnd
outer_inletRampBlocks = outerBlocks(-1*angle_WingBottom, zStartInletRamp, zStartCombustor, rInner_inletRamp, rOuterStart, rInnerEnd, side='RHS')
outerInletRamp_wingBottom = outer_inletRampBlocks.undersideBlock()
wingBottom_flatpaths = outer_inletRampBlocks.flatPaths
BLOCK_M3 = Block3D(outerInletRamp_wingBottom, label="BLOCK-M3", 
                   nni=nRadialWall, nnj=nLongInletRamp, nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

# Engine 4
#print 'm8', foreDiamonds_RHwingBottom.surfacePath0, foreDiamonds_RHwingBottom.surfacePath_outer0, centreEnd20, centreEnd_outer20
foreDiamonds_RHengine4, combustorTop3_engine4, combustorTop4_engine4 = foreDiamonds_RHwingBottom.rotate(-1*engineSmile1[0], -1*engineSmile1[1], rootFairing_flatPaths, face='BOTTOM')
BLOCK_1_2 = Block3D(foreDiamonds_RHengine4, label="BLOCK-1_2", 
	    	  nni=nLongInletRamp, nnj=nLongInletRamp, nnk=nCircumEngine,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_1_2.bc_list[WEST] = FixedTBC(Twall=wall_temp)

outerInletRamp_engine4 = outer_inletRampBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], wingBottom_flatpaths)
engine4_flatpaths = outer_inletRampBlocks.flatPaths
BLOCK_M4 = Block3D(outerInletRamp_engine4, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongInletRamp, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)
 
# Engine 3
foreDiamonds_RHengine3, combustorTop3_engine3, combustorTop4_engine3 = foreDiamonds_RHwingBottom.rotate(-1*engineSmile2[0], -1*engineSmile2[1], rootFairing_flatPaths, face='BOTTOM')
BLOCK_1_2 = Block3D(foreDiamonds_RHengine3, label="BLOCK-1_2", 
	    	  nni=nLongInletRamp, nnj=nLongInletRamp, nnk=nCircumEngine,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_1_2.bc_list[BOTTOM] = SlipWallBC(label='symmetry')
BLOCK_1_2.bc_list[WEST] = FixedTBC(Twall=wall_temp)

outerInletRamp_engine3 = outer_inletRampBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], engine4_flatpaths)
BLOCK_M5 = Block3D(outerInletRamp_engine3, label="BLOCK-M5", 
                   nni=nRadialWall , nnj=nLongInletRamp, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

# Wing
outerInletRamp_wing = wing(engine4_flatpaths)
BLOCK_W = Block3D(outerInletRamp_wing, label="BLOCK-M5", 
                   nni=nRadialWall , nnj=nLongInletRamp, nnk=nWing,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_W.bc_list[EAST] = SupInBC(inflow)
BLOCK_W.bc_list[WEST] = FixedTBC(Twall=wall_temp)

# LE Rounding
#edgeRoundingBottom = edgeRounding('RHS')
#edgeRounding2 = edgeRoundingBottom.block_LE(r_leadingEdge, foreDiamonds_RHwingBottom.surfacePath0, foreDiamonds_RHwingBottom.surfacePath_outer0, centreEnd20, centreEnd_outer20, angle_WingBottom, theta_leadingEdge=90.*np.pi/180.)
#centreEnd2 = edgeRoundingBottom.centreEnd
#centreEnd_outer2 = edgeRoundingBottom.centreEnd_outer
#vol = edgeRounding21.clone()
#BLOCK_2 = Block3D(edgeRounding2, label="BLOCK-2 - LE FAIRING ROUNDING", 
#	    	  nni=nNormalEdgeRounding, nnj=nNormalWingTop, nnk=nLongInlet,
#	    	  fill_condition=initial)


### DIAMONDS (5-7)
diamonds_RHwingBottom = diamonds(-1*angle_WingBottom, diamondRotation, side='RHS', face='BOTTOM')

# Right wing, underside
BLOCK_5 = Block3D(diamonds_RHwingBottom.diamond5(), label="BLOCK-5", 
	    	  nni=nElevonNeighbour , nnj=nLongInletRamp, nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 1, 0])
westFace5_wingBottom = diamonds_RHwingBottom.westFace5
southFace5_wingBottom = diamonds_RHwingBottom.southFace5
BLOCK_5.bc_list[TOP] = FixedTBC(Twall=wall_temp)

BLOCK_6 = Block3D(diamonds_RHwingBottom.diamond6(), label="BLOCK-6", 
	    	  nni=nLongInletRamp, nnj=nElevonNeighbour , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 1, 0])
northFace6_wingBottom = diamonds_RHwingBottom.northFace6
eastFace6_wingBottom = diamonds_RHwingBottom.eastFace6
BLOCK_6.bc_list[TOP] = FixedTBC(Twall=wall_temp)

BLOCK_7 = Block3D(diamonds_RHwingBottom.diamond7(), label="BLOCK-7", 
	    	  nni=nLongInletRamp, nnj=nLongInletRamp, nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 1, 0])
southFace7_wingBottom = diamonds_RHwingBottom.southFace7
eastFace7_wingBottom = diamonds_RHwingBottom.eastFace7
BLOCK_7.bc_list[TOP] = FixedTBC(Twall=wall_temp)

# Engine 4
BLOCK_5 = Block3D(diamonds_RHwingBottom.rotate(-1*engineSmile1[0], -1*engineSmile1[1], diamonds_RHwingBottom.flatPaths5), label="BLOCK-5", 
	    	  nni=nElevonNeighbour , nnj=nLongInletRamp, nnk=nCircumEngine,
	    	  fill_condition=initial)
westFace5_RHengine4 = diamonds_RHwingBottom.westFace5
southFace5_RHengine4 = diamonds_RHwingBottom.southFace5


BLOCK_6 = Block3D(diamonds_RHwingBottom.rotate(-1*engineSmile1[0], -1*engineSmile1[1], diamonds_RHwingBottom.flatPaths6), label="BLOCK-6", 
	    	  nni=nLongInletRamp, nnj=nElevonNeighbour , nnk=nCircumEngine,
	    	  fill_condition=initial)
northFace6_RHengine4 = diamonds_RHwingBottom.northFace6
eastFace6_RHengine4 = diamonds_RHwingBottom.eastFace6


BLOCK_7 = Block3D(diamonds_RHwingBottom.rotate(-1*engineSmile1[0], -1*engineSmile1[1], diamonds_RHwingBottom.flatPaths7), label="BLOCK-7", 
	    	  nni=nLongInletRamp, nnj=nLongInletRamp, nnk=nCircumEngine,
	    	  fill_condition=initial)
southFace7_RHengine4 = diamonds_RHwingBottom.southFace7
eastFace7_RHengine4 = diamonds_RHwingBottom.eastFace7


# Engine 3
BLOCK_5 = Block3D(diamonds_RHwingBottom.rotate(-1*engineSmile2[0], -1*engineSmile2[1], diamonds_RHwingBottom.flatPaths5), label="BLOCK-5", 
	    	  nni=nElevonNeighbour , nnj=nLongInletRamp, nnk=nCircumEngine,
	    	  fill_condition=initial)
westFace5_RHengine3 = diamonds_RHwingBottom.westFace5
southFace5_RHengine3 = diamonds_RHwingBottom.southFace5
BLOCK_5.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

BLOCK_6 = Block3D(diamonds_RHwingBottom.rotate(-1*engineSmile2[0], -1*engineSmile2[1], diamonds_RHwingBottom.flatPaths6), label="BLOCK-6", 
	    	  nni=nLongInletRamp, nnj=nElevonNeighbour , nnk=nCircumEngine,
	    	  fill_condition=initial)
northFace6_RHengine3 = diamonds_RHwingBottom.northFace6
eastFace6_RHengine3 = diamonds_RHwingBottom.eastFace6
BLOCK_6.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

BLOCK_7 = Block3D(diamonds_RHwingBottom.rotate(-1*engineSmile2[0], -1*engineSmile2[1], diamonds_RHwingBottom.flatPaths7), label="BLOCK-7", 
	    	  nni=nLongInletRamp, nnj=nLongInletRamp , nnk=nCircumEngine,
	    	  fill_condition=initial)
southFace7_RHengine3 = diamonds_RHwingBottom.southFace7
eastFace7_RHengine3 = diamonds_RHwingBottom.eastFace7
BLOCK_7.bc_list[BOTTOM] = SlipWallBC(label='symmetry')


### FORE DIAMONDS - COMBUSTOR (3-4)
# Right wing, underside
topFace = combustorTop3_wingBottom #evalFace(BLOCK_1_VOL, 'S', rsReversal=1)
bottomFace = westFace5_wingBottom #diamonds_RHwingTop.westFace5 #evalFace(BLOCK_5_VOL, 'W')
pyfunction_fore_3 = lambda r,s,t: surf2surf(r, s, t, bottomFace, topFace)
BLOCK_3_VOL = PyFunctionVolume(pyfunction_fore_3)
BLOCK_3 = Block3D(BLOCK_3_VOL, label="BLOCK-3", 
	    	  nni=nLongInletRamp , nnj=nNormalWingTop , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 1, 1, 0, 0])
BLOCK_3.bc_list[WEST] = FixedTBC(Twall=wall_temp)
BLOCK_7.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)

topFace = combustorTop4_wingBottom #= evalFace(BLOCK_6_VOL, 'N')
bottomFace = northFace6_wingBottom #diamonds_RHwingTop.northFace6 #evalFace(BLOCK_2_VOL, 'S')
pyfunction_fore_4 = lambda r,s,t: surf2surf(r, s, t, bottomFace, topFace)
BLOCK_4_VOL = PyFunctionVolume(pyfunction_fore_4)
BLOCK_4 = Block3D(BLOCK_4_VOL, label="BLOCK-4", 
	    	  nni=nLongInletRamp , nnj=nNormalWingTop , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 1, 0, 0, 0])  
BLOCK_4.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)

rOuterStart = outer_inletRampBlocks.rOuterEnd
zDiamonds = diamonds_RHwingBottom.Nodetwoz
outer_foreDiamondsBlocks = outerBlocks(-1*angle_WingBottom, zStartCombustor, zDiamonds, rInnerEnd, rOuterStart, diamonds_RHwingBottom.rOuter, side='RHS')
outerforeDiamonds_wingBottom = outer_foreDiamondsBlocks.undersideBlock()
wingBottom_flatpaths = outer_foreDiamondsBlocks.flatPaths
BLOCK_M3 = Block3D(outerforeDiamonds_wingBottom, label="BLOCK-M3", 
                   nni=nRadialWall, nnj=int(nLongCombust/2), nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

# Engine 4
topFace = combustorTop3_engine4 #evalFace(BLOCK_1_VOL, 'S', rsReversal=1)
bottomFace = westFace5_RHengine4 #diamonds_RHwingTop.westFace5 #evalFace(BLOCK_5_VOL, 'W')
pyfunction_fore_3 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_3_VOL = PyFunctionVolume(pyfunction_fore_3)
BLOCK_3 = Block3D(BLOCK_3_VOL, label="BLOCK-3", 
	    	  nni=nLongInletRamp , nnj=nCircumEngine , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_3.bc_list[WEST] = FixedTBC(Twall=wall_temp)

topFace = combustorTop4_engine4 #= evalFace(BLOCK_6_VOL, 'N')
bottomFace = northFace6_RHengine4 #diamonds_RHwingTop.northFace6 #evalFace(BLOCK_2_VOL, 'S')
pyfunction_fore_4 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_4_VOL = PyFunctionVolume(pyfunction_fore_4)
BLOCK_4 = Block3D(BLOCK_4_VOL, label="BLOCK-4", 
	    	  nni=nLongInletRamp , nnj=nCircumEngine , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)  

outerforeDiamonds_engine4 = outer_foreDiamondsBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], wingBottom_flatpaths)
engine4_flatpaths = outer_foreDiamondsBlocks.flatPaths
BLOCK_M4 = Block3D(outerforeDiamonds_engine4, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=int(nLongCombust/2), nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)

# Engine 3
topFace = combustorTop3_engine3 #evalFace(BLOCK_1_VOL, 'S', rsReversal=1)
bottomFace = westFace5_RHengine3 #diamonds_RHwingTop.westFace5 #evalFace(BLOCK_5_VOL, 'W')
pyfunction_fore_3 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_3_VOL = PyFunctionVolume(pyfunction_fore_3)
BLOCK_3 = Block3D(BLOCK_3_VOL, label="BLOCK-3", 
	    	  nni=nLongInletRamp , nnj=nCircumEngine , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])  
BLOCK_3.bc_list[WEST] = FixedTBC(Twall=wall_temp)
BLOCK_3.bc_list[SOUTH] = SlipWallBC(label='symmetry')

topFace = combustorTop4_engine3 #= evalFace(BLOCK_6_VOL, 'N')
bottomFace = northFace6_RHengine3 #diamonds_RHwingTop.northFace6 #evalFace(BLOCK_2_VOL, 'S')
pyfunction_fore_4 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_4_VOL = PyFunctionVolume(pyfunction_fore_4)
BLOCK_4 = Block3D(BLOCK_4_VOL, label="BLOCK-4", 
	    	  nni=nLongInletRamp , nnj=nCircumEngine , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
BLOCK_4.bc_list[SOUTH] = SlipWallBC(label='symmetry')

outerforeDiamonds_engine3 = outer_foreDiamondsBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], engine4_flatpaths)
BLOCK_M5 = Block3D(outerforeDiamonds_engine3, label="BLOCK-M5", 
                   nni=nRadialWall, nnj=int(nLongCombust/2), nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

# Wing
outerforeDiamonds_wing = wing(engine4_flatpaths)
BLOCK_W = Block3D(outerforeDiamonds_wing, label="BLOCK-M5", 
                   nni=nRadialWall, nnj=int(nLongCombust/2), nnk=nWing,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 1, 0, 0])  
BLOCK_W.bc_list[EAST] = SupInBC(inflow)
BLOCK_W.bc_list[WEST] = FixedTBC(Twall=wall_temp)

### OUTLET - FORE (12-15)
rears_wingBottom = rearBlocks(0., -1*angle_WingBottom, 'RHS')

# Right wing, underside
rears_wingBottom_12 = rears_wingBottom.undersideBlock(L_outlet-L_elevon, (xStartElevon - rInner_outlet), zStartOutlet, rInner_outlet, wallAngle=wallAngle_outlet, block='WALL')
rears_flatPaths12 = rears_wingBottom.flatPaths
combustorBottom12_wingBottom = rears_wingBottom.northFace
BLOCK_12 = Block3D(rears_wingBottom_12, label="BLOCK-12", 
	    	  nni=nElevonNeighbour, nnj=nLongOutletFore, nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 1, 0])
BLOCK_12.bc_list[TOP] = FixedTBC(Twall=wall_temp)
BLOCK_12.bc_list[WEST] = FixedTBC(Twall=wall_temp)

rears_wingBottom_13 = rears_wingBottom.undersideBlock(L_outlet-L_elevon, W_elevon/2., zStartOutlet, xStartElevon, block='13')
rears_flatPaths13 = rears_wingBottom.flatPaths
combustorBottom13_wingBottom = rears_wingBottom.northFace
BLOCK_13 = Block3D(rears_wingBottom_13, label="BLOCK-13", 
	    	  nni=nLongInletRamp, nnj=nLongOutletFore , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 1, 0])
BLOCK_13.bc_list[TOP] = FixedTBC(Twall=wall_temp)
   
rears_wingBottom_14 = rears_wingBottom.undersideBlock(L_outlet-L_elevon, W_elevon/2., zStartOutlet, xStartElevon+W_elevon/2., block='14')
rears_flatPaths14 = rears_wingBottom.flatPaths
combustorBottom14_wingBottom = rears_wingBottom.northFace
BLOCK_14 = Block3D(rears_wingBottom_14, label="BLOCK-14", 
	    	  nni=nLongInletRamp, nnj=nLongOutletFore , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 1, 0])
BLOCK_14.bc_list[TOP] = FixedTBC(Twall=wall_temp)  

rears_wingBottom_15 = rears_wingBottom.undersideBlock(L_outlet-L_elevon, '', zStartOutlet, xStartElevon+W_elevon, block='EDGE')
rears_flatPaths15 = rears_wingBottom.flatPaths
combustorBottom15_wingBottom = rears_wingBottom.northFace
surfacePath15 = rears_wingBottom.pE
surfacePath_outer15 = rears_wingBottom.pEtop
BLOCK_15 = Block3D(rears_wingBottom_15, label="BLOCK-15", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletFore , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 1, 0])
BLOCK_15.bc_list[TOP] = FixedTBC(Twall=wall_temp)    


# Engine 4
rears_engine4_12 = rears_wingBottom.rotate(-1*engineSmile1[0], -1*engineSmile1[1], rears_flatPaths12)
combustorBottom12_engine4 = rears_wingBottom.northFace
BLOCK_12 = Block3D(rears_engine4_12, label="BLOCK-12", 
	    	  nni=nElevonNeighbour, nnj=nLongOutletFore , nnk=nCircumEngine,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_12.bc_list[WEST] = FixedTBC(Twall=wall_temp)

rears_engine4_13 = rears_wingBottom.rotate(-1*engineSmile1[0], -1*engineSmile1[1], rears_flatPaths13)
combustorBottom13_engine4 = rears_wingBottom.northFace
BLOCK_13 = Block3D(rears_engine4_13, label="BLOCK-13", 
	    	  nni=nLongInletRamp, nnj=nLongOutletFore , nnk=nCircumEngine,
	    	  fill_condition=initial)
   
rears_engine4_14 = rears_wingBottom.rotate(-1*engineSmile1[0], -1*engineSmile1[1], rears_flatPaths14)
combustorBottom14_engine4 = rears_wingBottom.northFace
BLOCK_14 = Block3D(rears_engine4_14, label="BLOCK-14", 
	    	  nni=nLongInletRamp, nnj=nLongOutletFore , nnk=nCircumEngine,
	    	  fill_condition=initial)      

rears_engine4_15 = rears_wingBottom.rotate(-1*engineSmile1[0], -1*engineSmile1[1], rears_flatPaths15)
combustorBottom15_engine4 = rears_wingBottom.northFace
BLOCK_15 = Block3D(rears_engine4_15, label="BLOCK-15", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletFore , nnk=nCircumEngine,
	    	  fill_condition=initial)        
    
# Engine 3
rears_engine3_12 = rears_wingBottom.rotate(-1*engineSmile2[0], -1*engineSmile2[1], rears_flatPaths12)
combustorBottom12_engine3 = rears_wingBottom.northFace
BLOCK_12 = Block3D(rears_engine3_12, label="BLOCK-12", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletFore , nnk=nCircumEngine,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_12.bc_list[WEST] = FixedTBC(Twall=wall_temp)
BLOCK_12.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

rears_engine3_13 = rears_wingBottom.rotate(-1*engineSmile2[0], -1*engineSmile2[1], rears_flatPaths13)
combustorBottom13_engine3= rears_wingBottom.northFace
BLOCK_13 = Block3D(rears_engine3_13, label="BLOCK-13", 
	    	  nni=nLongInletRamp, nnj=nLongOutletFore , nnk=nCircumEngine,
	    	  fill_condition=initial)
BLOCK_13.bc_list[BOTTOM] = SlipWallBC(label='symmetry')
   
rears_engine3_14 = rears_wingBottom.rotate(-1*engineSmile2[0], -1*engineSmile2[1], rears_flatPaths14)
combustorBottom14_engine3 = rears_wingBottom.northFace
BLOCK_14 = Block3D(rears_engine3_14, label="BLOCK-14", 
	    	  nni=nLongInletRamp, nnj=nLongOutletFore , nnk=nCircumEngine,
	    	  fill_condition=initial)
BLOCK_14.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

rears_engine3_15 = rears_wingBottom.rotate(-1*engineSmile2[0], -1*engineSmile2[1], rears_flatPaths15)
combustorBottom15_engine3 = rears_wingBottom.northFace
BLOCK_15 = Block3D(rears_engine3_15, label="BLOCK-15", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletFore , nnk=nCircumEngine,
	    	  fill_condition=initial)  
BLOCK_15.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

### AFT DIAMONDS - COMBUSTOR (8-11)
# Right wing, underside
topFace = southFace5_wingBottom
bottomFace = combustorBottom12_wingBottom
pyfunction_rear_8 = lambda r,s,t: surf2surf(r, s, t, bottomFace, topFace)
BLOCK_8_VOL = PyFunctionVolume(pyfunction_rear_8)
BLOCK_8 = Block3D(BLOCK_8_VOL, label="BLOCK-8", 
	    	  nni=nElevonNeighbour , nnj=nNormalWingTop , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 1, 1, 0, 0])
BLOCK_8.bc_list[SOUTH] = FixedTBC(Twall=wall_temp)    
BLOCK_8.bc_list[WEST] = FixedTBC(Twall=wall_temp) 

topFace = southFace7_wingBottom
bottomFace = combustorBottom13_wingBottom
pyfunction_rear_9 = lambda r,s,t: surf2surf(r, s, t, bottomFace, topFace)
BLOCK_9_VOL = PyFunctionVolume(pyfunction_rear_9)
BLOCK_9 = Block3D(BLOCK_9_VOL, label="BLOCK-9", 
	    	  nni=nLongInlet , nnj=nNormalWingTop , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 1, 0, 0, 0])
BLOCK_9.bc_list[SOUTH] = FixedTBC(Twall=wall_temp) 

topFace = eastFace7_wingBottom
bottomFace = combustorBottom14_wingBottom
pyfunction_rear_10 = lambda r,s,t: surf2surf(r, s, t, bottomFace, topFace)
BLOCK_10_VOL = PyFunctionVolume(pyfunction_rear_10)
BLOCK_10 = Block3D(BLOCK_10_VOL, label="BLOCK-10", 
	    	  nni=nLongInlet , nnj=nNormalWingTop , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 1, 0, 0, 0])
BLOCK_10.bc_list[SOUTH] = FixedTBC(Twall=wall_temp) 

topFace = eastFace6_wingBottom
bottomFace = combustorBottom15_wingBottom
pyfunction_rear_11 = lambda r,s,t: surf2surf(r, s, t, bottomFace, topFace)
BLOCK_11_VOL = PyFunctionVolume(pyfunction_rear_11)
BLOCK_11 = Block3D(BLOCK_11_VOL, label="BLOCK-11", 
	    	  nni=nElevonNeighbour , nnj=nNormalWingTop , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 1, 0, 0, 0])
BLOCK_11.bc_list[SOUTH] = FixedTBC(Twall=wall_temp) 

# Outer blocks, Fore diamonds 3-4
rOuterStart = outer_foreDiamondsBlocks.rOuterEnd
rInnerEnd = rears_wingBottom.rOuterStart
zDiamonds = diamonds_RHwingBottom.Nodetwoz
rInner_diamonds = diamonds_RHwingBottom.rOuter
outer_aftDiamondsBlocks = outerBlocks(-1*angle_WingBottom, zDiamonds, zStartOutlet, rInner_diamonds, rOuterStart, rInnerEnd, side='RHS')
outeraftDiamonds_wingBottom = outer_aftDiamondsBlocks.undersideBlock()
wingBottom_flatpaths = outer_aftDiamondsBlocks.flatPaths
BLOCK_M3 = Block3D(outeraftDiamonds_wingBottom, label="BLOCK-M3", 
                   nni=nRadialWall , nnj=int(nLongCombust/2), nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

outeraftDiamonds_engine4 = outer_aftDiamondsBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], wingBottom_flatpaths)
engine4_flatpaths = outer_aftDiamondsBlocks.flatPaths
BLOCK_M4 = Block3D(outeraftDiamonds_engine4, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=int(nLongCombust/2), nnk=nCircumEngine,
                   fill_condition=initial) 
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)

outeraftDiamonds_engine3 = outer_aftDiamondsBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], engine4_flatpaths)
BLOCK_M5 = Block3D(outeraftDiamonds_engine3, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=int(nLongCombust/2), nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

# Wing
outeraftDiamonds_wing = wing(engine4_flatpaths)
BLOCK_W = Block3D(outeraftDiamonds_wing, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=int(nLongCombust/2), nnk=nWing,
                   fill_condition=initial)
BLOCK_W.bc_list[EAST] = SupInBC(inflow)
BLOCK_W.bc_list[WEST] = FixedTBC(Twall=wall_temp)

rOuterStart = outer_aftDiamondsBlocks.rOuterEnd
rInnerEnd = rears_wingBottom.rOuterEnd
rInnerStart = rears_wingBottom.rOuterStart
outer_elevon1Blocks = outerBlocks(-1*angle_WingBottom, zStartOutlet, zStartElevon, rInnerStart, rOuterStart, rInnerEnd, side='RHS')
outer_elevon1_wingBottom = outer_elevon1Blocks.undersideBlock()
wingBottom_flatpaths = outer_elevon1Blocks.flatPaths
BLOCK_M3 = Block3D(outer_elevon1_wingBottom, label="BLOCK-M3", 
                   nni=nRadialWall , nnj=nLongOutletFore, nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

# Engine 4
topFace = southFace5_RHengine4
bottomFace = combustorBottom12_engine4
pyfunction_rear_8 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_8_VOL = PyFunctionVolume(pyfunction_rear_8)
BLOCK_8 = Block3D(BLOCK_8_VOL, label="BLOCK-8", 
	    	  nni=nElevonNeighbour , nnj=nCircumEngine , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_8.bc_list[WEST] = FixedTBC(Twall=wall_temp)

topFace = southFace7_RHengine4
bottomFace = combustorBottom13_engine4
pyfunction_rear_9 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_9_VOL = PyFunctionVolume(pyfunction_rear_9)
BLOCK_9 = Block3D(BLOCK_9_VOL, label="BLOCK-9", 
	    	  nni=nLongInletRamp, nnj=nCircumEngine , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)

topFace = eastFace7_RHengine4
bottomFace = combustorBottom14_engine4
pyfunction_rear_10 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_10_VOL = PyFunctionVolume(pyfunction_rear_10)
BLOCK_10 = Block3D(BLOCK_10_VOL, label="BLOCK-10", 
	    	  nni=nLongInletRamp, nnj=nCircumEngine , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)

topFace = eastFace6_RHengine4
bottomFace = combustorBottom15_engine4
pyfunction_rear_11 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_11_VOL = PyFunctionVolume(pyfunction_rear_11)
BLOCK_11 = Block3D(BLOCK_11_VOL, label="BLOCK-11", 
	    	  nni=nElevonNeighbour , nnj=nCircumEngine , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)

 
outerelevon1_engine4 = outer_elevon1Blocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], wingBottom_flatpaths)
engine4_flatpaths = outer_elevon1Blocks.flatPaths
BLOCK_M4 = Block3D(outerelevon1_engine4, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongOutletFore, nnk=nCircumEngine,
                   fill_condition=initial)  
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)

# Engine 3
topFace = southFace5_RHengine3
bottomFace = combustorBottom12_engine3
pyfunction_rear_8 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_8_VOL = PyFunctionVolume(pyfunction_rear_8)
BLOCK_8 = Block3D(BLOCK_8_VOL, label="BLOCK-8", 
	    	  nni=nElevonNeighbour , nnj=nCircumEngine , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_8.bc_list[WEST] = FixedTBC(Twall=wall_temp)
BLOCK_8.bc_list[SOUTH] = SlipWallBC(label='symmetry')

topFace = southFace7_RHengine3
bottomFace = combustorBottom13_engine3
pyfunction_rear_9 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_9_VOL = PyFunctionVolume(pyfunction_rear_9)
BLOCK_9 = Block3D(BLOCK_9_VOL, label="BLOCK-9", 
	    	  nni=nLongInletRamp, nnj=nCircumEngine , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
BLOCK_9.bc_list[SOUTH] = SlipWallBC(label='symmetry')

topFace = eastFace7_RHengine3
bottomFace = combustorBottom14_engine3
pyfunction_rear_10 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_10_VOL = PyFunctionVolume(pyfunction_rear_10)
BLOCK_10 = Block3D(BLOCK_10_VOL, label="BLOCK-10", 
	    	  nni=nLongInletRamp, nnj=nCircumEngine , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
BLOCK_10.bc_list[SOUTH] = SlipWallBC(label='symmetry')

topFace = eastFace6_RHengine3
bottomFace = combustorBottom15_engine3
pyfunction_rear_11 = lambda r,s,t: surf2surf(r, s, t, topFace, bottomFace)
BLOCK_11_VOL = PyFunctionVolume(pyfunction_rear_11)
BLOCK_11 = Block3D(BLOCK_11_VOL, label="BLOCK-11", 
	    	  nni=nElevonNeighbour , nnj=nCircumEngine , nnk=int(nLongCombust/2),
	    	  fill_condition=initial)
BLOCK_11.bc_list[SOUTH] = SlipWallBC(label='symmetry')

outerelevon1_engine3 = outer_elevon1Blocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], engine4_flatpaths)
BLOCK_M5 = Block3D(outerelevon1_engine3, label="BLOCK-M4", 
                   nni=nRadialWall, nnj=nLongOutletFore, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

# Wing
outerelevon1_wing = wing(engine4_flatpaths)
BLOCK_W = Block3D(outerelevon1_wing, label="BLOCK-M4", 
                   nni=nRadialWall, nnj=nLongOutletFore, nnk=nWing,
                   fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_W.bc_list[EAST] = SupInBC(inflow)
BLOCK_W.bc_list[WEST] = FixedTBC(Twall=wall_temp)

### OUTLET - AFT - ELEVONS (16-19)
# Right wing, underside
rears_wingBottom_19 = rears_wingBottom.undersideBlock(L_elevon, '', zStartElevon, xStartElevon+W_elevon, 'EDGE', TEflag=1)
surfacePath19 = rears_wingBottom.pE; surfacePath_outer19 = rears_wingBottom.pEtop; TE_surfacePath19 = rears_wingBottom.pS; TE_surfacePath_outer19 = rears_wingBottom.pStop
rears_flatPaths19 = rears_wingBottom.flatPaths
BLOCK_19 = Block3D(rears_wingBottom_19, label="BLOCK-19", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletElevon , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 1, 0])
BLOCK_19.bc_list[TOP] = FixedTBC(Twall=wall_temp)

## LE Rounding AFT OUTLET ROW
#BLOCK_25_VOL, centreEnd, centreEnd_outer = edgeRounding(r_leadingEdge, surfacePath19, surfacePath_outer19, centreEnd, centreEnd_outer, angle_WingTop, theta_leadingEdge=90.*np.pi/180, endFlag=1, side='RHS')
#BLOCK_25 = Block3D(BLOCK_25_VOL, label="BLOCK-25 - LE ELEVON ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nLongOutletElevon,
#	    	  fill_condition=initial) 
#
#BLOCK_26_VOL, centreEnd, centreEnd_outer = edgeRounding(r_leadingEdge, TE_surfacePath19, TE_surfacePath_outer19, centreEnd, centreEnd_outer, angle_WingTop, theta_leadingEdge=90.*np.pi/180, endFlag=1, TEFlag=1, side='RHS')
#BLOCK_26 = Block3D(BLOCK_26_VOL, label="BLOCK-26 - TE ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nElevonNeighbour,
#	    	  fill_condition=initial) 


rears_wingBottom_18 = rears_wingBottom.undersideBlock(L_elevon, W_elevon/2., zStartElevon, xStartElevon+W_elevon/2., '18', TEflag=1)
TE_surfacePath18 = rears_wingBottom.pS; TE_surfacePath_outer18 = rears_wingBottom.pStop; 
rears_flatPaths18 = rears_wingBottom.flatPaths
BLOCK_18 = Block3D(rears_wingBottom_18, label="BLOCK-18", 
	    	  nni=nLongInletRamp, nnj=nLongOutletElevon , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 1, 0])
BLOCK_18.bc_list[TOP] = FixedTBC(Twall=wall_temp)

#BLOCK_27_VOL, centreEnd, centreEnd_outer = edgeRounding(r_leadingEdge, TE_surfacePath18, TE_surfacePath_outer18, centreEnd, centreEnd_outer, angle_WingTop, theta_leadingEdge=90.*np.pi/180, endFlag=1, TEFlag=1, side='RHS')
#BLOCK_27 = Block3D(BLOCK_27_VOL, label="BLOCK-25 - TE ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nLongInlet,
#	    	  fill_condition=initial) 




rears_wingBottom_17 = rears_wingBottom.undersideBlock(L_elevon, W_elevon/2., zStartElevon, xStartElevon, '17', TEflag=1)
TE_surfacePath17 = rears_wingBottom.pS; TE_surfacePath_outer17 = rears_wingBottom.pStop
rears_flatPaths17 = rears_wingBottom.flatPaths
BLOCK_17 = Block3D(rears_wingBottom_17, label="BLOCK-17", 
	    	  nni=nLongInletRamp, nnj=nLongOutletElevon , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 0, 1, 0])
BLOCK_17.bc_list[TOP] = FixedTBC(Twall=wall_temp)
      
#BLOCK_28_VOL, centreEnd, centreEnd_outer = edgeRounding(r_leadingEdge, TE_surfacePath17, TE_surfacePath_outer17, centreEnd, centreEnd_outer, angle_WingTop, theta_leadingEdge=90.*np.pi/180, endFlag=1, TEFlag=1, side='RHS')
#BLOCK_28 = Block3D(BLOCK_28_VOL, label="BLOCK-25 - TE ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nLongInlet,
#	    	  fill_condition=initial) 

width = W_elevon/3. + (L_outlet-L_elevon)*np.tan(-1*wallAngle_outlet)
rears_wingBottom_16 = rears_wingBottom.undersideBlock(L_elevon, width, zStartElevon, (rInner_outlet+(L_outlet-L_elevon)*np.tan(wallAngle_outlet)), wallAngle=wallAngle_outlet, block='WALL', TEflag=1)
TE_surfacePath16 = rears_wingBottom.pS; TE_surfacePath_outer16 = rears_wingBottom.pStop
rears_flatPaths16 = rears_wingBottom.flatPaths
BLOCK_16 = Block3D(rears_wingBottom_16, label="BLOCK-16", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletElevon , nnk=nNormalWingTop,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 1, 0])
BLOCK_16.bc_list[TOP] = FixedTBC(Twall=wall_temp)
BLOCK_16.bc_list[WEST] = FixedTBC(Twall=wall_temp)



#BLOCK_29_VOL, centreEnd, centreEnd_outer = edgeRounding(r_leadingEdge, TE_surfacePath16, TE_surfacePath_outer16, centreEnd, centreEnd_outer, angle_WingTop, theta_leadingEdge=90.*np.pi/180, endFlag=1, TEFlag=1, side='RHS')
#BLOCK_29 = Block3D(BLOCK_29_VOL, label="BLOCK-25 - TE ROUNDING", 
#	    	  nni=nNormalEdgeRounding , nnj=nNormalWingTop , nnk=nElevonNeighbour,
#	    	  fill_condition=initial) 


rears_wingBottom.undersideBlock(L_elevon, '', zStartElevon, xStartElevon+W_elevon, 'EDGE', TEflag=1)
rOuterStart = outer_elevon1Blocks.rOuterEnd
rInnerEnd = rears_wingBottom.rOuterEnd
rInnerStart = rears_wingBottom.rOuterStart
outer_elevon2Blocks = outerBlocks(-1*angle_WingBottom, zStartElevon, zStartboatTail, rInnerStart, rOuterStart, rInnerEnd, side='RHS')
outer_elevon2_wingBottom = outer_elevon2Blocks.undersideBlock()
wingBottom_flatpaths = outer_elevon2Blocks.flatPaths
BLOCK_M3 = Block3D(outer_elevon2_wingBottom, label="BLOCK-M3", 
                   nni=nRadialWall , nnj=nLongOutletElevon, nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

# Engine 4
rears_engine4_16 = rears_wingBottom.rotate(-1*engineSmile1[0], -1*engineSmile1[1], rears_flatPaths16)
BLOCK_16 = Block3D(rears_engine4_16, label="BLOCK-16", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletElevon , nnk=nCircumEngine,
	    	  fill_condition=initial)
             #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_16.bc_list[WEST] = FixedTBC(Twall=wall_temp)

rears_engine4_17 = rears_wingBottom.rotate(-1*engineSmile1[0], -1*engineSmile1[1], rears_flatPaths17)
BLOCK_17 = Block3D(rears_engine4_17, label="BLOCK-17", 
	    	  nni=nLongInletRamp, nnj=nLongOutletElevon , nnk=nCircumEngine,
	    	  fill_condition=initial)
        
        
rears_engine4_18 = rears_wingBottom.rotate(-1*engineSmile1[0], -1*engineSmile1[1], rears_flatPaths18)
BLOCK_18 = Block3D(rears_engine4_18, label="BLOCK-18", 
	    	  nni=nLongInletRamp, nnj=nLongOutletElevon , nnk=nCircumEngine,
	    	  fill_condition=initial)  

rears_engine4_19 = rears_wingBottom.rotate(-1*engineSmile1[0], -1*engineSmile1[1], rears_flatPaths19)
BLOCK_19 = Block3D(rears_engine4_19, label="BLOCK-19", 
	    	  nni=nElevonNeighbour , nnj=nLongOutletElevon , nnk=nCircumEngine,
	    	  fill_condition=initial)

outerelevon2_engine4 = outer_elevon2Blocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], wingBottom_flatpaths)
engine4_flatpaths = outer_elevon2Blocks.flatPaths
BLOCK_M4 = Block3D(outerelevon2_engine4, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongOutletElevon, nnk=nCircumEngine,
                   fill_condition=initial) 
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)


# Engine 3
rears_engine3_16 = rears_wingBottom.rotate(-1*engineSmile2[0], -1*engineSmile2[1], rears_flatPaths16)
BLOCK_16 = Block3D(rears_engine3_16, label="BLOCK-16", 
                   nni=nElevonNeighbour , nnj=nLongOutletElevon , nnk=nCircumEngine,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_16.bc_list[WEST] = FixedTBC(Twall=wall_temp)
BLOCK_16.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

rears_engine3_17 = rears_wingBottom.rotate(-1*engineSmile2[0], -1*engineSmile2[1], rears_flatPaths17)
BLOCK_17 = Block3D(rears_engine3_17, label="BLOCK-17", 
                   nni=nLongInletRamp, nnj=nLongOutletElevon , nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_17.bc_list[BOTTOM] = SlipWallBC(label='symmetry')
     
rears_engine3_18 = rears_wingBottom.rotate(-1*engineSmile2[0], -1*engineSmile2[1], rears_flatPaths18)
BLOCK_18 = Block3D(rears_engine3_18, label="BLOCK-18", 
                   nni=nLongInletRamp, nnj=nLongOutletElevon , nnk=nCircumEngine,
                   fill_condition=initial)  
BLOCK_18.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

rears_engine3_19 = rears_wingBottom.rotate(-1*engineSmile2[0], -1*engineSmile2[1], rears_flatPaths19)
BLOCK_19 = Block3D(rears_engine3_19, label="BLOCK-19", 
                   nni=nElevonNeighbour , nnj=nLongOutletElevon , nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_19.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

outerelevon2_engine3 = outer_elevon2Blocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], engine4_flatpaths)
BLOCK_M5 = Block3D(outerelevon2_engine3, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongOutletElevon, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

# Wing
outerelevon2_wing = wing(engine4_flatpaths)
BLOCK_M5 = Block3D(outerelevon2_wing, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongOutletElevon, nnk=nWing,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 1, 1, 0])
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[WEST] = FixedTBC(Twall=wall_temp)

### BOAT TAIL (20-23)
# Right wing, underside
boatTailBlocks = BTandWake(0., -1*angle_WingBottom, 'RHS', face='BOTTOM')
boatTail_wingBottom_20 = boatTailBlocks.undersideBlock(L_boatTail, (xStartElevon - rInner_boatTail), zStartboatTail, rInner_boatTail, wallAngle=wallAngle_boatTail, block='WALL')
boatTail_flatPaths20 = boatTailBlocks.flatPaths
combustorBottom20_wingBottom = boatTailBlocks.northFace
BLOCK_20 = Block3D(boatTail_wingBottom_20, label="BLOCK-20", 
                   nni=nElevonNeighbour , nnj=nLongBoatTail , nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_20.bc_list[WEST] = FixedTBC(Twall=wall_temp)

boatTail_wingBottom_21 = boatTailBlocks.undersideBlock(L_boatTail, W_elevon/2., zStartboatTail, xStartElevon, block='21')
boatTail_flatPaths21 = boatTailBlocks.flatPaths
combustorBottom21_wingBottom = boatTailBlocks.northFace
BLOCK_21 = Block3D(boatTail_wingBottom_21, label="BLOCK-20", 
                   nni=nLongInletRamp, nnj=nLongBoatTail , nnk=nNormalWingTop,
                   fill_condition=initial)
   
boatTail_wingBottom_22 = boatTailBlocks.undersideBlock(L_boatTail, W_elevon/2., zStartboatTail, xStartElevon+W_elevon/2., block='22')
boatTail_flatPaths22 = boatTailBlocks.flatPaths
combustorBottom22_wingBottom = boatTailBlocks.northFace
BLOCK_22 = Block3D(boatTail_wingBottom_22, label="BLOCK-22", 
                   nni=nLongInletRamp, nnj=nLongBoatTail , nnk=nNormalWingTop,
                   fill_condition=initial)

boatTail_wingBottom_23 = boatTailBlocks.undersideBlock(L_boatTail, '', zStartboatTail, xStartElevon+W_elevon, meshAngle=5.*np.pi/180, block='EDGE')
boatTail_flatPaths23 = boatTailBlocks.flatPaths
combustorBottom23_wingBottom = boatTailBlocks.northFace
surfacePath23 = boatTailBlocks.pE
surfacePath_outer23 = boatTailBlocks.pEtop
BLOCK_23 = Block3D(boatTail_wingBottom_23, label="BLOCK-23", 
                   nni=nElevonNeighbour , nnj=nLongBoatTail , nnk=nNormalWingTop,
                   fill_condition=initial)

rOuterStart = outer_elevon2Blocks.rOuterEnd
rInnerEnd = boatTailBlocks.rOuterEnd
rInnerStart = boatTailBlocks.rOuterStart
outer_boatTailBlocks = outerBlocks(-1*angle_WingBottom, zStartboatTail, zStartWake, rInnerStart, rOuterStart, rInnerEnd, side='RHS')
outer_boatTail_wingBottom = outer_boatTailBlocks.undersideBlock()
wingBottom_flatpaths = outer_boatTailBlocks.flatPaths
BLOCK_M3 = Block3D(outer_boatTail_wingBottom, label="BLOCK-M3", 
                   nni=nRadialWall, nnj=nLongBoatTail, nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)

# Engine 4
boatTail_engine4_20 = boatTailBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], boatTail_flatPaths20, block='WALL')
BLOCK_20 = Block3D(boatTail_engine4_20, label="BLOCK-20", 
                   nni=nElevonNeighbour , nnj=nLongBoatTail , nnk=nCircumEngine,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_20.bc_list[WEST] = FixedTBC(Twall=wall_temp)

boatTail_engine4_21 = boatTailBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], boatTail_flatPaths21)
BLOCK_21 = Block3D(boatTail_engine4_21, label="BLOCK-21", 
                   nni=nLongInletRamp, nnj=nLongBoatTail , nnk=nCircumEngine,
                   fill_condition=initial)
        
        
boatTail_engine4_22 = boatTailBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], boatTail_flatPaths22)
BLOCK_22 = Block3D(boatTail_engine4_22, label="BLOCK-2", 
                   nni=nLongInletRamp, nnj=nLongBoatTail , nnk=nCircumEngine,
                   fill_condition=initial)  

boatTail_engine4_23 = boatTailBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], boatTail_flatPaths23)
BLOCK_23= Block3D(boatTail_engine4_23, label="BLOCK-23", 
                  nni=nElevonNeighbour , nnj=nLongBoatTail , nnk=nCircumEngine,
                  fill_condition=initial)

outerboatTail_engine4 = outer_boatTailBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], wingBottom_flatpaths)
engine4_flatpaths = outer_boatTailBlocks.flatPaths
BLOCK_M4 = Block3D(outerboatTail_engine4, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongBoatTail, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)

# Engine 3
boatTail_engine3_20 = boatTailBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], boatTail_flatPaths20, block='WALL')
BLOCK_20 = Block3D(boatTail_engine3_20, label="BLOCK-20", 
                   nni=nElevonNeighbour , nnj=nLongBoatTail , nnk=nCircumEngine,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 1, 0, 0])
BLOCK_20.bc_list[WEST] = FixedTBC(Twall=wall_temp)
BLOCK_20.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

boatTail_engine3_21 = boatTailBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], boatTail_flatPaths21)
BLOCK_21 = Block3D(boatTail_engine3_21, label="BLOCK-21", 
                   nni=nLongInletRamp, nnj=nLongBoatTail , nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_21.bc_list[BOTTOM] = SlipWallBC(label='symmetry')
    
boatTail_engine3_22 = boatTailBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], boatTail_flatPaths22)
BLOCK_22 = Block3D(boatTail_engine3_22, label="BLOCK-22", 
                   nni=nLongInletRamp, nnj=nLongBoatTail , nnk=nCircumEngine,
                   fill_condition=initial)  
BLOCK_22.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

boatTail_engine3_23 = boatTailBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], boatTail_flatPaths23)
BLOCK_23 = Block3D(boatTail_engine3_23, label="BLOCK-23", 
                   nni=nElevonNeighbour , nnj=nLongBoatTail , nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_23.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

outerboatTail_engine3 = outer_boatTailBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], engine4_flatpaths)
BLOCK_M5 = Block3D(outerboatTail_engine3, label="BLOCK-M4", 
                   nni=nRadialWall, nnj=nLongBoatTail, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[BOTTOM] = SlipWallBC(label='symmetry')

# Wing
boatTail_wing_20 = wing(boatTail_flatPaths20, block='WALL')
BLOCK_W20 = Block3D(boatTail_wing_20, label="BLOCK-20", 
                   nni=nElevonNeighbour , nnj=nLongBoatTail , nnk=nWing,
                   fill_condition=initial)
                   #xforce_list = [1, 0, 0, 1, 0, 0])
BLOCK_W20.bc_list[WEST] = FixedTBC(Twall=wall_temp)
BLOCK_W20.bc_list[NORTH] = FixedTBC(Twall=wall_temp)


boatTail_wing_21 = wing(boatTail_flatPaths21)
BLOCK_W21 = Block3D(boatTail_wing_21, label="BLOCK-21", 
                   nni=nLongInletRamp, nnj=nLongBoatTail , nnk=nWing,
                   fill_condition=initial)
                   #xforce_list = [1, 0, 0, 0, 0, 0])
BLOCK_W21.bc_list[NORTH] = FixedTBC(Twall=wall_temp)
        
boatTail_wing_22 = wing(boatTail_flatPaths22)
BLOCK_W22 = Block3D(boatTail_wing_22, label="BLOCK-22", 
                   nni=nLongInletRamp, nnj=nLongBoatTail , nnk=nWing,
                   fill_condition=initial)
                   #xforce_list = [1, 0, 0, 0, 0, 0])
BLOCK_W22.bc_list[NORTH] = FixedTBC(Twall=wall_temp)

boatTail_wing_23 = wing(boatTail_flatPaths23)
BLOCK_W23 = Block3D(boatTail_wing_23, label="BLOCK-23", 
                   nni=nElevonNeighbour , nnj=nLongBoatTail , nnk=nWing,
                   fill_condition=initial)
                   #xforce_list = [1, 0, 0, 0, 0, 0])
BLOCK_W23.bc_list[NORTH] = FixedTBC(Twall=wall_temp)

outerboatTail_wing = wing(engine4_flatpaths)
BLOCK_M5 = Block3D(outerboatTail_wing, label="BLOCK-M4", 
                   nni=nRadialWall, nnj=nLongBoatTail, nnk=nWing,
                   fill_condition=initial)
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)

### WAKE (24-27)
# Right wing, underside
wakeBlocks = BTandWake(0., -1*angle_WingBottom, 'RHS', face='BOTTOM')
wake_wingBottom_24 = wakeBlocks.undersideBlock(L_aftFarField, (xStartElevon - r_bt), zStartWake, r_bt, wallAngle=0., block='WALL')
wake_flatPaths24 = wakeBlocks.flatPaths
combustorBottom24_wingBottom = wakeBlocks.northFace
BLOCK_24 = Block3D(wake_wingBottom_24, label="BLOCK-24", 
	    	  nni=nElevonNeighbour, nnj=nLongWake, nnk=nNormalWingTop,
	    	  fill_condition=initial)
BLOCK_24.bc_list[SOUTH] = ExtrapolateOutBC()

wake_wingBottom_25 = wakeBlocks.undersideBlock(L_aftFarField, W_elevon/2., zStartWake, xStartElevon, block='25')
wake_flatPaths25 = wakeBlocks.flatPaths
combustorBottom25_wingBottom = wakeBlocks.northFace
BLOCK_25 = Block3D(wake_wingBottom_25, label="BLOCK-25", 
	    	  nni=nLongInletRamp, nnj=nLongWake, nnk=nNormalWingTop,
	    	  fill_condition=initial)
BLOCK_25.bc_list[SOUTH] = ExtrapolateOutBC()
   
wake_wingBottom_26 = wakeBlocks.undersideBlock(L_aftFarField, W_elevon/2., zStartWake, xStartElevon+W_elevon/2., block='26')
wake_flatPaths26 = wakeBlocks.flatPaths
combustorBottom26_wingBottom = wakeBlocks.northFace
BLOCK_26 = Block3D(wake_wingBottom_26, label="BLOCK-26", 
	    	  nni=nLongInletRamp, nnj=nLongWake, nnk=nNormalWingTop,
	    	  fill_condition=initial)      
BLOCK_26.bc_list[SOUTH] = ExtrapolateOutBC()

wake_wingBottom_27 = wakeBlocks.undersideBlock(L_aftFarField, '', zStartWake, xStartElevon+W_elevon, meshAngle=5.*np.pi/180, block='EDGE')
wake_flatPaths27 = wakeBlocks.flatPaths
combustorBottom27_wingBottom = wakeBlocks.northFace
surfacePath27 = wakeBlocks.pE
surfacePath_outer27 = boatTailBlocks.pEtop
BLOCK_27 = Block3D(wake_wingBottom_27, label="BLOCK-27", 
                   nni=nElevonNeighbour, nnj=nLongWake , nnk=nNormalWingTop,
                   fill_condition=initial)   
BLOCK_27.bc_list[SOUTH] = ExtrapolateOutBC()

rOuterStart = outer_boatTailBlocks.rOuterEnd
rInnerEnd = wakeBlocks.rOuterEnd
rInnerStart = wakeBlocks.rOuterStart
outer_wakeBlocks = outerBlocks(-1*angle_WingBottom, zStartWake, zStartWake-L_aftFarField, rInnerStart, rOuterStart, rInnerEnd, side='RHS')
outer_wake_wingBottom = outer_wakeBlocks.undersideBlock()
wingBottom_flatpaths = outer_wakeBlocks.flatPaths
BLOCK_M3 = Block3D(outer_wake_wingBottom, label="BLOCK-M3", 
                   nni=nRadialWall, nnj=nLongWake, nnk=nNormalWingTop,
                   fill_condition=initial)
BLOCK_M3.bc_list[EAST] = SupInBC(inflow)
BLOCK_M3.bc_list[SOUTH] = ExtrapolateOutBC()

# Engine 4
wake_engine4_24 = wakeBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], wake_flatPaths24, block='WALL')
BLOCK_24 = Block3D(wake_engine4_24, label="BLOCK-24", 
                   nni=nElevonNeighbour, nnj=nLongWake, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_24.bc_list[SOUTH] = ExtrapolateOutBC()

wake_engine4_25 = wakeBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], wake_flatPaths25)
BLOCK_25 = Block3D(wake_engine4_25, label="BLOCK-25", 
                   nni=nLongInletRamp, nnj=nLongWake, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_25.bc_list[SOUTH] = ExtrapolateOutBC()   
        
wake_engine4_26 = wakeBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], wake_flatPaths26)
BLOCK_26 = Block3D(wake_engine4_26, label="BLOCK-26", 
                   nni=nLongInletRamp, nnj=nLongWake, nnk=nCircumEngine,
                   fill_condition=initial)  
BLOCK_26.bc_list[SOUTH] = ExtrapolateOutBC()

wake_engine4_27 = wakeBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], wake_flatPaths27)
BLOCK_27 = Block3D(wake_engine4_27, label="BLOCK-27",
                  nni=nElevonNeighbour, nnj=nLongWake, nnk=nCircumEngine,
                  fill_condition=initial)
BLOCK_27.bc_list[SOUTH] = ExtrapolateOutBC()

outerwake_engine4 = outer_wakeBlocks.rotate(-1*engineSmile1[0], -1*engineSmile1[1], wingBottom_flatpaths)
engine4_flatpaths = outer_wakeBlocks.flatPaths
BLOCK_M4 = Block3D(outerwake_engine4, label="BLOCK-M4", 
                   nni=nRadialWall , nnj=nLongWake, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M4.bc_list[EAST] = SupInBC(inflow)
BLOCK_M4.bc_list[SOUTH] = ExtrapolateOutBC()

# Engine 3
wake_engine3_24 = wakeBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], wake_flatPaths24, block='WALL')
BLOCK_24 = Block3D(wake_engine3_24, label="BLOCK-24", 
                   nni=nElevonNeighbour, nnj=nLongWake, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_24.bc_list[BOTTOM] = SlipWallBC(label='symmetry')
BLOCK_24.bc_list[SOUTH] = ExtrapolateOutBC()

wake_engine3_25 = wakeBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], wake_flatPaths25)
BLOCK_25 = Block3D(wake_engine3_25, label="BLOCK-25", 
                   nni=nLongInletRamp, nnj=nLongWake, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_25.bc_list[BOTTOM] = SlipWallBC(label='symmetry')
BLOCK_25.bc_list[SOUTH] = ExtrapolateOutBC()

wake_engine3_26 = wakeBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], wake_flatPaths26)
BLOCK_26 = Block3D(wake_engine3_26, label="BLOCK-26", 
                   nni=nLongInletRamp, nnj=nLongWake, nnk=nCircumEngine,
                   fill_condition=initial)  
BLOCK_26.bc_list[BOTTOM] = SlipWallBC(label='symmetry')
BLOCK_26.bc_list[SOUTH] = ExtrapolateOutBC()

wake_engine3_27 = wakeBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], wake_flatPaths27)
BLOCK_27 = Block3D(wake_engine3_27, label="BLOCK-27", 
                   nni=nElevonNeighbour , nnj=nLongWake, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_27.bc_list[BOTTOM] = SlipWallBC(label='symmetry')
BLOCK_27.bc_list[SOUTH] = ExtrapolateOutBC()

outerwake_engine3 = outer_boatTailBlocks.rotate(-1*engineSmile2[0], -1*engineSmile2[1], engine4_flatpaths)
BLOCK_M5 = Block3D(outerwake_engine3, label="BLOCK-M4", 
                   nni=nRadialWall, nnj=nLongWake, nnk=nCircumEngine,
                   fill_condition=initial)
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[BOTTOM] = SlipWallBC(label='symmetry')
BLOCK_M5.bc_list[SOUTH] = ExtrapolateOutBC()

# Wing
wake_wing_24 = wing(wake_flatPaths24, block='WALL')
BLOCK_W24 = Block3D(wake_wing_24, label="BLOCK-24", 
                   nni=nElevonNeighbour, nnj=nLongWake, nnk=nWing,
                   fill_condition=initial)
BLOCK_W24.bc_list[SOUTH] = ExtrapolateOutBC()

wake_wing_25 = wing(wake_flatPaths25)
BLOCK_W25 = Block3D(wake_wing_25, label="BLOCK-25", 
                   nni=nLongInletRamp, nnj=nLongWake, nnk=nWing,
                   fill_condition=initial)
BLOCK_W25.bc_list[SOUTH] = ExtrapolateOutBC()

wake_wing_26 = wing(wake_flatPaths26)
BLOCK_W26 = Block3D(wake_wing_26, label="BLOCK-26", 
                   nni=nLongInletRamp, nnj=nLongWake, nnk=nWing,
                   fill_condition=initial)  
BLOCK_W26.bc_list[SOUTH] = ExtrapolateOutBC()

wake_wing_27 = wing(wake_flatPaths27)
BLOCK_W27 = Block3D(wake_wing_27, label="BLOCK-27", 
                   nni=nElevonNeighbour , nnj=nLongWake, nnk=nWing,
                   fill_condition=initial)
BLOCK_W27.bc_list[SOUTH] = ExtrapolateOutBC()

outerwake_wing = wing(engine4_flatpaths)
BLOCK_M5 = Block3D(outerwake_wing, label="BLOCK-M4", 
                   nni=nRadialWall, nnj=nLongWake, nnk=nWing,
                   fill_condition=initial)
BLOCK_M5.bc_list[EAST] = SupInBC(inflow)
BLOCK_M5.bc_list[SOUTH] = ExtrapolateOutBC()

### WAKE CENTRE (C0-C12)
wakeCentre = BTandWake(0., -1*angle_WingBottom, 'RHS', face='BOTTOM')
wakeCentre.centreBlocks()

BLOCK_C0 = Block3D(wakeCentre.centre0, label="BLOCK-C0", 
                   nni=nCircumTail, nnj=nRadialWake_Out, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C0.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C0.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_C0.bc_list[WEST] = SlipWallBC(label='symmetry')

BLOCK_C1 = Block3D(wakeCentre.centre1, label="BLOCK-C1", 
                   nni=nCircumFuselage , nnj=nRadialWake_Out, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C1.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C1.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
 
BLOCK_C2 = Block3D(wakeCentre.centre2, label="BLOCK-C2", 
                   nni=nRadialWake_Out , nnj=nNormalWingTop, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C2.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C2.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

BLOCK_C3 = Block3D(wakeCentre.centre3, label="BLOCK-C3", 
                   nni=nRadialWake_Out , nnj=nWing, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C3.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C3.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

BLOCK_C4 = Block3D(wakeCentre.centre4, label="BLOCK-C4", 
                   nni=nRadialWake_Out , nnj=nNormalWingTop, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C4.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C4.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

BLOCK_C5 = Block3D(wakeCentre.centre5, label="BLOCK-C5", 
                   nni=nRadialWake_Out, nnj=nCircumEngine, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C5.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C5.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

BLOCK_C6 = Block3D(wakeCentre.centre6, label="BLOCK-C6", 
                   nni=nCircumEngine , nnj=nRadialWake_Out, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C6.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C6.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_C6.bc_list[WEST] = SlipWallBC(label='symmetry')

BLOCK_C7 = Block3D(wakeCentre.centre7, label="BLOCK-C7", 
                   nni=nCircumEngine, nnj=nNormalWingTop, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C7.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C7.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_C7.bc_list[WEST] = SlipWallBC(label='symmetry')

BLOCK_C8 = Block3D(wakeCentre.centre8, label="BLOCK-C8", 
                   nni=nCircumEngine, nnj=nNormalWingTop, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C8.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C8.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
          
BLOCK_C9 = Block3D(wakeCentre.centre9, label="BLOCK-C9", 
                   nni=nCircumTail, nnj=nWing, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C9.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C9.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_C9.bc_list[WEST] = SlipWallBC(label='symmetry')

BLOCK_C10 = Block3D(wakeCentre.centre10, label="BLOCK-C7", 
                   nni=nCircumEngine, nnj=nWing, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C10.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C10.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)

BLOCK_C11 = Block3D(wakeCentre.centre11, label="BLOCK-C8", 
                   nni=nCircumEngine, nnj=nNormalWingTop, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C11.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C11.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
                  
BLOCK_C12 = Block3D(wakeCentre.centre12, label="BLOCK-C9", 
                   nni=nCircumTail, nnj=nNormalWingTop, nnk=nLongWake,
                   fill_condition=initial)
                   #xforce_list = [0, 0, 0, 0, 0, 1])
BLOCK_C12.bc_list[TOP] = ExtrapolateOutBC()
BLOCK_C12.bc_list[BOTTOM] = FixedTBC(Twall=wall_temp)
BLOCK_C12.bc_list[WEST] = SlipWallBC(label='symmetry')

#------------------------------  End of File ---------------------------------#

