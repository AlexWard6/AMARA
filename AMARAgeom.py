# AMARAgeom.py
# A. WARD 2016

"""
This is the master file and controls geometry flight conditions and cell discretisations. 
This is the file to give to e3prep.py/e3mpi.exe and also executes AMARAlibrary.py and 
AMARAblock.py.

It may be convenient to use this file as a record of the specific geometry for 
a given simulation.
    
Be careful of global variables.    
    
"""

import numpy as np
import math
import time

#---------- Global Data
gdata.title = "AMARA"; print gdata.title
gdata.dimensions = 3
gdata.axisymmetric_flag = 0
global phi_sweep; global diamondRotation

#------------------------------ Flow conditions -------------------------------#
# Select ideal air
select_gas_model(model='ideal gas', species=['air'])

g_gas = 1.4                      # Ideal Air, ratio of specific heats
R_gas = 287.058			 # Gas constant

# TEST 1 (alpha = 4.86)
flightConditions = [1943.04, 0., -1797.462*np.sin(4.86*np.pi/180), -1797.462*np.cos(4.86*np.pi/180), 223.32, 0., 1797.462, 299.577]
# TEST 2 (alpha = 4.2)
#flightConditions = [1263.70, 0., -2110.29*np.sin(6.*np.pi/180), -2110.29*np.cos(6.*np.pi/180), 226.15, 0., 2110.29, 301.47]
# TEST 3 (alpha = 3.93)
#flightConditions = [793.85, 0., -3650.88*np.sin(6.*np.pi/180), -3650.88*np.cos(6.*np.pi/180), 230.33, 0., 3650.88, 304.24]
# TEST 4 (alpha = 5.)
#flightConditions = [26436.27, 0., -598.92*np.sin(6.*np.pi/180), -598.92*np.cos(6.*np.pi/180), 223.15, 0., 598.92, 298.455]



# Freestream
p = flightConditions[0]     # pressure [Pa]
u = flightConditions[1]     # x velocity [m/s]
v = flightConditions[2]     # y velocity [m/s]
w = flightConditions[3]     # y velocity [m/s]
T = flightConditions[4]     # temperature [K]
alpha = flightConditions[5] # AoA [rads]
vel = flightConditions[6]   # freestream velocity [m/s]
a = flightConditions[7]     # sound speed [m/s]
M = vel/a                   # Mach number

# Initialise freestream with above parameters
inflow  = FlowCondition(p=p, u=u, v=v, w=w, T=T)

# Define initial flow at low temperature and pressure
#initial = FlowCondition(p=p*0.1, u=0., v=0. , T=T*0.1)
initial = inflow


newSim = 1
# If we are starting a set of new flow conditions:
if newSim:
    # Initialise freestream with above parameters
    inflow  = FlowCondition(p=p, u=u, v=v, w=w, T=T)
else:
    print 'uh oh'    
    old_flow = ExistingSolution(rootName, solutionworkDir, nblock, tindx,
                                dimensions=3, assume_same_grid=1, zipFiles=1,
                                add_velocity=Vector(0., 0., 0.))

# Temperature of vehicle
wall_temp = 300.


#--------------------------------- GEOMETRY ----------------------------------#

#---------- Parametric Variables
# Forebody
r_fb = 1.05              # Radius of centrebody/forebody base/boat tail base [m]
phi_fb = 5./180*np.pi    # Half angle of forebody ramp [rads]
r_nose = 0.05		   # Radius of spherical nose [m] # (0.1)
AR = 1.0                 # Aspect ratio 1 <= AR <= 2 (1)

C_a_long = 0.0          # Longitudinal camber
C_a_trans = 0.0        # Transverse camber

# Wings
halfspan = 4.675		       # Wing half span [m] (4.6745)

r_leadingEdge = 0.03
theta_leadingEdge = 85.*np.pi/180
r_root = 0.05

# Engine
L_engine = 6. #(6)
inletFraction = 0.100
inletRampFraction = 0.151
combustorFraction = 0.574
outletFraction = 0.175
# Ratios estimated from Jazra (2010) and Wise (2014)
assert inletFraction+inletRampFraction+combustorFraction+outletFraction == 1., "Sum of engine length fractions is not 1."

# These are the angles the engine parts make wrt to horizontal
wallAngle_inlet = 35.*np.pi/180
wallAngle_inletRamp = 5.*np.pi/180
wallAngle_combustor = 0.*np.pi/180

# Boat Tail
wallAngle_boatTail = -14.*np.pi/180   # Boat tail angle  
L_boatTail = 2.*r_fb                   # Length boat tail r_fb < L < 3*r_fb

#---------- Define mesh
phiTip = phi_fb*1.  		# angle phi for N nodes, define square nose block
phiRing = phiTip*2.5       	# angle phi for M nodes, define ring
phiJoin = 90./180*np.pi - phi_fb # angle phi for M nodes, define join

t_nose = 1.5     		# Thickness of mesh in Radial direction [m], nose
t_noseFF = 1.     		# Thickness of mesh in Radial direction [m], nose
t_join = t_nose     		# Thickness of mesh in Radial direction [m], join
t_ramp = t_nose*math.sin(phiJoin)	# Thickness of mesh in radial direction [m], ramp

phi_mesh = 15./180*np.pi     # Angle forebody mesh makes wrt z axis
phi_meshWing = phi_mesh/5.   # Angle wing mesh makes wrt z axis

L_aftFarField = 3.

diamondRotation = 0.*np.pi/180.


phiJoinOuter = 90./180*np.pi - phi_mesh # angle phi for M nodes, define join

angle_Tip0 = 0. *np.pi/180
angle_Tip1 = 90. *np.pi/180
angle_Tip2 = 180. *np.pi/180
angle_Tip3 = -90. *np.pi/180

angle_0 = 0. *np.pi/180
angle_WingTop = 7.5*np.pi/180
angle_WingBottom = angle_WingTop
angle_Tail = 60. *np.pi/180
angle_90 = 90. *np.pi/180
angle_180 = 180. *np.pi/180
angle_270 = 270. *np.pi/180

angle_Ring0 = np.arctan(0.5)
angle_Join0 = 30. *np.pi/180

inletSmile = 30. *np.pi/180
enginesSmile = angle_180 - 2.* angle_WingBottom
inletSmile = enginesSmile/4.

angle_Engine1 = -1* (angle_WingBottom + inletSmile/2)
angle_Engine2 = angle_Engine1 - inletSmile
angle_Engine3 = angle_Engine2 - inletSmile
angle_Engine4 = angle_Engine3 - inletSmile

engineSmile1 = (angle_Engine1-inletSmile/2, angle_Engine1+inletSmile/2)
engineSmile2 = (angle_Engine2-inletSmile/2, angle_Engine2+inletSmile/2)
engineSmile3 = (angle_Engine3-inletSmile/2, angle_Engine3+inletSmile/2)
engineSmile4 = (angle_Engine4-inletSmile/2, angle_Engine4+inletSmile/2)



#---------- Dependent Stuff
# Forebody 
r_join = r_nose*np.cos(phi_fb) 			      # 2D radius at join-ramp
L_rampFull = (r_fb - r_join)/np.tan(phi_fb)		# Length of ramp [m]
L_fb = L_rampFull + r_nose*(1 - np.sin(phi_fb))	# Length of forebody [m]


# Inlet
L_inlet = inletFraction*L_engine  # Length of the inlet [m]
L_inletRamp = inletRampFraction*L_engine # length of the inlet ramp [m]
rInner_inlet = (L_rampFull - L_inlet - L_inletRamp)*np.tan(phi_fb) + r_join
meshAngle_inlet = phi_mesh
rOuterStartInlet = rInner_inlet + t_ramp + np.tan(phi_mesh)*(L_rampFull - L_inlet - L_inletRamp)
zStartInlet = -1.*L_rampFull + L_inlet + L_inletRamp + np.cos(phiJoin)*r_nose

# Inlet ramp
rInner_inletRamp = rInner_inlet + L_inlet*np.tan(wallAngle_inlet)
meshAngle_inletRamp = phi_mesh
rOuterStartInletRamp = rOuterStartInlet + L_inlet*np.tan(meshAngle_inlet)
zStartInletRamp = zStartInlet - L_inlet

dr_ramp = (L_rampFull - L_inlet - L_inletRamp)*np.tan(phi_fb)
L_ramp = L_rampFull - L_inlet - L_inletRamp
xStartInletRamp = L_ramp*np.tan(phi_fb)
xStartInletRamp = L_ramp*np.tan(phi_fb)

# Combustor
rInner_combustor = rInner_inletRamp + L_inletRamp*np.tan(wallAngle_inletRamp)
L_combustor = combustorFraction*L_engine
meshAngle_combustor = phi_mesh
rOuterStartCombustor = rOuterStartInletRamp + L_inletRamp*np.tan(meshAngle_inletRamp)
zStartCombustor = zStartInletRamp - L_inletRamp

# Outlet
rInner_outlet = rInner_combustor + L_combustor*np.tan(wallAngle_combustor)
L_outlet = outletFraction*L_engine
wallAngle_outlet = -1*np.arctan((rInner_outlet - r_fb)/L_outlet)
meshAngle_outlet = phi_mesh
rOuterStartOutlet = rOuterStartCombustor + L_combustor* np.tan(meshAngle_combustor)
zStartOutlet = zStartCombustor - L_combustor

# Boat tail
rInner_boatTail = r_fb
meshAngle_boatTail = phi_mesh
rOuterStartboatTail = rOuterStartOutlet + L_outlet* np.tan(meshAngle_outlet)
zStartboatTail = zStartOutlet - L_outlet
r_bt = (r_fb-abs(np.tan(wallAngle_boatTail)*L_boatTail)) # radius of boat tail (truncated face)

# Wake
zStartWake = zStartboatTail - L_boatTail
r_wakeCentre = r_bt/2.


# Wings
L_wing = L_inlet + L_inletRamp + L_combustor

A_wing = 0.5*L_wing*(halfspan-r_fb)       # Area of the elevons
L_elevon = 0.15*halfspan
W_elevon = 0.2*halfspan
A_elevon = 0.15*A_wing

theta_trailingEdge = theta_leadingEdge

zStartElevon = zStartboatTail + L_elevon
xStartElevon = rInner_outlet + W_elevon/3.

xStartInlet = (L_ramp+L_inletRamp)*np.tan(phi_fb)
xStartInletRamp = L_ramp*np.tan(phi_fb)

phi_sweep = np.pi*90./180. - np.arctan((halfspan - rInner_inletRamp)/(L_engine-L_inlet))	# Wing sweep angle [rads]
if phi_sweep > 80.*np.pi/180:
    print 'WARNING: High wing sweep, check grid quality and think about reducing cfl.'
    flag=1

theta_trailingEdge = theta_leadingEdge

zStartElevon = zStartboatTail + L_elevon
xStartElevon = rInner_outlet + W_elevon/3.

L_vehicle = r_nose + L_rampFull + L_combustor + L_boatTail

#-------------------------- CELLS, CLUSTERING  -------------------------------#
#---------- Define Cells
factor = 1.5
# radial cells out from wall
nRadialWall = int(round(14*factor)) #13
# radial cells out from wall
nRadialFF = int(round(2*factor))
# Cells normal to wing top surface
# Cells normal to wing top surface
nNormalWingTop = int(round(2*factor))
# Cells normal to wing bottom surface
nNormalWingBottom = int(round(2*factor))
# Cells across thickness of wing
nWing = int(round(2))
# Circumferential cells over Engines
nCircumEngine = int(round(10*factor))
# Radial cells in wake blocks 0-6
nRadialWake_Out = int(round(5*factor))
# Radial cells in Nose blocks 0-6
nRadialNose_Out = int(round(4*factor))
# Circumferential cells over tail
nCircumTail = nCircumEngine
# Circumferential cells over fuselage
nCircumFuselage = nCircumEngine
# longitudinal cells in join
nLongJoin = int(round(nNormalWingBottom*(phiJoin - phiTip)/phiTip))
# longitudinal cells in ramp
nLongRamp = int(2.2*round(nLongJoin * L_ramp/np.cos(phi_fb)/((2.*np.pi*t_join)*(phiJoin - phiRing)/2.*np.pi)))
# longitudinal cells inlet
nLongInlet = max([int(round(nLongRamp * L_inlet/L_ramp)), 2])
# longitudinal cells inlet ramp
nLongInletRamp = int(round(nLongInlet * L_inletRamp/L_inlet))
nLongTemp = (nLongInlet + nLongInletRamp)/2.
nLongInlet = max([int(round(nLongTemp)), 2])
nLongInletRamp = nLongInlet
# longitudinal cells combustor
nLongCombust = int(round(nLongInlet * L_combustor/L_inlet))
# longitudinal cells outlet
nLongOutlet = int(round(nLongInlet * L_outlet/L_inlet))
# longitudinal cells outlet, fore of elevons
nLongOutletFore = max([int(round(nLongOutlet * (L_outlet-L_elevon)/L_outlet)), 2])
# longitudinal cells outlet, over elevons
nLongOutletElevon = int(round(nLongOutlet * (L_elevon)/L_outlet))
# longitudinal cells boat tail
nLongBoatTail = int(round(nLongInlet * L_boatTail/L_inlet))
# longitudinal cells wake
nLongWake = int(round(nLongInlet * L_aftFarField/L_inlet))
# Transverse cells either side of the elevons
nElevonNeighbour = max([int(round(2.0*nLongInlet * (xStartElevon-r_fb)/W_elevon)), 2])



#---------- Clustering
# Normal to the wall
cf_wallNormal = HypertanClusterFunction(1, 0.5)
cf_wallNormal = None

# Call the functions
execfile('AMARAlibrary.py', globals())

# Build the blocks
execfile('AMARAblock.py', globals())

# Connect blocks and set adjacent BCs
identify_block_connections()

## Write log file
#f = open('amara.log', 'w+')
#f.write(time.strftime("Run %d-%m_%H:%M:%S.py"))
#f.write('Geometry')
#f.write('r_nose')
#f.write('phi_fb')
#f.write('r_fb')
#f.write('halfspan')
#
#f.write('Discretisation factor = %d')

print 'phi_sweep', phi_sweep*180./np.pi
print 'wallAngle_inletRamp', wallAngle_inletRamp*180./np.pi
print 'Block angle', 90.- (phi_sweep*180./np.pi) - wallAngle_inletRamp*180./np.pi 

### History cells
z = matingFace.eval(0, 0).z - L_ramp
HistoryLocation(0., 0., r_nose)	# Stagnation point
HistoryLocation(0., r_fb, z)	# Top Shoulder
HistoryLocation(-1.*r_fb, 0., z)	      # Side Shoulder
HistoryLocation(0., -1.*r_fb, z)	# Bottom Shoulder
HistoryLocation(0., -1.*r_fb, zStartboatTail)	# Start of Boat Tail
HistoryLocation(0., 1.*r_fb, zStartboatTail)	# Start of Boat Tail
HistoryLocation(0., 0., zStartWake)	# Wake

#-------------------------------- GLOBAL DATA --------------------------------#

gdata.viscous_flag = 0			# 0 inviscid
gdata.cfl = 0.2				# Starting CFL number
gdata.stringent_cfl=1
gdata.cfl_count=8
gdata.flux_calc = ADAPTIVE		# 
#gdata.gasdynamic_update_scheme = "classic-rk3"
flowtime = L_vehicle/vel		# Characteristic flow time
gdata.max_time = 2.*flowtime		# 0.2 ms ~= 5 flowtimes
gdata.max_step = 1e10			# 
gdata.dt = 1.0e-8			#
gdata.dt_plot = gdata.max_time/20.	# Write 100 flow solutions
gdata.dt_history = 1e-5	

#------------------------------  End of File ---------------------------------#
