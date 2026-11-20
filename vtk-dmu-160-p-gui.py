#!/usr/bin/env python3
import hal
import math
import sys
from vtk_vismach import *

# get current working directory
cwd =  os.getcwd()
# define path for machine model files
path_stl = cwd + "/vismach-stl-files/"
# define path for tool model files
path_tool_stl = cwd + "/vismach-stl-files/tools/"

try: # Expect files in working directory
    # NOTE to run this file as standalone python script absolute paths might have to be used to find the stl files
    # create the work piece from file
    work_piece = ReadPolyData('work_piece_1.stl', path_stl).set_group('work')
    EGO_BC = ReadPolyData('160Pbase.stl', path_stl).set_group('machine')
    EGO_B = ReadPolyData('160PB.stl', path_stl).set_group('machine')
    EGO_C = ReadPolyData('160PC.stl', path_stl)
    EGO_X = ReadPolyData('160PX.stl', path_stl).set_group('machine')
    EGO_Y = ReadPolyData('160PY.stl', path_stl).set_group('machine')
    EGO_Z = ReadPolyData('160PZ.stl', path_stl).set_group('machine')
except Exception as detail:
    print(detail)
    raise SystemExit('Vismach requires 3d files in working directory')

#for setting in sys.argv[1:]: exec(setting)
options = sys.argv[1:]

c = hal.component('vtk-dmu-160-p-gui')
# nutation-angle
c.newpin("nutation-angle", hal.Type.REAL, hal.HAL_IN)
# geometric offsets in the spindle-rotary-assembly
c.newpin("pivot-y", hal.Type.REAL, hal.HAL_IN)
c.newpin("pivot-z", hal.Type.REAL, hal.HAL_IN)
c.newpin("offset-x", hal.Type.REAL, hal.HAL_IN)
c.newpin("offset-z", hal.Type.REAL, hal.HAL_IN)
# rot-point offsets distances from pivot point ( spindle AB)
# to rotation axis ( table C)
c.newpin("rot-axis-x", hal.Type.REAL, hal.HAL_IN)
c.newpin("rot-axis-y", hal.Type.REAL, hal.HAL_IN)
# work piece show/hide
c.newpin("hide-work-piece", hal.HAL_BOOL, hal.HAL_IN)
# spindle body show/hide
c.newpin("hide-spindle-body", hal.HAL_BOOL, hal.HAL_IN)
# scale coordinate system indicators
c.newpin("scale-coords", hal.Type.REAL, hal.HAL_IN)

c.newpin("g68-work-rotation", hal.Type.REAL, hal.HAL_OUT)
c.newpin("g68-work-origin-x", hal.Type.REAL, hal.HAL_OUT)
c.newpin("g68-work-origin-y", hal.Type.REAL, hal.HAL_OUT)
c.newpin("g68-work-origin-z", hal.Type.REAL, hal.HAL_OUT)
c.ready()

# Machine zero as measured from center surface of the rotary c table
machine_zero_x = -1000
machine_zero_y =  1000
machine_zero_z =  1000
# how the axis letters map to the joints
pin_axis_x = 'joint.0.pos-fb'
pin_axis_y = 'joint.1.pos-fb'
pin_axis_z = 'joint.2.pos-fb'
pin_rotary_b = 'joint.3.pos-fb'
pin_rotary_c = 'joint.4.pos-fb'

# function
work_rotation_old = None
g5x_rotation_old = None
g92_g68_offset_old = None

def myfunction():
    def mult_matx(A, B):
        rowsA = len(A)
        colsB = len(B[0])
        result = [[0] * colsB for i in range(rowsA)]
        for i in range(rowsA):
            for j in range(colsB):
                for k in range(len(B)):
                    result[i][j] += A[i][k] * B[k][j]
        return result

    global work_rotation_old
    s.poll()
    work_rotation = [[ s.g68_kins_rotation[0], s.g68_kins_rotation[1], s.g68_kins_rotation[2] ],
                     [ s.g68_kins_rotation[3], s.g68_kins_rotation[4], s.g68_kins_rotation[5] ],
                     [ s.g68_kins_rotation[6], s.g68_kins_rotation[7], s.g68_kins_rotation[8] ]]
    if work_rotation != work_rotation_old:
        # calculate rotation around Z from the rotation matrix
        c["g68-work-rotation"] = -degrees(atan2(s.g68_kins_rotation[3], s.g68_kins_rotation[0]))
    work_rotation_old = work_rotation
    global g5x_rotation_old, g92_g68_offset_old, g92_ijk_offset
    g92_g68_offset = [[ s.g92_offset[0] + s.g68_offset[0],
                        s.g92_offset[1] + s.g68_offset[1],
                        s.g92_offset[2] + s.g68_offset[2] ]]
    g5x_rotation = (s.rotation_yz, s.rotation_zx, s.rotation_xy)
    if (g5x_rotation != g5x_rotation_old) or (g92_g68_offset != g92_g68_offset_old):
        # rotate g92 to match the G10 L2 I J K(R) offsets values
        rx = -radians(s.rotation_yz)
        Rx = [[1,       0,        0],
            [0, cos(rx), -sin(rx)],
            [0, sin(rx),  cos(rx)]]
        g92_i = mult_matx(g92_g68_offset, Rx)

        ry = -radians(s.rotation_zx)
        Ry = [[ cos(ry), 0, sin(ry)],
            [       0, 1,       0],
            [-sin(ry), 0, cos(ry)]]
        g92_ij = mult_matx(g92_i, Ry)

        rz = -radians(s.rotation_xy)
        Rz = [[cos(rz), -sin(rz), 0],
            [sin(rz),  cos(rz), 0],
            [      0,        0, 1]]
        g92_ijk = mult_matx(g92_ij, Rz)
        g92_ijk_offset = (g92_ijk[0][0], g92_ijk[0][1], g92_ijk[0][2])
    g5x_rotation_old = g5x_rotation
    g92_g68_offset_old = g92_g68_offset
    c["g68-work-origin-x"] = g92_ijk_offset[0] + s.g5x_offset[0] - c["rot-axis-x"]
    c["g68-work-origin-y"] = g92_ijk_offset[1] + s.g5x_offset[1] - c["rot-axis-y"]
    c["g68-work-origin-z"] = g92_ijk_offset[2] + s.g5x_offset[2] + machine_zero_z
# /function

# Create an indicator for the machine reference coordinates
machine_coords = Axes(c,'200*{scale-coords}')
machine_coords = Translate([machine_coords], machine_zero_x, 0, machine_zero_z)

# start toolside
# Create the tooltip tracker (tool control point)
tooltip = Capture('tool')
# create indicator for tool coordinates
tool_coords = Axes(c,'50*{scale-coords}')
tool_coords = Collection([tool_coords,
                          tooltip
                        ])
# create a visual indicator for the position of the control point
tool_coords = Translate([tool_coords], hal,
                           '-{motion.tooloffset.x}',
                           '-{motion.tooloffset.y}',
                           '-{motion.tooloffset.z}')
# create an indicator for the tool-offsets when spindle-body is hidden
ind_tool_offset_x = CylinderX(hal,'-{motion.tooloffset.x}',2)
ind_tool_offset_x = Translate([ind_tool_offset_x],hal,0,0,'-{motion.tooloffset.z}')
ind_tool_offset_y = CylinderY(hal,'-{motion.tooloffset.y}',2)
ind_tool_offset_y = Translate([ind_tool_offset_y],hal,
                                                    '-{motion.tooloffset.x}',
                                                    0,
                                                    '-{motion.tooloffset.z}')
ind_tool_offset_z = CylinderZ(hal,'-{motion.tooloffset.z}',2)
ind_tool_offset = Collection([ind_tool_offset_x,
                              ind_tool_offset_y,
                              ind_tool_offset_z])
ind_tool_offset = Color([ind_tool_offset],1,0,1,1)
tool_cylinder = Collection([
                CylinderZ(hal,'-{motion.tooloffset.z}', '0.5*{halui.tool.diameter}'),
                # this indicates the spindle nose when no tool-offset is active
                CylinderZ(-0.1, 0),
                ]).set_group('tool')
tool_cylinder = Color([tool_cylinder],1,0,1,1)
tool = Collection([
       tool_coords,
       ind_tool_offset,
       tool_cylinder
       ])
tool = Translate([tool],c,0,'-{pivot-y}','-{pivot-z}')
# create spindle head
EGO_B = Color([EGO_B],0.7,0.7,0,1)
# rotate the nutation joint to the nutation angle, 90° should have the nutation axis in the horizontal plane
EGO_B = Rotate([EGO_B],-90,-1,0,0)
EGO_B = Rotate([EGO_B],c,'nutation-angle',-1,0,0)
spindle_assembly = Collection([
                   tool,
                   EGO_B
                   ])
# create HAL-link for b-axis rotational joint'
spindle_assembly = Rotate([spindle_assembly],
                hal,
                'joint.3.pos-fb',
                0,
                'sin(radians( {vtk-dmu-160-p-gui.nutation-angle} ))',
                'cos(radians( {vtk-dmu-160-p-gui.nutation-angle} ))'
                )
# Z carriage
EGO_Z = Color([EGO_Z],1,0.5,0,1)
EGO_Z = Translate([EGO_Z], 0, 0, -759.5)
spindle_z = Collection([
             spindle_assembly,
             EGO_Z
             ])
# move by the values set for the pivot lengths so the vismach origin is in the center of the spindle nose
spindle_z = Translate([spindle_z],c,0,'pivot-y','pivot-z')
# spindle head and z-slide move with z
spindle_z = Translate([spindle_z],hal,0,0,'joint.2.pos-fb')
# move spindle_z to z-home position
spindle_z = Translate([spindle_z], 0, 0, machine_zero_z)
# x carriage
EGO_X = Color([EGO_X],0.6,0.8,0.3,1)
spindle_xz = Collection([
             spindle_z,
             EGO_X
             ])
# spindle head, z-slide and x-slide move with x
spindle_xz = Translate([spindle_xz],hal,'joint.0.pos-fb',0,0)
# move spindle_xz to x-home position
spindle_xz = Translate([spindle_xz], machine_zero_x, 0, 0)
# end toolside

# start workside
# Create the work tracker (top center of the rotary table)
work = Capture('work')

# work_plane
# create the g68.n work_plane using position and orientation from the status channel
work_plane_data = ( s,"g68_rotation[0]", "g68_rotation[3]", "g68_rotation[6]",
                      "g68_rotation[2]", "g68_rotation[5]", "g68_rotation[8]",
                      "g68_offset[0]", "g68_offset[1]", "g68_offset[2]")

work_plane_coords = Axes(c,'300*{scale-coords}')
work_plane = Color([Plane(300)],c,1,0,1,0.3)
work_plane = Collection([work_plane,
                         work_plane_coords
                         ])
work_plane = MatrixTransform([work_plane],*work_plane_data)
work_plane = Translate([work_plane],s,"g92_offset[0]", "g92_offset[1]", "g92_offset[2]")
work_plane = Rotate([work_plane],s,'rotation_yz',1,0,0)
work_plane = Rotate([work_plane],s,'rotation_zx',0,1,0)
work_plane = Rotate([work_plane],s,'rotation_xy',0,0,1)
work_plane = Translate([work_plane],s,"g5x_offset[0]", "g5x_offset[1]", "g5x_offset[2]")
work_plane = Translate([work_plane], machine_zero_x, 0, machine_zero_z)
work_plane = Translate([work_plane], 0, machine_zero_y,0)
work_plane = Rotate([work_plane],c, 'g68-work-rotation', 0,0,-1)
arc_work_plane = Arrow_circular(c, "g68-work-origin-x","g68-work-origin-y","g68-work-origin-z",
                                   0,0,0,
                                   0,0,1,
                                   "-{g68-work-rotation}",
                                   3)
work_plane = Collection ([work_plane, arc_work_plane])
work_plane = Translate([work_plane], 0, -machine_zero_y,0)
work_plane = Translate([work_plane],s,0,"-{joint_position[1]}",0)
work_plane = Scale([work_plane],s,True,'g68_active',1,0)
# /work_plane

# work_coords
# create work coordinates that represent G52/G92 offset
g92_coords = Axes(100)
line_g68_offset = Arrow(s,0,0,0,"g68_offset[0]","g68_offset[1]","g68_offset[2]",3)
g92_coords = Collection([g92_coords, line_g68_offset])
g92_coords = Translate([g92_coords],s,"g92_offset[0]", "g92_offset[1]", "g92_offset[2]")
line_g92_offset = Arrow(s,0,0,0,"g92_offset[0]","g92_offset[1]","g92_offset[2]",3)
g92_coords = Collection([g92_coords, line_g92_offset])
# create work coordinates that represent work_offset offset
work_coords = Axes(100)
work_coords = Collection([work_coords, g92_coords])
# rotate to match the G10 L2 I J K(R) offsets values
work_coords = Rotate([work_coords],s,'rotation_yz',1,0,0)
work_coords = Rotate([work_coords],s,'rotation_zx',0,1,0)
work_coords = Rotate([work_coords],s,'rotation_xy',0,0,1)
work_coords = Translate([work_coords],s,"g5x_offset[0]","{g5x_offset[1]}-{joint_position[1]}","g5x_offset[2]")
line_work_offset = Arrow(s,0,0,0,"g5x_offset[0]","{g5x_offset[1]}-{joint_position[1]}","g5x_offset[2]",3)
work_coords = Collection([work_coords, line_work_offset])
work_coords = Translate([work_coords], machine_zero_x, 0, machine_zero_z)
# /work_coords


# Create geometry for the work piece
work_piece = Translate([Box(600,600,600)],0,0,300).set_group('work')
# Create rotary table
EGO_C = Color([EGO_C],0.1,0.7,0.9,1)
rotary_table_c = Collection([
                 work_piece,
                 EGO_C,
                 work
                 ])
# Rotary table and work roatae with axis C'
rotary_table_c = Rotate([rotary_table_c],hal,'joint.4.pos-fb',0,0,-1)
# y-carriage that carries the rotary table
EGO_Y = Color([EGO_Y],0.2,0.2,0.2,1)
table = Collection([
        rotary_table_c,
        EGO_Y
        ])
# Table moves with y axis
table = Translate([table],hal,0,'-{joint.1.pos-fb}',0,0).set_group('table')
# move table to y-home position
table = Translate([table], 0, -machine_zero_y, 0)
#/work-side

# Create machine base
base = Color([EGO_BC],0.3,0.3,0.3,1)



model = Collection([
        machine_coords,
        spindle_xz,
        table,
        base,
        work_coords,
        work_plane
        ])

#hud
myhud = Hud(color='mint',opacity=0.4) # This will always be displayed
myhud.add_txt('DMU-160-P')
myhud.add_txt('---------')
myhud.add_txt('')
myhud.add_txt('Kinematic Mode: IDENTITY',0)
myhud.add_txt('Kinematic Mode: TCP',1)
myhud.add_txt('')
myhud.add_txt('TWP-Orientation:')
myhud.add_txt('   X             Z')
myhud.add_var('{:6.3f}        {:6.3f}',s,["g68_rotation[0]","g68_rotation[2]"])
myhud.add_var('{:6.3f}        {:6.3f}',s,["g68_rotation[3]","g68_rotation[5]"])
myhud.add_var('{:6.3f}        {:6.3f}',s,["g68_rotation[6]","g68_rotation[8]"])
myhud.show_tag_eq_var_offs(hal,'motion.kins-type')
myhud.extra_text_enable = True
#/hud


main(options,
     c,
     model,
     huds=[myhud],
     function = myfunction,
     window_width=1400,
     window_height=1000,
     window_title = 'Vtk_Vismach Tutorial')

