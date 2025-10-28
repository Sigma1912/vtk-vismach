#!/usr/bin/env python3

from vtk_vismach import *
import hal
from math import sin, cos, radians, degrees
import sys



try: # Expect files in working directory
    # NOTE to run this file as standalone python script absolute paths might have to be used to find the stl files
    # create the work piece from file
    work_piece = ReadPolyData("./vismach-stl-files/work_piece_1.stl")
    EGO_BC = ReadPolyData("./vismach-stl-files/160Pbase.stl")
    EGO_B = ReadPolyData("./vismach-stl-files/160PB.stl")
    EGO_C = ReadPolyData("./vismach-stl-files/160PC.stl")
    EGO_X = ReadPolyData("./vismach-stl-files/160PX.stl")
    EGO_Y = ReadPolyData("./vismach-stl-files/160PY.stl")
    EGO_Z = ReadPolyData("./vismach-stl-files/160PZ.stl")
except Exception as detail:
    print(detail)
    raise SystemExit("vtk-dmu-160-p-gui requires stl files in working directory")

#for setting in sys.argv[1:]: exec(setting)

c = hal.component("vtk-dmu-160-p-gui")
# nutation-angle
c.newpin("nutation_angle", hal.HAL_FLOAT, hal.HAL_IN)
# virtual tool rotation
c.newpin("virtual_rotation", hal.HAL_FLOAT, hal.HAL_IN)
# geometric offsets in the spindle-rotary-assembly
c.newpin("pivot_y", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("pivot_z", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("offset_x", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("offset_z", hal.HAL_FLOAT, hal.HAL_IN)
# rot-point offsets distances from pivot point ( spindle AB)
# to rotation axis ( table C)
c.newpin("rot_axis_x", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("rot_axis_y", hal.HAL_FLOAT, hal.HAL_IN)
# active work offset values (ie g54,g55...)
c.newpin("work_offset_x", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("work_offset_y", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("work_offset_z", hal.HAL_FLOAT, hal.HAL_IN)
# selected kinematics
c.newpin("kinstype_select", hal.HAL_FLOAT, hal.HAL_IN)
# work piece show/hide
c.newpin("hide_work_piece_1", hal.HAL_BIT, hal.HAL_IN)
# spindle body show/hide
c.newpin("hide_spindle_body", hal.HAL_BIT, hal.HAL_IN)
# work piece show/hide
c.newpin("hide_somethingelse", hal.HAL_BIT, hal.HAL_IN)
# scale coordinate system indicators
c.newpin("scale_coords", hal.HAL_FLOAT, hal.HAL_IN)
# twp pins
c.newpin("twp_status", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("twp_defined", hal.HAL_BIT, hal.HAL_IN)
c.newpin("twp_active", hal.HAL_BIT, hal.HAL_IN)
# the origin of the twp plane display needs to be independent of the work offsets
# because those offsets change as the kinematic switches between world and tool
c.newpin("twp_ox_world", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("twp_oy_world", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("twp_oz_world", hal.HAL_FLOAT, hal.HAL_IN)
# origin of the twp
c.newpin("twp_ox", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("twp_oy", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("twp_oz", hal.HAL_FLOAT, hal.HAL_IN)
# normal vector defining the temporary work plane as set by user with G68.2
c.newpin("twp_zx", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("twp_zy", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("twp_zz", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("twp_xx", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("twp_xy", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("twp_xz", hal.HAL_FLOAT, hal.HAL_IN)
# twp plane rotation, This is only used to show euler rotations of the plane in vismach
# Currently not used anymore
c.newpin("rot_order", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("rot_th1", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("rot_th2", hal.HAL_FLOAT, hal.HAL_IN)
c.newpin("rot_th3", hal.HAL_FLOAT, hal.HAL_IN)
c.ready()


# used to rotate parts around the nutation axis
class HalNutate(Collection):
    def update(self):
        th,x,y,z = self.coords()
        # override x,y,z
        x = 0
        y = sin(radians(c.nutation_angle))
        z = cos(radians(c.nutation_angle))
        self.transformation = vtk.vtkTransform()
        self.transformation.PreMultiply()
        self.transformation.RotateWXYZ(th,x,y,z)

    def transform(self):
        self.SetUserTransform(self.transformation)


# Machine zero as measured from center surface of the rotary c table
machine_zero_x = -1000
machine_zero_y =  1000
machine_zero_z =  1000

# Create an indicator for the machine reference coordinates
machine_axes = Axes(scale=200)
machine_axes = Translate([machine_axes], machine_zero_x, 0, machine_zero_z)

# start toolside
# Create the tooltip tracker (tool control point)
tooltip = Capture()
# Create an indicator for the tool coordinates
tool_axes = Axes(scale=100)
tool_shape = Collection([
                CylinderZ(hal,"motion.tooloffset.z", ('halui.tool.diameter', 0.5)),
                # this indicates the spindle nose when no tool-offset is active
                CylinderZ(-0.1, 0),
                ])
tool = Collection([
                    tooltip,
                    tool_axes,
                    tool_shape
                    ])
tool = Rotate([tool],c,"virtual_rotation",0,0,1)
tool = Translate([tool],hal,0,0,("motion.tooloffset.z",-1))
tool = Translate([tool],c,0,("pivot_y",-1),("pivot_z",-1))
tool = Color([tool],(1,0,1),1)
# create spindle head
EGO_B = Color([EGO_B],(0.7,0.7,0),1)
# rotate the nutation joint to the nutation angle, 90° should have the nutation axis in the horizontal plane
EGO_B = Rotate([EGO_B],-90,-1,0,0)
EGO_B = Rotate([EGO_B],c,"nutation_angle",-1,0,0)

spindle_assembly = Collection([
                   tool,
                   EGO_B,
                   #ind_pivot_point,
                   #ind_pivot_z,
                   #ind_pivot_y
                   ])

# create HAL-link for b-axis rotational joint"
spindle_assembly = HalNutate([spindle_assembly],hal,"joint.3.pos-fb",0,-1,0)
# Z carriage
EGO_Z = Color([EGO_Z],(1,0.5,0),1)
EGO_Z = Translate([EGO_Z], 0, 0, -759.5)
spindle_z = Collection([
             spindle_assembly,
             EGO_Z
             ])
# move by the values set for the pivot lengths so the vismach origin is in the center of the spindle nose
spindle_z = Translate([spindle_z],c,0,"pivot_y","pivot_z")
# spindle head and z-slide move with z
spindle_z = Translate([spindle_z],hal,0,0,"joint.2.pos-fb")
# move spindle_z to z-home position
spindle_z = Translate([spindle_z], 0, 0, machine_zero_z)
# x carriage
EGO_X = Color([EGO_X],(0.6,0.8,0.3),1)
spindle_xz = Collection([
             spindle_z,
             EGO_X
             ])
# spindle head, z-slide and x-slide move with x
spindle_xz = Translate([spindle_xz],hal,"joint.0.pos-fb",0,0)
# move spindle_xz to x-home position
spindle_xz = Translate([spindle_xz], machine_zero_x, 0, 0)
# end toolside

# start workside
# Create the work tracker (top center of the rotary table)
work = Capture()
# Create an indicator for the tool coordinates
work_axes = Axes(scale=100)
# Create an indicator for the defined Tilted Work Plane (TWP)
work_plane_defined =  HalGridFromNormalAndDirection(c,
        "twp_ox", "twp_oy", "twp_oz",
        "twp_xx", "twp_xy", "twp_xz",
        "twp_zx", "twp_zy", "twp_zz",
        300
        )
# for twp-defined = true, we show the plane in gray
work_plane_defined = Color([work_plane_defined],(0.7,0.7,0.7),0.3)
work_plane_defined = HalShow([work_plane_defined],c,0,"twp_active")
# Create an indicator for the active Tilted Work Plane (TWP), same as active but with different color
work_plane_active =  HalGridFromNormalAndDirection(c,
        "twp_ox", "twp_oy", "twp_oz",
        "twp_xx", "twp_xy", "twp_xz",
        "twp_zx", "twp_zy", "twp_zz",
        300
        )
# for twp-active = true, we show the plane in pink
work_plane_active = Color([work_plane_active],(1,0,1),0.3)
work_plane_active = HalShow([work_plane_active],c,1,"twp_active")
# Create a coordinate system for the twp-plane
work_plane_coords =  HalCoordsFromNormalAndDirection(c,
        "twp_ox", "twp_oy", "twp_oz",
        "twp_xx", "twp_xy", "twp_xz",
        "twp_zx", "twp_zy", "twp_zz",
        300
        )
work_plane = Collection([work_plane_defined, work_plane_active, work_plane_coords])
# make the work_plane hidable
work_plane = HalShow([work_plane],c,1,"twp_defined")
# Create geometry for the work piece
work_piece = BoxCentered(600,600,600)
work_piece = Translate([work_piece],0,0,300)
# Make the work piece hidable
work_piece = HalShow([work_piece],c,True,"hide_work_piece_1",0,1)
# Create rotary table
EGO_C = Color([EGO_C],(0.1,0.7,0.9),1)
rotary_table_c = Collection([
                 work_piece,
                 EGO_C,
                 work_plane,
                 work_axes,
                 work
                 ])
# Rotary table and work roatae with axis C"
rotary_table_c = Rotate([rotary_table_c],hal,"joint.4.pos-fb",0,0,-1)
# y-carriage that carries the rotary table
EGO_Y = Color([EGO_Y],(0.2,0.2,0.2),1)
table = Collection([
        rotary_table_c,
        EGO_Y
        ])
# Table moves with y axis
#table = HalTranslate_orig([table],c,"axis_y",0,-1,0)
table = Translate([table],hal,0,("joint.1.pos-fb",-1),0)
# move table to y-home position
table = Translate([table], 0, -machine_zero_y, 0)
#/work-side

# Create machine base
base = Color([EGO_BC],(0.3,0.3,0.3),1)

arrow = ArrowOriented(hal, machine_zero_x, 0,-machine_zero_z,0,("vismach.work_pos_y",1),0,50)
#arrow = Translate_orig([arrow], machine_zero_x, 0, machine_zero_z)

model = Collection([
        machine_axes,
        spindle_xz,
        table,
        base,
        arrow,
        CylinderX(hal,"joint.0.pos-fb",50),
        Box(hal,"joint.1.pos-fb",0,0,100,100,-100),
        Sphere(0,0,0,5),
        Line(hal,("joint.1.pos-fb",-1),100,100,-1000,-1000,1000,2),
        ArrowOriented(hal,0,0,0,"vismach.work_pos_x","vismach.work_pos_y","vismach.work_pos_z",50),
        CylinderOriented(hal,"joint.1.pos-fb",100,100,-1000,-1000,1000,50),
        ArrowOriented(c,0,0,0,"twp_ox_world","twp_oy_world","twp_oz_world",20),
        ])

#hud
myhud = Hud(color='mint') # This will always be displayed
myhud.add_txt("DMU-160-P")
myhud.add_txt("---------")
myhud.add_txt("")
myhud.add_txt("Kinematic Mode: IDENTITY",0)
myhud.add_txt("Kinematic Mode: TCP",1)
myhud.add_txt("Kinematic Mode: TOOL",2)
myhud.add_txt("")
myhud.add_txt("TWP-Orientation Vector X:")
myhud.add_pin("Xx: {:8.3f}",c,"twp_xx")
myhud.add_pin("Xy: {:8.3f}",c,"twp_xy")
myhud.add_pin("Xz: {:8.3f}",c,"twp_xz")
myhud.add_txt("")
myhud.add_txt("TWP-Orientation Vector Z:")
myhud.add_pin("Zx: {:8.3f}",c,"twp_zx")
myhud.add_pin("Zy: {:8.3f}",c,"twp_zy")
myhud.add_pin("Zz: {:8.3f}",c,"twp_zz")
myhud.show_tag_eq_pin_offs(hal,"motion.switchkins-type")

myhud2= Hud(c,"kinstype_select",2,"tomato", 0.4) # This is displayed when kintype is 2
myhud2.add_txt("")
myhud2.add_txt("")
myhud2.add_txt("")
myhud2.add_txt("")
myhud2.add_txt("")
myhud2.add_txt("")
myhud2.add_pin("    {:8.3f}",c,"twp_xx")
myhud2.add_pin("    {:8.3f}",c,"twp_xy")
myhud2.add_pin("    {:8.3f}",c,"twp_xz")
myhud2.add_txt("")
myhud2.add_txt("")
myhud2.add_pin("    {:8.3f}",c,"twp_zx")
myhud2.add_pin("    {:8.3f}",c,"twp_zy")
myhud2.add_pin("    {:8.3f}",c,"twp_zz")
myhud2.extra_text_enable = True
#/hud

main(c,model, tooltip, work, huds=[myhud,myhud2], window_width=1400, window_height=1000, window_title = "Vtk_Vismach Tutorial")
