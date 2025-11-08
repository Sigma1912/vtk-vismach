#!/usr/bin/env python3

import sys
import re
import tempfile
import shutil
import os
import math
import vtk

from collections import OrderedDict

from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QApplication

import linuxcnc
import gcode

from .basecanon import StatCanon

COLOR_MAP = {
    'traverse': QColor(200, 35, 35, 255),
    'arcfeed': QColor(110, 110, 255, 255),
    'feed': QColor(210, 210, 255, 255),
    'dwell': QColor(0, 0, 255, 255),
    'user': QColor(0, 100, 255, 255)
}


class AxesActor(vtk.vtkAxesActor):
    def __init__(self, *args):
        super(AxesActor, self).__init__()
        if  True:
            self.axes_length = 1.5
        else:
            self.axes_length = 0.375
        transform = vtk.vtkTransform()
        transform.Translate(0.0, 0.0, 0.0)  # Z up
        self.SetUserTransform(transform)
        self.AxisLabelsOff()
        self.SetShaftTypeToLine()
        self.SetTipTypeToCone()

        # Lathe modes
        if False:
            self.SetTotalLength(self.axes_length, 0, self.axes_length)
        else:
            self.SetTotalLength(self.axes_length, self.axes_length, self.axes_length)


class PathActor(vtk.vtkActor):
    def __init__(self, *args):
        super(PathActor, self).__init__()
        self.origin_index = None
        self.origin_cords = None
        if True : #self._datasource.isMachineMetric():
            self.axes_length = 1
        else:
            self.axes_length = .25
        self.axes_actor = AxesActor()
        if False: #self._datasource.isMachineLathe():
            self.axes_actor.SetTotalLength(self.axes_length, 0, self.axes_length)
        else:
            self.axes_actor.SetTotalLength(self.axes_length, self.axes_length, self.axes_length)
        # Create a vtkUnsignedCharArray container and store the colors in it
        self.colors = vtk.vtkUnsignedCharArray()
        self.colors.SetNumberOfComponents(4)
        self.points = vtk.vtkPoints()
        self.lines = vtk.vtkCellArray()
        self.poly_data = vtk.vtkPolyData()
        self.data_mapper = vtk.vtkPolyDataMapper()

    def set_origin_index(self, index):
        self.origin_index = index

    def get_origin_index(self):
        return self.origin_index

    def set_orgin_coords(self, *cords):
        self.origin_cords = cords

    def get_origin_coords(self):
        return self.origin_cords

    def get_axes_actor(self):
        return self.axes_actor


class VtkCanon(StatCanon):
    def __init__(self, colors=COLOR_MAP, *args, **kwargs):
        super(VtkCanon, self).__init__(*args, **kwargs)
        self.stat = linuxcnc.stat()
        self.stat.poll()
        self.path_colors = colors
        self.path_actors = OrderedDict()
        self.path_points = OrderedDict()
        self.initial_wcs_offsets = OrderedDict()
        self.paths_start_point = OrderedDict()
        self.paths_angle_point = OrderedDict()
        self.paths_end_point = OrderedDict()
        self.path_start_point = OrderedDict()
        self.paths_angle_points = OrderedDict()
        self.path_end_point = OrderedDict()
        self.active_wcs_index = 0
        #self.active_wcs_index = self._datasource.getActiveWcsIndex()
        self.active_rotation = 0
        #self.active_rotation = self._datasource.getRotationOfActiveWcs()
        #g5x = (0,0,0,0,0,0,0,0,0)
        #g5x = self._datasource.getActiveWcsOffsets()
        s = self.stat
        g5x = (s.g5x_offset[0],s.g5x_offset[1],s.g5x_offset[2],0,0,0,0,0,0)
        print(f" G5x offsets = {g5x}")
        print(f"XY Rotation = {self.active_rotation}")
        # ensure Canon has correct starting offsets per var file
        super().set_g5x_offset(self.active_wcs_index, g5x[0],g5x[1],g5x[2],g5x[3],g5x[4],g5x[5],g5x[6],g5x[7],g5x[8])
        self.foam_z = 0.0
        self.foam_w = 0.0

    def comment(self, comment):
        print("G-code Comment: {}".format(comment))
        items = comment.lower().split(',', 1)
        if len(items) > 0 and items[0] in ['axis', 'backplot']:
            cmd = items[1].strip()
            if cmd == "hide":
                self.suppress += 1
            elif cmd == "show":
                self.suppress -= 1
            elif cmd == 'stop':
                print("Backplot generation aborted.")
                raise KeyboardInterrupt
            elif cmd.startswith("xy_z_pos"):
                self.foam_z = float(cmd.split(',')[1])

            elif cmd.startswith("uv_z_pos"):
                self.foam_w = float(cmd.split(',')[1])

    def message(self, msg):
        print("G-code Message: {}".format(msg))

    def set_g5x_offset(self, index, x, y, z, a, b, c, u, v, w):
        # ensure the passed values get set on 'self' via super
        print("----------------------------------")
        print("--------- set_g5x_offset ---------")
        print("----------------------------------")
        super().set_g5x_offset(index, x, y, z, a, b, c, u, v, w)
        new_wcs = index - 1  # this index counts also G53 so we need to do -1
        print("---------received wcs change: {}".format(new_wcs))
        print("--------- wcs values: x, y, z, a, b, c, u, v, w")
        print(f"--------- wcs values: {x}, {y}, {z}, {a}, {b}, {c}, {u}, {v}, {w}")
        try:
            if new_wcs not in list(self.path_points.keys()):
                self.path_actors[new_wcs] = PathActor()
                self.path_points[new_wcs] = list()
                self.initial_wcs_offsets[new_wcs] = (x, y, z, a, b, c, u, v, w)
        except Exception as e:
            print(f"CANON ERROR {e}")
        self.active_wcs_index = new_wcs

    def set_xy_rotation(self, rotation):
        self.rotation_xy = 0.0
        theta = math.radians(0.0)
        self.rotation_cos = math.cos(theta)
        self.rotation_sin = math.sin(theta)

    def add_path_point(self, line_type, start_point, end_point):
        #print("----------------------------------")
        #print("--------- add_path_point ---------")
        #print("----------------------------------")
        # As the points come through with the active wcs offsets baked in
        # remove them to allow vtk setusertransforms to work correctly.
        # These transforms apply wcs offsets for us in VTK
        adj_start_point = list(start_point)
        adj_end_point = list(end_point)
        # check to see if active wcs is in the path_actor list.
        if self.active_wcs_index not in list(self.path_actors.keys()):
            #self.path_actors[self.active_wcs_index] = PathActor(self._datasource)
            self.path_actors[self.active_wcs_index] =[]
        for count,value in enumerate(self.initial_wcs_offsets[self.active_wcs_index]):
            adj_start_point[count] -= value
            adj_end_point[count] -= value
        line = [tuple(adj_start_point), tuple(adj_end_point)]
        self.path_points.get(self.active_wcs_index).append((line_type, line))

    def draw_lines(self):
        # Used to draw the lines of the loaded program
        print("------------------------------")
        print("--------- draw_lines ---------")
        print("------------------------------")
        print("--------- path points size: {}".format(sys.getsizeof(self.path_points)))
        print("--------- path points length: {}".format(len(self.path_points)))
        #unit_factor = 25.4 if self._datasource.isMachineMetric() else 1
        unit_factor = 25.4 if True else 1
        paths_count = 0
        prev_wcs_index = 0
        for wcs_index, data in self.path_points.items():
            path_actor = self.path_actors.get(wcs_index)
            if path_actor is not None:
                first_point = False
                angle_point = None
                last_point = None
                point_count = 0
                for line_type, line_data in data:
                    start_point = line_data[0]
                    end_point = line_data[1]
                    if len(self.path_actors) > 1:
                        # skip rapids from original path offsets. This actually means previous wcs
                        # >1 path_actors means more than one g5x in use in the file.
                        if (paths_count > 0) and (point_count == 0) and (line_type == "traverse"):
                            continue
                        if point_count == 0:
                            position = [start_point[0] * unit_factor,
                                        start_point[1] * unit_factor,
                                        start_point[2] * unit_factor]
                            # Get start point for a transition line between different WCS
                            self.path_start_point[prev_wcs_index] = position
                    path_actor.points.InsertNextPoint(end_point[0] * unit_factor,
                                                        end_point[1] * unit_factor,
                                                        end_point[2] * unit_factor)
                    path_actor.points.InsertNextPoint(start_point[0] * unit_factor,
                                                        start_point[1] * unit_factor,
                                                        start_point[2] * unit_factor)
                    path_actor.colors.InsertNextTypedTuple(self.path_colors.get(line_type).getRgb())
                    line = vtk.vtkLine()
                    line.GetPointIds().SetId(0, point_count)
                    line.GetPointIds().SetId(1, point_count + 1)
                    path_actor.lines.InsertNextCell(line)
                    point_count += 2
                    last_point = end_point
                if (len(self.path_actors) > 1) and (last_point is not None):
                    # Store the last point of the part as first point of the rapid line
                    position = [last_point[0] * unit_factor,
                                last_point[1] * unit_factor,
                                last_point[2] * unit_factor]
                    # Get end point for a transition line between different WCS
                    self.path_end_point[wcs_index] = position
                else:
                    self.path_end_point[wcs_index] = None
                # free up memory, lots of it for big files
                self.path_points[wcs_index].clear()
                QApplication.processEvents()
                path_actor.poly_data.SetPoints(path_actor.points)
                path_actor.poly_data.SetLines(path_actor.lines)
                path_actor.poly_data.GetCellData().SetScalars(path_actor.colors)
                path_actor.data_mapper.SetInputData(path_actor.poly_data)
                path_actor.data_mapper.Update()
                path_actor.SetMapper(path_actor.data_mapper)
            paths_count += 1
            prev_wcs_index = wcs_index
        print("----------------------------------")
        print("--------- draw_lines END ---------")
        print("----------------------------------")

    def get_path_actors(self):
        return self.path_actors

    # Methods get_offsets_start_point and get_offsets_end_point provide
    # the start and end points for a transition line between different WCS
    def get_offsets_start_point(self):
        return self.path_start_point

    def get_offsets_end_point(self):
        return self.path_end_point

    def get_foam(self):
        return self.foam_z, self.foam_w



class VtkPreview(vtk.vtkAssembly):
    def __init__(self):
        super().__init__()
        self.work = (0,0,0)
        self.tool = (0,0,0)
        self.first_update_done = False
        self.transformation = vtk.vtkTransform()
        self.SetUserTransform(self.transformation)
        self.canon = None
        self.stat = linuxcnc.stat()
        self.stat.poll()
        self._current_file = self.stat.file
        inifile = os.environ.get('INI_FILE_NAME', '/dev/null')
        self.ini = linuxcnc.ini(inifile)
        self.config_dir = os.path.dirname(inifile)
        temp = self.ini.find("EMCIO", "RANDOM_TOOLCHANGER")
        self.random = int(temp or 0)
        temp = self.ini.find("DISPLAY", "GEOMETRY") or 'XYZ'
        self.geometry = temp.upper()
        temp = self.ini.find("DISPLAY", "LATHE") or "false"
        self.lathe_option = temp.lower() in ["1", "true", "yes"]
        temp = self.ini.find("RS274NGC", "PARAMETER_FILE") or "linuxcnc.var"
        self.parameter_file = os.path.join(self.config_dir, temp)
        self.temp_parameter_file = os.path.join(self.parameter_file + '.temp')

    def update(self):
        s = self.stat
        s.poll()
        if  s.file and (self._current_file != s.file):
            print('self._current_file, s.file: ',self._current_file, s.file)
            if self.canon:
                print('cleanup the scene')
                # Cleanup the scene, remove any previous actors if any.
                # Do this for each WCS.
                for wcs_index, path_actor in list(self.path_actors.items()):
                    print('wcs_index ',wcs_index)
                    print('Removing part')
                    axes = path_actor.get_axes_actor()
                    self.RemovePart(path_actor)
                    self.RemovePart(axes)
                self.path_actors.clear()
            self.canon = VtkCanon()
            self.load(s.file)
            self.canon.draw_lines()
            self.path_actors = self.canon.get_path_actors()
            self._current_file = s.file
            for wcs_index, path_actor in list(self.path_actors.items()):
                print('wcs_index ',wcs_index)
                axes = path_actor.get_axes_actor()
                offset = []
                for axis in range(0,9):
                    offset.append(s.g5x_offset[axis])
                path_actor.SetPosition(offset[:3])
                axes.SetPosition(offset[:3])
                self.AddPart(path_actor)
                self.AddPart(axes)
        if not self.first_update_done:
            print('self.work, self.tool: ',self.work, self.tool)
            self.SetPosition(self.work)
            self.AddPosition(self.tool)
            self.first_update_done = True


    def load(self, filename=None):
        if self.canon is None:
            return
        filename = filename or self.last_filename
        if filename is None:
            self.stat.poll()
            filename = self.stat.file
        if filename is None or not os.path.isfile(filename):
            self.canon = None
            print("3D plot", "Can't load backplot, invalid file: {}".format(filename))
        self.last_filename = filename
        if os.path.exists(self.parameter_file):
            shutil.copy(self.parameter_file, self.temp_parameter_file)
        self.canon.parameter_file = self.temp_parameter_file
        # Some initialization g-code to set the units and optional user code
        unitcode = "G%d" % (20 + (self.stat.linear_units == 1))
        initcode = self.ini.find("RS274NGC", "RS274NGC_STARTUP_CODE") or ""
        # THIS IS WHERE IT ALL HAPPENS: load_preview will execute the code,
        # call back to the canon with motion commands, and record a history
        # of all the movements.
        try:
            result, seq = gcode.parse(filename, self.canon, unitcode, initcode)
            if result > gcode.MIN_ERROR:
                msg = gcode.strerror(result)
                fname = os.path.basename(filename)
                print("3D plot", "Error in {} line {}\n{}".format(fname, seq - 1, msg))
        except KeyboardInterrupt:
            # probably raised by an (AXIS, stop) comment in the G-code file abort generating the backplot
            pass
        except Exception as e:
            print(f"CANON ERROR {e}")
        # clean up temp var file and the backup
        os.unlink(self.temp_parameter_file)
        os.unlink(self.temp_parameter_file + '.bak')
