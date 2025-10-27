# This is a modified version of 'vismach.py' to visualize Tilted Work Plane (TWP)
# Author: David mueller
# email: mueller_david@hotmail.com

import hal
import vtk  # this should be changed to only import the used modules
import os

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import Qt
from PyQt5.QtCore import QTimer


def update_passed_args(self, v):
    if isinstance(getattr(self,v), str):
        # class variable of the same name is a string (ie a pinname)
        halpin = getattr(self,v)
        if self.comp:
            # if the component has been passed then we need to get the value using that
            var = self.comp[halpin]
        else:
            # if the comp variable is None then we need to get the value through hal
            var = hal.get_value(halpin)
    else:
        # class variable of the same name is not a string (ie a constant value)
        var = getattr(self,v)
    return var


class CoordsBase(vtk.vtkActor):
    def __init__(self, *args):
        if args and isinstance(args[0], hal.component):
            self.comp = args[0]
            args = args[1:]
        elif args and args[0] == None:
            self.comp = None
            args = args[1:]
        else:
            self.comp = None
        self._coords = args
        if hasattr(self, "create"):
            self.create()

    def coords(self):
        return list(map(self._coord, self._coords))

    def _coord(self, v):
        if isinstance(v, str) and self.comp:
            return self.comp[v]
        elif isinstance(v, str) and not self.comp:
            return hal.get_value(v)
        else:
            return v


class Box(CoordsBase):
    def create(self, *args):
        x1, y1, z1, x2, y2, z2 = self.coords()
        if x1 > x2:
            tmp = x1
            x1 = x2
            x2 = tmp
        if y1 > y2:
            tmp = y1
            y1 = y2
            y2 = tmp
        if z1 > z2:
            tmp = z1
            z1 = z2
            z2 = tmp
        self.cube = vtk.vtkCubeSource()
        self.cube.SetXLength(x2-x1)
        self.cube.SetYLength(y2-y1)
        self.cube.SetZLength(z2-z1)
        self.cube.Update()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.cube.GetOutput())
        self.SetMapper(mapper)
        self.SetPosition(x1,y1,z1)
        self.AddPosition((x2-x1)/2,(y2-y1)/2,(z2-z1)/2)


    def update(self):
        x1, y1, z1, x2, y2, z2 = self.coords()
        if x1 > x2:
            tmp = x1
            x1 = x2
            x2 = tmp
        if y1 > y2:
            tmp = y1
            y1 = y2
            y2 = tmp
        if z1 > z2:
            tmp = z1
            z1 = z2
            z2 = tmp
        self.cube.SetXLength(x2-x1)
        self.cube.SetYLength(y2-y1)
        self.cube.SetZLength(z2-z1)
        self.cube.Update()
        self.SetPosition(x1,y1,z1)
        self.AddPosition((x2-x1)/2,(y2-y1)/2,(z2-z1)/2)


# specify the width in X and Y, and the height in Z
# the box is centered on the origin
class BoxCentered(CoordsBase):
    def create(self, *args):
        xw, yw, zw = self.coords()
        self.cube = vtk.vtkCubeSource()
        self.cube.SetXLength(xw)
        self.cube.SetYLength(yw)
        self.cube.SetZLength(zw)
        self.cube.Update()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.cube.GetOutput())
        self.SetMapper(mapper)

    def update(self):
        xw, yw, zw = self.coords()
        self.cube.SetXLength(xw)
        self.cube.SetYLength(yw)
        self.cube.SetZLength(zw)
        self.cube.Update()

# specify the width in X and Y, and the height in Z
# the box is centered on the origin
class Sphere(CoordsBase):
    def create(self, *args):
        x, y, z, r = self.coords()
        self.sphere = vtk.vtkSphereSource()
        self.sphere.SetRadius(r)
        self.sphere.Update()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.sphere.GetOutput())
        self.SetMapper(mapper)
        self.SetPosition(x,y,z)

    def update(self):
        x, y, z, r = self.coords()
        self.sphere.SetRadius(r)
        self.sphere.Update()
        self.SetPosition(x,y,z)


# Create cylinder along Y axis (default direction for vtkCylinderSource)
class CylinderY(CoordsBase):
    def create(self, *args):
        length, radius = self.coords()
        self.cylinder = vtk.vtkCylinderSource()
        self.cylinder.SetRadius(radius)
        self.cylinder.SetHeight(length)
        self.cylinder.SetResolution(10)
        self.cylinder.Update()
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputData(self.cylinder.GetOutput())
        self.SetMapper(self.mapper)

    def update(self):
        length, radius = self.coords()
        self.cylinder.SetRadius(radius)
        self.cylinder.SetHeight(length)
        self.SetUserTransform(vtk.vtkTransform())
        self.GetUserTransform().Translate(0,length/2,0)
        self.cylinder.Update()


# Create cylinder along Z axis
class CylinderZ(CylinderY):
    def create(self, *args):
        super().create(self)
        self.RotateX(90)

    def update(self):
        length, radius = self.coords()
        self.cylinder.SetRadius(radius)
        self.cylinder.SetHeight(length)
        self.SetUserTransform(vtk.vtkTransform())
        self.GetUserTransform().Translate(0,0,length/2)
        self.cylinder.Update()


# Creates a 3d cylinder from (xs,ys,zs) to (xe,ye,ze)
class CylinderOriented(CoordsBase):
    def create(self):
        self.update()

    def update(self):
        xs, ys, zs, xe, ye, ze, radius = self.coords()
        resolution = 10
        # Create a cylinder (cylinders are created along Y axis by default)
        # Cylinder center is in the middle of the cylinder
        cylinderSource = vtk.vtkCylinderSource()
        cylinderSource.SetRadius(radius)
        cylinderSource.SetResolution(resolution)
        # Generate a random start and end point
        startPoint = [xs, ys, zs]
        endPoint = [xe, ye, ze]
        # Compute a basis
        normalizedX = [0] * 3
        normalizedY = [0] * 3
        normalizedZ = [0] * 3
        # The X axis is a vector from start to end
        vtk.vtkMath.Subtract(endPoint, startPoint, normalizedX)
        length = vtk.vtkMath.Norm(normalizedX)
        vtk.vtkMath.Normalize(normalizedX)
        vtk.vtkMath.Cross(normalizedX, [0,0,1], normalizedZ)
        vtk.vtkMath.Normalize(normalizedZ)
        # The Y axis is Z cross X
        vtk.vtkMath.Cross(normalizedZ, normalizedX, normalizedY)
        matrix = vtk.vtkMatrix4x4()
        # Create the direction cosine matrix
        matrix.Identity()
        for i in range(0, 3):
            matrix.SetElement(i, 0, normalizedX[i])
            matrix.SetElement(i, 1, normalizedY[i])
            matrix.SetElement(i, 2, normalizedZ[i])
        # Apply the transforms
        transform = vtk.vtkTransform()
        transform.Translate(startPoint)  # translate to starting point
        transform.Concatenate(matrix)  # apply direction cosines
        transform.RotateZ(-90.0)  # align cylinder to x axis
        transform.Scale(1.0, length, 1.0)  # scale along the height vector
        transform.Translate(0, .5, 0)  # translate to start of cylinder
        # mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(cylinderSource.GetOutputPort())
        self.SetUserMatrix(transform.GetMatrix())
        self.SetMapper(mapper)


# Creates a 3d arrow pointing from (xs,ys,zs) to (xe,ye,ze)
class ArrowOriented(CoordsBase):
    def create(self):
        self.update()

    def update(self):
        xs, ys, zs, xe, ye, ze, radius = self.coords()
        resolution = 10
        # Create arrow (arrows are created along the X axis by default)
        cylinderSource = vtk.vtkArrowSource()
        cylinderSource.SetShaftRadius(radius)
        #cylinderSource.SetTipLength(radius*10/length)
        cylinderSource.SetTipRadius(radius*3)
        cylinderSource.SetTipResolution(resolution)
        # Generate a random start and end point
        startPoint = [xs, ys, zs]
        endPoint = [xe, ye, ze]
        # Compute a basis
        normalizedX = [0] * 3
        normalizedY = [0] * 3
        normalizedZ = [0] * 3
        # The X axis is a vector from start to end
        vtk.vtkMath.Subtract(endPoint, startPoint, normalizedX)
        length = vtk.vtkMath.Norm(normalizedX)
        if length < 0.1: length = 1
        cylinderSource.SetTipLength(radius*2/length)
        vtk.vtkMath.Normalize(normalizedX)
        vtk.vtkMath.Cross(normalizedX, [0,0,1], normalizedZ)
        vtk.vtkMath.Normalize(normalizedZ)
        # The Y axis is Z cross X
        vtk.vtkMath.Cross(normalizedZ, normalizedX, normalizedY)
        matrix = vtk.vtkMatrix4x4()
        # Create the direction cosine matrix
        matrix.Identity()
        for i in range(0, 3):
            matrix.SetElement(i, 0, normalizedX[i])
            matrix.SetElement(i, 1, normalizedY[i])
            matrix.SetElement(i, 2, normalizedZ[i])
        # Apply the transforms
        transform = vtk.vtkTransform()
        transform.Translate(startPoint)  # translate to starting point
        transform.Concatenate(matrix)  # apply direction cosines
        transform.Scale(length, 1.0, 1.0)  # scale along the height vector
        transform.Translate(0, .5, 0)  # translate to start of cylinder
        # mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(cylinderSource.GetOutputPort())
        self.SetUserMatrix(transform.GetMatrix())
        self.SetMapper(mapper)


# function to crate the reader for 3d geometry files
def make_reader(file_name):
    path, extension = os.path.splitext(file_name)
    extension = extension.lower()
    if extension == '.ply':
        reader = vtk.vtkPLYReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.vtp':
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.obj':
        reader = vtk.vtkOBJReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.stl':
        reader = vtk.vtkSTLReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.vtk':
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.g':
        reader = vtk.vtkBYUReader()
        reader.SetGeometryFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    else:
        print("ReadPolyData Error: Unable to read file ", file_name)
    return reader


# Loads 3D geometry on startup
class ReadPolyData(vtk.vtkActor):
    def __init__(self, geometry_file):
        mapper = vtk.vtkPolyDataMapper()
        if not os.path.isfile(geometry_file):
            print("Vtk_Vismach HalReadPolyData Error: Unable to read file ", geometry_file)
            # create a dummy sphere instead
            sphereSource = vtk.vtkSphereSource()
            sphereSource.SetCenter(0.0, 0.0, 0.0)
            sphereSource.SetRadius(0.5)
            mapper.SetInputConnection(sphereSource.GetOutputPort())
        else:
            # create from stl file
            reader = make_reader(geometry_file)
            mapper.SetInputConnection(reader.GetOutputPort())
        self.SetMapper(mapper)
        self.SetUserTransform(vtk.vtkTransform())
        # Avoid visible backfaces on Linux with some video cards like intel
        # From: https://stackoverflow.com/questions/51357630/vtk-rendering-not-working-as-expected-inside-pyqt?rq=1#comment89720589_51360335
        self.GetProperty().SetBackfaceCulling(1)


# Dynamically loads 3D geometry depending on the self.comp[self.var] value
class HalReadPolyData(vtk.vtkActor):
    def __init__(self, comp, var, directory):
        self.comp = comp            # instance of the halcomponent used in the model
        self.var = var              # name of the file (integer or float)
        self.directory = directory  # directory holding the 3D geometry files
        self.SetUserTransform(vtk.vtkTransform())

    def update(self):
        for v in ['var']:
            # create variable from list and update from class variables of the same name
            globals()[v] = update_passed_args(self, v)
        self.tool_nr = var
        self.filename = self.directory + str(self.tool_nr) + ".stl"
        # Create a mapper
        mapper = vtk.vtkPolyDataMapper()
        if not os.path.isfile(self.filename):
            print("Vtk_Vismach HalReadPolyData Error: Unable to read file ", self.filename)
            # create a dummy sphere instead
            sphereSource = vtk.vtkSphereSource()
            sphereSource.SetCenter(0.0, 0.0, 0.0)
            sphereSource.SetRadius(0.5)
            mapper.SetInputConnection(sphereSource.GetOutputPort())
        else:
            # create from stl file
            reader = make_reader(self.filename)
            mapper.SetInputConnection(reader.GetOutputPort())
        self.SetMapper(mapper)
        # Avoid visible backfaces on Linux with some video cards like intel
        # From: https://stackoverflow.com/questions/51357630/vtk-rendering-not-working-as-expected-inside-pyqt?rq=1#comment89720589_51360335
        self.GetProperty().SetBackfaceCulling(1)


# Draws a polyline showing the path of 'tooltip' with respect to 'work'
class Plotter(vtk.vtkActor):
    def __init__(self, comp, work, tooltip, clear, color="magenta"):
        self.comp = comp            # instance of the halcomponent used in the model
        self.tooltip = tooltip      # Capture object with '.matrix' holding the current transformation tool->world
        self.work = work            # Capture object with '.matrix' holding the current transformation work->world
        self.clear = clear          # halpin that clears the backplot
        self.color = color          # color of backplot in eiter nomalized RGB or one of vtkNamedColors
        # We initialize at the origin, this is cleared and set to the actual
        # machine reference position during the 1. update loop
        self.setup_points([0,0,0])
        self.initial_run = True

    def setup_points(self, pos):
        self.index = 0
        self.num_points = 2
        self.points = vtk.vtkPoints()
        self.points.InsertNextPoint(pos)
        self.lines = vtk.vtkCellArray()
        self.lines.InsertNextCell(1)  # number of points
        self.lines.InsertCellPoint(0)
        self.lines_poligon_data = vtk.vtkPolyData()
        self.lines_poligon_data.SetPoints(self.points)
        self.lines_poligon_data.SetLines(self.lines)
        colors = vtk.vtkNamedColors()
        self.GetProperty().SetColor(colors.GetColor3d(self.color))
        self.GetProperty().SetLineWidth(2.5)
        self.GetProperty().SetOpacity(0.5)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.lines_poligon_data)
        mapper.Update()
        self.SetMapper(mapper)

    def update(self):
        tool2world = vtk.vtkTransform()
        tool2world.Concatenate(self.tooltip.current_matrix)
        #world2tool = vtk.vtkTransform()
        #world2tool = tool2world.GetInverse()
        work2world = vtk.vtkTransform()
        work2world.Concatenate(self.work.current_matrix)
        world2work = vtk.vtkTransform()
        world2work = work2world.GetInverse()
        plot_transform = vtk.vtkTransform()
        plot_transform.Concatenate(world2work) # work position [0,0,0] > World [0,y,0]
        plot_transform.Concatenate(tool2world) # tool position [0,0,0] > World [x,0,z]
        x = plot_transform.GetMatrix().GetElement(0,3)
        y = plot_transform.GetMatrix().GetElement(1,3)
        z = plot_transform.GetMatrix().GetElement(2,3)
        current_position = (x, y, z)
        clear = self.comp[self.clear] if self.comp else hal.get_value(self.clear)
        if self.comp[self.clear] or self.initial_run:
            self.points.Reset()
            self.setup_points([x,y,z])
            self.initial_run = False
        self.index += 1
        self.points.InsertNextPoint(current_position)
        self.points.Modified()
        self.lines.InsertNextCell(self.num_points)
        self.lines.InsertCellPoint(self.index - 1)
        self.lines.InsertCellPoint(self.index)
        self.lines.Modified()
        self.SetUserTransform(vtk.vtkTransform())
        self.GetUserTransform().Concatenate(work2world)
        # Calculate tool2work to experiment
        self.tool2work = vtk.vtkTransform()
        self.tool2work.Concatenate(world2work)
        self.tool2work.Concatenate(tool2world) # tool position [0,0,0] > Work [x,y,z]



# create (invisible) actor that can be used to track combined transformation to world coordinates
class Capture(vtk.vtkActor):
    def __init__(self):
        self.SetUserTransform(vtk.vtkTransform())
        self.current_matrix = self.GetMatrix()
        self.tracked_parts = [self]

    def update(self):
        self.current_matrix = self.GetMatrix()            # store the total transformation from this cycle
        self.SetUserTransform(vtk.vtkTransform()) # reset tranform for next update cycle


# draw a line from point_1 to point_2, thickness is optional (defaults to 5)
class HalLine(vtk.vtkActor):
    def __init__(self, comp, x1var, y1var, z1var, x2var, y2var, z2var, stretch=1, r=1):
        self.SetUserTransform(vtk.vtkTransform())
        self.comp = comp
        self.x1var = x1var
        self.y1var = y1var
        self.z1var = z1var
        self.x2var = x2var
        self.y2var = y2var
        self.z2var = z2var
        self.stretch = stretch
        self.r = r
        # we initialize with both points (0,0,0) and get the actual values in the update loop
        (x1,y1,z1) = (x2,y2,z2) = (0,0,0)
        # create the line
        self.lineSource = vtk.vtkLineSource()
        self.lineSource.SetPoint1(x1,y1,z1)
        self.lineSource.SetPoint2(x2,y2,z2)
        # Visualize
        colors = vtk.vtkNamedColors()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(self.lineSource.GetOutputPort())
        self.SetMapper(mapper)
        self.GetProperty().SetLineWidth(r)
        self.GetProperty().SetColor(colors.GetColor3d("Silver"))

    def update(self):
        # update the two points, P0 and P1
        for v in ['x1var','y1var','z1var','x2var','y2var','z2var']:
            # create variable from list and update from class variables of the same name
            globals()[v] = update_passed_args(self, v)
        s = self.stretch
        self.lineSource.SetPoint1(x1var,y1var,z1var)
        self.lineSource.SetPoint2(s*x2var,s*y2var,s*z2var)


# Create a trihedron indicating coordinate orientation
class Axes(vtk.vtkAxesActor):
    def __init__(self, scale=50, radius_factor=0.5, label_x="X", label_y="Y", label_z="Z",):
        self.SetUserTransform(vtk.vtkTransform())
        self.SetXAxisLabelText(label_x) #FIXME Labels are not shown
        self.SetYAxisLabelText(label_y)
        self.SetZAxisLabelText(label_z)
        self.SetShaftTypeToCylinder()
        self.SetCylinderRadius(radius_factor*self.GetCylinderRadius())
        self.SetConeRadius(radius_factor*self.GetConeRadius())
        self.SetScale(scale,scale,scale)


# draw a grid defined by it's normal vector(zx,zy,zz) and x-direction vector(xx, xy, xz)
# optional s to define the half-width from the origin (ox,oy,oz)
class HalGridFromNormalAndDirection(vtk.vtkAssembly):
    def __init__(self, comp, origin, vector_x, vector_z, s=100, r=2, scale_origin=(1,1,1)):
        self.SetUserTransform(vtk.vtkTransform())
        self.comp = comp
        self.ox = origin[0]
        self.oy = origin[1]
        self.oz = origin[2]
        self.xx = vector_x[0]
        self.xy = vector_x[1]
        self.xz = vector_x[2]
        self.zx = vector_z[0]
        self.zy = vector_z[1]
        self.zz = vector_z[2]
        self.s = s
        self.r = r
        self.s_ox  = scale_origin[0]
        self.s_oy  = scale_origin[1]
        self.s_oz  = scale_origin[2]
        self.grid()

    def cross(self, a, b):
        return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

    def grid(self):
        s = self.s
        for i in range (-s,s+10,10):
            # create line in X direction
            line = vtk.vtkActor()
            line.SetUserTransform(vtk.vtkTransform())
            lineSource = vtk.vtkLineSource()
            lineSource.SetPoint1(-s, i, 0)
            lineSource.SetPoint2( s, i, 0)
            # Visualize line in X direction
            colors = vtk.vtkNamedColors()
            mapper1 = vtk.vtkPolyDataMapper()
            mapper1.SetInputConnection(lineSource.GetOutputPort())
            line.SetMapper(mapper1)
            line.GetProperty().SetLineWidth(self.r)
            line.GetProperty().SetColor(colors.GetColor3d("Silver"))
            self.AddPart(line)
            # create line in Y direction
            line = vtk.vtkActor()
            lineSource = vtk.vtkLineSource()
            lineSource.SetPoint1( i,-s, 0)
            lineSource.SetPoint2( i, s, 0)
            # Visualize line in Y direction
            colors = vtk.vtkNamedColors()
            mapper2 = vtk.vtkPolyDataMapper()
            mapper2.SetInputConnection(lineSource.GetOutputPort())
            line.SetMapper(mapper2)
            line.GetProperty().SetLineWidth(self.r)
            line.GetProperty().SetColor(colors.GetColor3d("Silver"))
            self.AddPart(line)

    def update(self):
        # create variables and get values (either 'self.var', 'self.comp[self.var]' or 'hal.get_value(var)')
        for v in ['ox','oy','oz','xx','xy','xz','zx','zy','zz','s_ox','s_oy','s_oz',]:
            # create variable from list and update from class variables of the same name
            globals()[v] = update_passed_args(self, v)
        s = self.s
        vo = [ox, oy, oz]
        vx = [xx, xy, xz]
        vz = [zx, zy, zz]
        # calculate the missing y vector
        vy = [yx, yy, yz] = self.cross(vz,vx)
        matrix = [[ xx, yx, zx, ox*s_ox],
                  [ xy, yy, zy, oy*s_oy],
                  [ xz, yz, zz, oz*s_oz],
                  [  0,  0,  0,       1]]
        transform_matrix = vtk.vtkMatrix4x4()
        for column in range (0,4):
            for row in range (0,4):
                transform_matrix.SetElement(column, row, matrix[column][row])
        self.SetUserMatrix(transform_matrix)


# draw a coordinate system defined by it's normal vector(zx,zy,zz) and x-direction vector(xx, xy, xz)
# optional r to define the thickness of the cylinders
class HalCoordsFromNormalAndDirection(vtk.vtkAxesActor):
    def __init__(self, comp, origin, vector_x, vector_z,
                 scale=100, radius_factor=0.5, scale_origin=(1,1,1), label_prefix=""):
        self.SetUserTransform(vtk.vtkTransform())
        self.comp = comp
        self.ox = origin[0]
        self.oy = origin[1]
        self.oz = origin[2]
        self.xx = vector_x[0]
        self.xy = vector_x[1]
        self.xz = vector_x[2]
        self.zx = vector_z[0]
        self.zy = vector_z[1]
        self.zz = vector_z[2]
        self.s_ox  = scale_origin[0]
        self.s_oy  = scale_origin[1]
        self.s_oz  = scale_origin[2]
        self.create_axes(scale,radius_factor,label_prefix)

    def cross(self, a, b):
        return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

    def create_axes(self, scale=50, radius_factor=0.5, label_prefix="", label_x="X", label_y="Y", label_z="Z",):
        self.SetUserTransform(vtk.vtkTransform())
        self.SetXAxisLabelText(label_prefix + label_x) #FIXME axis labels not showing
        self.SetYAxisLabelText(label_prefix + label_y)
        self.SetZAxisLabelText(label_prefix + label_z)
        self.SetShaftTypeToCylinder()
        self.SetCylinderRadius(radius_factor*self.GetCylinderRadius())
        self.SetConeRadius(radius_factor*self.GetConeRadius())
        self.SetScale(scale,scale,scale)

    def update(self):
        # create variables and get values (either 'self.var', 'self.comp[self.var]' or 'hal.get_value(var)')
        for v in ['ox','oy','oz','xx','xy','xz','zx','zy','zz','s_ox','s_oy','s_oz']:
            # create variable from list and update from class variables of the same name
            globals()[v] = update_passed_args(self, v)
        s_ox = self.s_ox
        s_oy = self.s_oy
        s_oz = self.s_oz
        vo = [ox, oy, oz]
        vx = [xx, xy, xz]
        vz = [zx, zy, zz]
        # calculate the missing y vector
        vy = [yx, yy, yz] = self.cross(vz,vx)
        matrix = [[ xx, yx, zx, ox*s_ox],
                  [ xy, yy, zy, oy*s_oy],
                  [ xz, yz, zz, oz*s_oz],
                  [  0,  0,  0,       1]]
        transform_matrix = vtk.vtkMatrix4x4()
        for column in range (0,4):
            for row in range (0,4):
                transform_matrix.SetElement(column, row, matrix[column][row])
        self.SetUserMatrix(transform_matrix)


class Collection(vtk.vtkAssembly):
    def __init__(self, parts):
        self.SetUserTransform(vtk.vtkTransform())
        for part in parts:
            self.AddPart(part)
            if hasattr(part, "tracked_parts"):
                if not hasattr(self, "tracked_parts"):
                    self.tracked_parts = []
                self.tracked_parts += part.tracked_parts

    def capture(self):
        if hasattr(self, "tracked_parts"):
            if hasattr(self, "transformation"):
                for tracked_part in self.tracked_parts:
                    tracked_part.GetUserTransform().Concatenate(self.transformation)


class Color(Collection):
    def __init__(self, parts, color, opacity=1):
        super().__init__(parts)
        self.color = color # either name string (eg "red","magenta") or normalized RGB as tuple (eg (0.7,0.7,0.1))
        self.opacity = opacity
        def find_actors(parts):
            colors = vtk.vtkNamedColors()
            if hasattr(parts, "GetParts"):
                for item in parts.GetParts():
                    if isinstance(item, vtk.vtkActor):
                        if not isinstance(self.color, tuple):
                            item.GetProperty().SetColor(colors.GetColor3d(self.color))
                        else:
                            item.GetProperty().SetColor(self.color)
                        item.GetProperty().SetOpacity(self.opacity)
                        #item.GetProperty().EdgeVisibilityOn()
                        #item.GetProperty().SetRepresentationToWireframe()
                    elif isinstance(item, vtk.vtkAssembly):
                        find_actors(item)
            else:
                if not isinstance(self.color, tuple):
                    parts.GetProperty().SetColor(colors.GetColor3d(self.color))
                else:
                    parts.GetProperty().SetColor(self.color)
                parts.GetProperty().SetOpacity(self.opacity)
        for part in self.GetParts():
            find_actors(part)


class Translate(Collection):
    def __init__(self, parts, x, y, z):
        super().__init__(parts)
        self.transformation = vtk.vtkTransform()
        self.transformation.Translate(x,y,z)

    def transform(self):
        self.SetUserTransform(self.transformation)


class HalTranslate(Collection):
    def __init__(self, parts, comp, var, x, y, z):
        super().__init__(parts)
        self.comp = comp
        self.var = var
        self.x = x
        self.y = y
        self.z = z

    def update(self):
        for v in ['var']:
            # create variable from list and update from class variables of the same name
            globals()[v] = update_passed_args(self, v)
        self.transformation = vtk.vtkTransform()
        self.transformation.Translate(var*self.x, var*self.y, var*self.z)

    def transform(self):
        self.SetUserTransform(self.transformation)



# translates an object using a variable translation vector, use scale=-1 to change direction
class HalVectorTranslate(Collection):
    def __init__(self, parts, comp, xvar, yvar, zvar, scale=1):
        super().__init__(parts)
        self.comp = comp
        self.xvar = xvar
        self.yvar = yvar
        self.zvar = zvar
        self.current_scale  = scale

    def update(self):
        for v in ['xvar','yvar','zvar']:
            # create variable from list and update from class variables of the same name
            globals()[v] = update_passed_args(self, v)
        self.transformation = vtk.vtkTransform()
        self.transformation.Translate(self.current_scale*xvar, self.current_scale*yvar, self.current_scale*zvar)

    def transform(self):
        self.SetUserTransform(self.transformation)


class Rotate(Collection):
    def __init__(self, parts, th, x, y, z):
        super().__init__(parts)
        self.transformation = vtk.vtkTransform()
        self.transformation.PreMultiply()
        self.transformation.RotateWXYZ(th,x,y,z)

    def transform(self):
        self.SetUserTransform(self.transformation)


class HalRotate(Collection):
    def __init__(self, parts, comp, var, scale, x, y, z):
        super().__init__(parts)
        self.SetUserTransform(vtk.vtkTransform())
        self.comp = comp
        self.var = var
        self.current_scale = scale
        self.x = x
        self.y = y
        self.z = z

    def update(self):
        for v in ['var']:
            # create variable from list and update from class variables of the same name
            globals()[v] = update_passed_args(self, v)
        th = self.current_scale * var
        self.transformation  = vtk.vtkTransform()
        self.transformation.PreMultiply()
        self.transformation.RotateWXYZ(th,self.x, self.y, self.z)

    def transform(self):
        self.SetUserTransform(self.transformation)


class HalEulerRotate(Collection):
    def __init__(self, part, comp, th1, th2, th3, order=123):
        self.AddPart(part)
        self.SetUserTransform(vtk.vtkTransform())
        self.comp = comp
        self.order = order
        self.th1 = th1
        self.th2 = th2
        self.th3 = th3

    def update(self):
        for v in ['order','th1','th2','th3']:
            # create variable from list and update from class variables of the same name
            globals()[v] = update_passed_args(self, v)
        try:
            order = str(int(self.comp[self.order]))
        except:
            order = str(int(self.order))
        if order == '131':
            rotation1 = (th1, 1, 0, 0)
            rotation2 = (th2, 0, 0, 1)
            rotation3 = (th3, 1, 0, 0)
        elif order =='121':
            rotation1 = (th1, 1, 0, 0)
            rotation2 = (th2, 0, 1, 0)
            rotation3 = (th3, 1, 0, 0)
        elif order =='212':
            rotation1 = (th1, 0, 1, 0)
            rotation2 = (th2, 1, 0, 0)
            rotation3 = (th3, 0, 1, 0)
        elif order =='232':
            rotation1 = (th1, 0, 1, 0)
            rotation2 = (th2, 0, 0, 1)
            rotation3 = (th3, 0, 1, 0)
        elif order =='323':
            rotation1 = (th1, 0, 0, 1)
            rotation2 = (th2, 0, 1, 0)
            rotation3 = (th3, 0, 0, 1)
        elif order =='313':
            rotation1 = (th1, 0, 0, 1)
            rotation2 = (th2, 1, 0, 0)
            rotation3 = (th3, 0, 0, 1)
        elif order =='123':
            rotation1 = (th1, 1, 0, 0)
            rotation2 = (th2, 0, 1, 0)
            rotation3 = (th3, 0, 0, 1)
        elif order =='132':
            rotation1 = (th1, 1, 0, 0)
            rotation2 = (th2, 0, 0, 1)
            rotation3 = (th3, 0, 1, 0)
        elif order =='213':
            rotation1 = (th1, 0, 1, 0)
            rotation2 = (th2, 1, 0, 0)
            rotation3 = (th3, 0, 0, 1)
        elif order =='231':
            rotation1 = (th1, 0, 1, 0)
            rotation2 = (th2, 0, 0, 1)
            rotation3 = (th3, 1, 0, 0)
        elif order =='321':
            rotation1 = (th1, 0, 0, 1)
            rotation2 = (th2, 0, 1, 0)
            rotation3 = (th3, 1, 0, 0)
        elif order =='312':
            rotation1 = (th1, 0, 0, 1)
            rotation2 = (th2, 1, 0, 0)
            rotation3 = (th3, 0, 1, 0)
        euler_transform = vtk.vtkTransform()
        euler_transform.RotateWXYZ(*rotation1)
        euler_transform.RotateWXYZ(*rotation2)
        euler_transform.RotateWXYZ(*rotation3)
        self.SetUserMatrix(euler_transform.GetMatrix())


# shows an object if const=var and hides it otherwise, behavior can be changed
# using the optional arguments for scalefactors when true or false
class HalShow(Collection):
    def __init__(self, parts, comp, const, var, scaleby_true=1, scaleby_false=0):
        super().__init__(parts)
        self.comp = comp
        if isinstance(const, list):
             self.const = const
        else:
             self.const = [const]
        self.var = var
        self.current_scaleby_true = scaleby_true
        self.current_scaleby_false = scaleby_false

    def update(self):
        s_t = self.current_scaleby_true
        s_f = self.current_scaleby_false
        for v in ['var']:
            # create variable from list and update from class variables of the same name
            globals()[v] = update_passed_args(self, v)
        if var in self.const:
            self.SetScale(s_t,s_t,s_t)
        else:
            self.SetScale(s_f,s_f,s_f)


# Create a text overlay (HUD)
# color can be either name string (eg "red","magenta") or normalized RGB as tuple (eg (0.7,0.7,0.1))
class Hud(vtk.vtkActor2D):
    def __init__(self, comp, var, const, color, opacity=1, font_size=20, line_spacing=1):
        self.comp = comp
        self.var = var
        self.const = const
        self.strs = []
        self.hud_lines = []
        self.show_tags = []
        self.hide_alls = []
        self.extra_text_enable = False
        self.extra_text = None
        self.textMapper = vtk.vtkTextMapper()
        tprop = self.textMapper.GetTextProperty()
        tprop.SetLineSpacing(line_spacing)
        tprop.SetFontSize(font_size)
        tprop.SetFontFamilyAsString("Courier")
        tprop.SetJustificationToLeft()
        tprop.SetVerticalJustificationToTop()
        colors = vtk.vtkNamedColors()
        if not isinstance(color, tuple):
            tprop.SetColor(colors.GetColor3d(color))
        else:
            tprop.SetColor(color)
        tprop.SetOpacity(opacity)
        self.SetMapper(self.textMapper)
        self.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
        self.GetPositionCoordinate().SetValue(0.05, 0.95)

    # displays a string, optionally a tag or list of tags can be assigned
    def add_txt(self, string, tag=None):
        self.hud_lines += [[str(string), None, tag]]

    # displays a formatted pin value (can be embedded in a string)
    def add_pin(self, string, pin=None, tag=None):
        self.hud_lines += [[str(string), pin, tag]]

    # shows all lines with the specified tags if the pin value = val
    def show_tags_if_pin_eq_val(self, tag, pin, val=True):
        self.show_tags += [[tag, pin, val]]

    # shows all lines with a tag equal to the pin value + offset
    def show_tag_eq_pin_offs(self, pin, offs=0):
        self.show_tags += [[pin, None, offs]]

    # hides the complete hud if the pin value is equal to val
    def hide_all(self, pin, val=True):
        self.hide_alls += [[pin, val]]

    # update the lines in the hud using the lists created above
    def update(self):
        for v in ['var']:
            # create variable from list and update from class variables of the same name
            globals()[v] = update_passed_args(self, v)
        hide_hud = 0 if var == self.const else 1

        strs = []
        show_list = [None]
        # check if hud should be hidden
        for a in self.hide_alls:
            if hal.get_value(a[0]) == a[1]:
                hide_hud = 1
        if hide_hud == 0:
            # create list of all line tags to be shown
            for b in self.show_tags:
                if b[1] == None: # show_tag_eq_pin_offs
                    tag = int(hal.get_value(b[0]) + b[2])
                else: # show_tags_if_pin_eq_val
                    if  hal.get_value(b[1]) == b[2]:
                        tag = b[0]
                if not isinstance(tag, list):
                    tag = [tag]
                show_list = show_list + tag
            # build the strings
            for c in self.hud_lines:
                if not isinstance(c[2], list):
                    c[2] = [c[2]]
                if any(item in c[2] for item in show_list):
                    if c[1] == None: # txt
                        strs += [c[0]]
                    else: # pin
                        strs += [c[0].format(hal.get_value(c[1]))]
        combined_string = ""
        for string in strs:
            combined_string += (string + "\n")
        if self.extra_text_enable and self.extra_text and not hide_hud:
            combined_string += self.extra_text
        self.textMapper.SetInput(combined_string)






def main(comp,
         model, tooltip, work, huds,
         window_title="Vtk-Vismach", window_width=600, window_height=300,
         camera_azimuth=-50, camera_elevation=30,
         background_rgb = (0.2, 0.3, 0.4)):

    # create a separate hal component and create a pin to clear the backplot
    vcomp = hal.component("vismach")
    vcomp.newpin("plotclear",hal.HAL_BIT,hal.HAL_IN)
    vcomp.ready()
    # create the backplot to be added to the renderer
    backplot = Plotter(vcomp, work, tooltip, "plotclear")
    # Event loop to periodically update the model
    def update():
        def get_actors_to_update(objects):
            for item in objects.GetParts():
                if hasattr(item, "update"):
                    item.update()
                if hasattr(item, "transform"):
                    item.transform()
                if hasattr(item, "capture"):
                    item.capture()
                if isinstance(item, vtk.vtkAssembly):
                    get_actors_to_update(item)
        get_actors_to_update(model)
        backplot.update()
        t2w = backplot.tool2work.GetMatrix()
        r1=("{:8.3f} {:8.3f} {:8.3f} {:8.3f}".format(t2w.GetElement(0,0), t2w.GetElement(0,1), t2w.GetElement(0,2), t2w.GetElement(0,3)))
        r2=("{:8.3f} {:8.3f} {:8.3f} {:8.3f}".format(t2w.GetElement(1,0), t2w.GetElement(1,1), t2w.GetElement(1,2), t2w.GetElement(1,3)))
        r3=("{:8.3f} {:8.3f} {:8.3f} {:8.3f}".format(t2w.GetElement(2,0), t2w.GetElement(2,1), t2w.GetElement(2,2), t2w.GetElement(2,3)))
        for hud in huds:
            hud.extra_text = "\ntool2work Matrix:"+"\n"+r1+"\n"+r2+"\n"+r3
            hud.update()
        vtkWidget.GetRenderWindow().Render()

    app = Qt.QApplication([])
    # Qt Window
    mainWindow = Qt.QMainWindow()
    mainWindow.setMinimumHeight(window_height)
    mainWindow.setMinimumWidth(window_width)
    mainWindow.setWindowTitle(window_title)
    # A renderer and render window
    renderer = vtk.vtkRenderer()
    renderer.AddActor(model)
    renderer.AddActor(backplot)
    if huds:
        if not isinstance(huds, list):
            huds = [huds]
        for hud in huds:
            renderer.AddActor(hud)
    renderer.SetBackground(*background_rgb)
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    # An interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderWindow)
    # A widget to indicate the orientation of the model in the window
    def trihedron():
        axes = vtk.vtkAxesActor()
        axes.SetShaftTypeToCylinder()
        axes.SetXAxisLabelText("X")
        axes.SetYAxisLabelText("Y")
        axes.SetZAxisLabelText("Z")
        axes.SetCylinderRadius(0.5 * axes.GetCylinderRadius())
        return axes
    orientation_marker = vtk.vtkOrientationMarkerWidget()
    orientation_marker.SetOrientationMarker(trihedron())
    #orientation_marker.InteractiveOn()
    # Camera setup
    camera = renderer.GetActiveCamera()
    # We want  Z pointing up and the X pointing right
    camera_fp = camera.GetFocalPoint()
    camera_t = vtk.vtkTransform()
    camera_t.Translate(camera_fp[0], camera_fp[1], camera_fp[2])
    camera_t.RotateX(90)
    camera_t.RotateY(90)
    camera_t.Translate(-camera_fp[0], -camera_fp[1], -camera_fp[2])
    camera.ApplyTransform(camera_t)
    camera.Azimuth(camera_azimuth)
    camera.Elevation(camera_elevation)
    renderer.ResetCamera()
    # Put everything in the Qt window
    vtkWidget = QVTKRenderWindowInteractor(mainWindow)
    vtkWidget.GetRenderWindow().AddRenderer(renderer)
    mainWindow.setCentralWidget(vtkWidget)
    # Set interactor style and initialize
    interactor = vtkWidget.GetRenderWindow().GetInteractor()
    orientation_marker.SetInteractor(interactor)
    orientation_marker.EnabledOn()
    interactor_style = vtk.vtkInteractorStyleTrackballCamera()
    interactor.SetInteractorStyle(interactor_style)
    interactor.Initialize()
    # NOTE
    # We really only use Qt because we need a timer outside of VTK. Due to a vtk bug we cannot use
    # the vtk timer as it stops reporting when we interact with the window (eg rotating the scene)
    # Set up a Qt timer to create update events
    timer = QTimer()
    timer.timeout.connect(update)
    timer.start(100)
    # Show Qt window
    mainWindow.show()
    app.exec_()
