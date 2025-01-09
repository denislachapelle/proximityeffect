#DL241009
# two round wires model in python.
# 2d geometric model.
# to be used MFEM proximity2d.cpp to simulate electromagnetic proximity effect.
#
#define the following physical
#contour, conductor.

import gmsh
import sys
import math


gmsh.initialize()

gmsh.model.add("tworoundwires2d_py")

radius = 0.016
spacing = 0.07945

#wire 1
w=0
gmsh.model.occ.addCircle(-spacing/2, 0, 0, radius, 1)
gmsh.model.occ.addCurveLoop([1], 1)
gmsh.model.occ.addPlaneSurface([1], 1)

#wire 2
w=1
gmsh.model.occ.addCircle(spacing/2, 0, 0, radius, 5)
gmsh.model.occ.addCurveLoop([5], 5)
gmsh.model.occ.addPlaneSurface([5], 5)


#domain limit
gmsh.model.occ.addCircle(0, 0, 0, 3.5*(radius+spacing), 6)
gmsh.model.occ.addCurveLoop([6], 6)
gmsh.model.occ.addPlaneSurface([6, -1, -5], 6)

#
#synchronize prior to add physical group.
#
gmsh.model.occ.synchronize()
#
#add physical groups.
#
gmsh.model.addPhysicalGroup(2, [1], 1, name="wire_1")
gmsh.model.addPhysicalGroup(2, [5], 2, name="wire_2")
gmsh.model.addPhysicalGroup(2, [6], 3, name="air")
gmsh.model.addPhysicalGroup(1, [6], 4, name="aircontour")
#gmsh.model.addPhysicalGroup(1, [1], 5, name="wire_1_contour")




# We can then generate a 2D mesh...
gmsh.option.setNumber("Mesh.Algorithm", 6)

gmsh.model.mesh.generate(2)
gmsh.model.mesh.refine()
gmsh.model.mesh.refine()



# glvis can read mesh version 2.2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# ... and save it to disk
gmsh.write("tworoundwires2d.msh")

# start gmsh
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#before leaving.
gmsh.finalize()

