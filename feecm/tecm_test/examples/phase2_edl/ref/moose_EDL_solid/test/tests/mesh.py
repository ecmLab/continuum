import gmsh
import sys

# Initialize the gmsh library
gmsh.initialize()

# Create a new model
gmsh.model.add("1D_Mesh_Varying_Density")

# Parameters
Length = 38.0       # Total length of the line (1000 units)
Lc_left = 1 / 1e7   # Characteristic length at the left side (finer mesh)
Lc_right = 1 / 1e4  # Characteristic length at the right side (coarser mesh)
n_points = 100000     # Number of points for the mesh (adjust for refinement)
progression = (Lc_right / Lc_left) ** (1 / (n_points - 1))  # Geometric progression

# Define points at both ends of the line with explicit tags
p1 = gmsh.model.geo.addPoint(0, 0, 0, Lc_left, tag=1)    # Left side at x = 0
p2 = gmsh.model.geo.addPoint(Length, 0, 0, Lc_right, tag=2)  # Right side at x = 1000

# Create a line between the two points with an explicit tag
line = gmsh.model.geo.addLine(p1, p2, tag=1)

# Synchronize the model before applying physical groups
gmsh.model.geo.synchronize()

# Apply transfinite meshing along the line with progression
gmsh.model.geo.mesh.setTransfiniteCurve(line, n_points)
gmsh.model.geo.mesh.setTransfiniteCurve(line, n_points, coef=progression)

# # Add a physical group for the 1D line (this is the main domain)
# gmsh.model.addPhysicalGroup(1, [line], name="1D_Line")

# # Add physical groups for the boundaries (sidesets in MOOSE)
# # For MOOSE, boundaries must be physical entities of dimension 1 in a 1D mesh
# gmsh.model.addPhysicalGroup(1, [p1], name="left")
# gmsh.model.addPhysicalGroup(1, [p2], name="right")

# Synchronize again after adding physical groups
gmsh.model.geo.synchronize()

# Generate the mesh
gmsh.model.mesh.generate(1)

# Set the element order to 2 (second-order Lagrange)
# gmsh.model.mesh.setOrder(2)

# Optionally save the mesh to a file (e.g., .msh format)
# gmsh.write("line-2e2-L500.msh")
gmsh.write("line-1e5-L38.msh")

# Optionally run the GUI to visualize the mesh
# gmsh.fltk.run()

# Finalize the gmsh library
gmsh.finalize()