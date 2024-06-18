# Render something with Beautiful Atoms.
# To use this script, run something like:
#
# # To disabled possibly conflicting packages outside of the conda env.
# PYTHONNOUSERSITE=1
# blender -b -P batoms-renderer.py
#
# Where 'blender' is the path to the Beautiful Atoms Blender executable.

import sys
import argparse
import itertools
import bpy
import yaml
import math
from pathlib import Path
import logging

import ase.io
from batoms import Batoms
from batoms.utils.butils import object_mode

def add_molecule(
        cube_file,
        name,
        rotations = None,
        visible = True,
        isovalues = None,
        isotype = "both",
        primary_color = [1, 0.058, 0.0, 0.55],
        secondary_color = [0.1, 0.1, 0.9, 0.55],
        style = "default"
        ):
    """
    """
    rotations = [] if rotations is None else rotations
    isovalues = [] if isovalues is None else isovalues
    
    surface_settings = []
    for isovalue in isovalues:
        if isotype in ["both", "positive"]:
            surface_settings.append({'level': isovalue, 'color': primary_color, 'material_style': style})
        
        if isotype in ["both", "negative"]:
            surface_settings.append({'level': -isovalue, 'color': secondary_color, 'material_style': style})
    
    
    # Load the input data.
    cube = ase.io.read(cube_file, format="cube", read_data=True, full_output=True)
    
    # The centre of the cube is often offset, fix that by shifting the atoms.
    cube["atoms"].translate(-cube["origin"][0:3])
    
    # Get the mol object.
    mol = Batoms(name, from_ase = cube["atoms"])
    
    # Set some look and feel options.    
    # Change molecule style.
    mol.model_style = 1
    
    # Hide cell boundaries.
    mol.cell.hide = True
    
    # Colour tuning.
    # Carbon to 'black'.
    try:
        mol["C"].color = (0.095, 0.095, 0.095, 1)
    except AttributeError:
        pass
    try:
        mol["B"].color = (1.0, 0.396, 0.468, 1)
    except AttributeError:
        pass
    
    if not visible:
        mol.hide = True
    
    # Add volumes.
    if len(surface_settings) != 0:
        mol.volumetric_data['surface'] = cube['data']
        
        for index, settings in enumerate(surface_settings):
            mol.isosurface.settings[index+1] = settings
        
        mol.isosurface.draw()
    
    # Now move the entire molecule (isosurface and all) back to the origin.
    mol.translate(cube["origin"][0:3])
    
    # Fix the origin point so we can still rotate properly.
    # For some reason, this code moves the bond objects to a new location?
#     object_mode()
#     bpy.ops.object.select_all(action='DESELECT')
#     mol.obj.select_set(True)
#     bpy.context.view_layer.objects.active = mol.obj
#     bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
    
    # If we have any rotations, apply those.
    for axis, angle in rotations:
        # Convert to degree.
        # Not used.
        degree = angle * (180/math.pi)
        
        if axis == 0:
            axis = "x"
        
        elif axis == 1:
            axis = "y"
        
        elif axis == 2:
            axis = "z"
        
        else:
            raise ValueError("Unknown rotation axis '{}'".format(axis))
        
        #mol.rotate(degree, axis)
        
        # Taken from batoms source (collection.py: 51)
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        mol.obj.select_set(True)
        bpy.context.view_layer.objects.active = mol.obj
        bpy.ops.transform.rotate(value=angle, orient_axis=axis.upper(),
                                 center_override = (0,0,0))
    
    return mol

# Adapted from https://blender.stackexchange.com/questions/5898/how-can-i-create-a-cylinder-linking-two-points-with-python?newreg=f372ba9448694f5b97879a6dab963cee
def draw_primitive(start, end, radius, mesh_type, color, collection = None):
    """
    """
    dx = end[0] - start[0]
    dy = end[1] - start[1]
    dz = end[2] - start[2]  
    dist = math.sqrt(dx**2 + dy**2 + dz**2)
    
    # First draw the shape.
    if mesh_type == "cylinder":
        bpy.ops.mesh.primitive_cylinder_add(
            radius = radius, 
            depth = dist,
            vertices = 32,
            location = (dx/2 + start[0], dy/2 + start[1], dz/2 + start[2]) 
        )
    elif mesh_type == "cone":
        bpy.ops.mesh.primitive_cone_add(
            radius1 = radius,
            depth = dist,
            vertices = 32,
            location = (dx/2 + start[0], dy/2 + start[1], dz/2 + start[2])
        )
    else:
        raise ValueError("Unknown mesh type '{}'".format(mesh_type))
        
    # Get a reference to the object we just made.
    obj = bpy.context.active_object
    
    phi = math.atan2(dy, dx)
    theta = math.acos(dz/dist) 
    
    bpy.context.object.rotation_euler[1] = theta 
    bpy.context.object.rotation_euler[2] = phi
    
    # Get material
    mat = bpy.data.materials.new(name="Material")
    
    # assign to 1st material slot
    obj.active_material = mat
        
    mat.use_nodes = True
    tree = mat.node_tree
    nodes = tree.nodes
    bsdf = nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = color
    bsdf.inputs["Metallic"].default_value = 0.1
    bsdf.inputs["Specular"].default_value = 0.2
    bsdf.inputs["Roughness"].default_value = 0.2
    mat.diffuse_color = color
    
    # If we've been asked to, asign our new object to a given collection.
    # First unlink from any old collections
    for coll in obj.users_collection:
        # Unlink the object
        coll.objects.unlink(obj)

    # Link each object to the target collection
    collection.objects.link(obj)

    
def draw_arrow(start, end, radius, color, split = 0.9, collection = None):
    # Decide what proportion of the total vector to dedicate to the arrow stem and head.
    dx = end[0] - start[0]
    dy = end[1] - start[1]
    dz = end[2] - start[2]
    dist = math.sqrt(dx**2 + dy**2 + dz**2)
    
    join = (dx * split + start[0], dy * split + start[1], dz * split + start[2])
    draw_primitive(start, join, radius, "cylinder", color, collection = collection)
    draw_primitive(join, end, radius*2, "cone", color, collection = collection)


def main():
    parser = argparse.ArgumentParser(
        prog='Beautiful Atoms Renderer',
        description='Render images with BAtoms')
    
    parser.add_argument("cube_file", help = "Path to the cube file to read")
    parser.add_argument("output", help = "File to write to")
    parser.add_argument("--second_cube", help = "Optional second cube file to read additional isosurface data from", default = None)
    parser.add_argument("--isovalues", help = "List of isovalues to render", nargs = "*", type = float, default = [])
    parser.add_argument("--isotype", help = "Whether to render positive, negative or both isosurfaces for each isovalue", choices = ["positive", "negative", "both"], default = "both")
    parser.add_argument("--isocolor", help = "The colouring method to use for isosurfaces", choices = ["sign", "cube"], default = "sign")
    parser.add_argument("--primary-color", help = "RGBA for one of the colors to use for isosurfaces", type = float, nargs = 4, default = [0.1, 0.1, 0.9, 0.65])
    parser.add_argument("--secondary-color", help = "RGBA for the other color to use for isosurfaces", type = float, nargs = 4, default = [1, 0.058, 0.0, 0.65])
    parser.add_argument("--style", help = "Material style for isosurfaces", choices = ('default', 'metallic', 'plastic', 'ceramic', 'mirror'), default = "default")
    parser.add_argument("--cpus", help = "Number of parallel CPUs to use for rendering", type = int, default = 1)
    parser.add_argument("--use-gpu", help = "Whether to enable GPU rendering", action = "store_true")
    parser.add_argument("--orientation", help = "The orientation to render from, as x, y, z values", nargs = 3, type = float, default = [0, 0, 1])
    parser.add_argument("--resolution", help = "The output resolution in px", type = int, default = 1024)
    parser.add_argument("--render-samples", help = "The maximum number of render samples, more generally results in higher quality but longer render times", type = int, default = 256)
    parser.add_argument("--rotations", help = "A list of rotations (in JSON) to rotate the molecule to a given alignment. The first item in each list item is the axis to rotate about (0=x, 1=y, 2=z), the second is the angle to rotate by (in radians)", nargs = "*", default = [])
    parser.add_argument("--dipoles", help = "Draw dipoles from a list of the following data (in JSON): 0) start coord, 1) end coord, 2) RGBA color information", nargs = "*", default = [])
    parser.add_argument("--alpha", help = "Override the opacity value for all molecule objects (but not dipoles) to this value, useful for showing dipole arrows more clearly", default = None, type = float)
    parser.add_argument("--perspective", help = "The perspective mode, either orthographic or perspective", default = "perspective", choices = ["perspective", "orthographic"])
    parser.add_argument("--padding", help = "Padding", type = float, default = 1.0)
    
    # Both blender and python share the same command line arguments.
    # They are separated by double dash ('--'), everything before is for blender,
    # everything afterwards is for python (except for the first argument, wich is
    # the program name, which is for both).
    if "--" in sys.argv:
        python_argv = sys.argv[sys.argv.index("--") +1:]
    else:
        python_argv = []
    
    args = parser.parse_args(python_argv)
    
    # Batoms or blender will silently set the extension to png if it's not already.
    # This is surprising, so stop now before that happens.
    if Path(args.output).suffix.lower() != ".png":
        raise ValueError("Output location must have a .png extension")
    
    if args.rotations is not None:
        rotations = [yaml.safe_load(rotation) for rotation in args.rotations]
    
    # Remove the starting cube object.
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()

    # Load the input data.
    mol = add_molecule(
        args.cube_file,
        name = "molecule",
        visible = True,
        rotations = rotations,
        isovalues = args.isovalues,
        isotype = args.isotype,
        primary_color = args.primary_color,
        secondary_color = args.secondary_color if args.isocolor == "sign" else args.primary_color,
        style = args.style
    )
    
    # Uncomment to show atom labels.
    # Needs some tweaking to appear in render (viewport only by default).
    #mol.show_label = 'species'
    
    # If we have a second cube, add that too.
    if args.second_cube is not None:
        mol2 = add_molecule(
            args.second_cube,
            name = "molecule2",
            visible = True,
            rotations = rotations,
            isovalues = args.isovalues,
            isotype = args.isotype,
            primary_color = args.primary_color if args.isocolor == "sign" else args.secondary_color,
            secondary_color = args.secondary_color,
            style = args.style
        )
    
    if args.alpha:
        # Set all materials transparent.
        for material in bpy.data.materials:
            try:
                material.node_tree.nodes['Principled BSDF'].inputs['Alpha'].default_value = args.alpha
            except Exception as e:
                pass
        
    
    # Draw any dipoles.
    if args.dipoles is not None:
        
        dipoles = [yaml.safe_load(dipole) for dipole in args.dipoles]
        for start_coord, end_coord, rgba in dipoles:
            draw_arrow(start_coord, end_coord, 0.08, rgba, collection = mol.coll)
        
    
    # Setup rendering settings.
#     mol.render.engine = 'workbench'
#     mol.render.engine = 'eevee'
    mol.render.engine = 'cycles'
    mol.render.resolution = [args.resolution, args.resolution]
    # Set up cycles for good quality rendering.
    # Prevents early end to rendering (forces us to use the actual number of samples).
    bpy.context.scene.cycles.use_adaptive_sampling = False
    # Quality control, more = better and slower.
    bpy.context.scene.cycles.samples = args.render_samples
    # Post-processing to remove noise, works well for coloured backgrounds, useless for transparency.
    bpy.context.scene.cycles.use_denoising = True
    # Ray-tracing options
    bpy.context.scene.cycles.max_bounces = 48
    bpy.context.scene.cycles.transparent_max_bounces = 24
    if args.use_gpu:
        bpy.context.scene.cycles.device = "GPU"
    
    # Use maximum compression.
    bpy.context.scene.render.image_settings.compression = 1000
    
    # Change light intensity.
    bpy.data.lights["batoms_light_Default"].node_tree.nodes["Emission"].inputs[1].default_value = 0.5 #0.3
    #mol.render.lights["Default"].energy=10
    
    # Change view mode.
    if args.perspective == "perspective":
        mol.render.camera.type = "PERSP"
    
    else:
        mol.render.camera.type = "ORTHO"
    
    # Enable to add an outline.
    #bpy.context.scene.render.use_freestyle = True
    
    # We have plenty of memory to play with, use one tile.
    bpy.context.scene.cycles.tile_x = args.resolution
    bpy.context.scene.cycles.tile_y = args.resolution
    bpy.context.scene.cycles.tile_size = args.resolution

    # Performance options.
    bpy.context.scene.render.threads_mode = 'FIXED'
    bpy.context.scene.render.threads = args.cpus
    
    mol.get_image(viewport = args.orientation, output = args.output, padding = args.padding)
#     # Move the camera.
#     mol.render.camera.location = (100,0,0)
#     mol.render.camera.look_at = mol.get_center_of_geometry()
#     bpy.ops.object.select_all(action='DESELECT')
#     for obj in mol.coll.objects[:]:
#         obj.select_set(True)
#     #bpy.ops.view3d.camera_to_view_selected()
    
    return 0
    
# If we've been invoked as a program, call main().    
if __name__ == '__main__':
    try:
        sys.exit(main())
    
    except Exception as e:
        logging.error("Erro", exc_info = True)
        sys.exit(1)