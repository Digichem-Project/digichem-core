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
        ):
    """
    """
    rotations = [] if rotations is None else rotations
    isovalues = [] if isovalues is None else isovalues
    
    surface_settings = []
    for isovalue in isovalues:
        if isotype in ["both", "positive"]:
            surface_settings.append({'level': isovalue, 'color': primary_color})
        
        if isotype in ["both", "negative"]:
            surface_settings.append({'level': -isovalue, 'color': secondary_color})
    
    
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
    parser.add_argument("--primary-color", help = "RGBA for one of the colors to use for isosurfaces", nargs = 4, default = [1, 0.058, 0.0, 0.55])
    parser.add_argument("--secondary-color", help = "RGBA for the other color to use for isosurfaces", nargs = 4, default = [0.1, 0.1, 0.9, 0.55])
    parser.add_argument("--cpus", help = "Number of parallel CPUs to use for rendering", type = int, default = 1)
    parser.add_argument("--orientation", help = "The orientation to render from, as x, y, z values", nargs = 3, type = float, default = [0, 0, 1])
    parser.add_argument("--resolution", help = "The output resolution in px", type = int, default = 1024)
    parser.add_argument("--render-samples", help = "The maximum number of render samples, more generally results in higher quality but longer render times", type = int, default = 256)
    parser.add_argument("--rotations", help = "A list of rotations (in JSON) to rotate the molecule to a given alignment. The first item in each list item is the axis to rotate about (0=x, 1=y, 2=z), the second is the angle to rotate by (in radians)", default = None)
    
    # Both blender and python share the same command line arguments.
    # They are separated by double dash ('--'), everything before is for blender,
    # everything afterwards is for python (except for the first argument, wich is
    # the program name, which is for both).
    if "--" in sys.argv:
        python_argv = sys.argv[sys.argv.index("--") +1:]
    else:
        python_argv = []
    
    args = parser.parse_args(python_argv)
    
    # Parse rotations.
    if args.rotations is not None:
        rotations = yaml.safe_load(args.rotations)
    
    else:
        rotations = None
    
    # Load the input data.
    mol = add_molecule(
        args.cube_file,
        name = "molecule",
        visible = True,
        rotations = rotations,
        isovalues = args.isovalues,
        isotype = args.isotype,
        primary_color = args.primary_color,
        secondary_color = args.secondary_color if args.isocolor == "sign" else args.primary_color
    )
        
    # If we have a second cube, add that too.
    if args.second_cube is not None:
        mol2 = add_molecule(
            args.second_cube,
            name = "molecule2",
            visible = False,
            rotations = rotations,
            isovalues = args.isovalues,
            isotype = args.isotype,
            primary_color = args.primary_color if args.isocolor == "sign" else args.secondary_color,
            secondary_color = args.secondary_color
        )
        
    
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
    
    # Change light intensity.
    bpy.data.lights["batoms_light_Default"].node_tree.nodes["Emission"].inputs[1].default_value = 0.3
    #mol.render.lights["Default"].energy=10
    
    # Enable to add an outline.
    #bpy.context.scene.render.use_freestyle = True
    
    # Performance options.
    bpy.context.scene.render.threads_mode = 'FIXED'
    bpy.context.scene.render.threads = args.cpus
    
    mol.get_image(viewport = args.orientation, output = args.output)
    
    return 0
    
# If we've been invoked as a program, call main().    
if __name__ == '__main__':
    main()
    #sys.exit(main())