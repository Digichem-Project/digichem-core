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

import ase.io
from batoms import Batoms

def main():
    parser = argparse.ArgumentParser(
        prog='Beautiful Atoms Renderer',
        description='Render images with BAtoms')
    
    parser.add_argument("cube_file", help = "Path to the cube file to read")
    parser.add_argument("output", help = "File to write to")
    parser.add_argument("--isovalues", help = "List of isovalues to render", nargs = "*", type = float, default = [])
    parser.add_argument("--cpus", help = "Number of parallel CPUs to use for rendering", type = int, default = 1)
    parser.add_argument("--orientation", help = "The orientation to render from, as x, y, z values", nargs = 3, type = float, default = [0, 0, 1])
    
    # Both blender and python share the same command line arguments.
    # They are separated by double dash ('--'), everything before is for blender,
    # everything afterwards is for python (except for the first argument, wich is
    # the program name, which is for both).
    if "--" in sys.argv:
        python_argv = sys.argv[sys.argv.index("--") +1:]
    else:
        python_argv = []
    
    args = parser.parse_args(python_argv)
    
    # Load the input data.
    cube = ase.io.read(args.cube_file, format="cube", read_data=True, full_output=True)
    cube["atoms"].translate(-cube["origin"][0:3])
    mol = Batoms("molecule", from_ase = cube["atoms"])
    
    mol.render.lights["Default"].energy=10
    
    # Change molecule style.
    mol.model_style = 1
    
    # Hide cell boundaries.
    mol.cell.hide = True
    
    # Add volumes.
    mol.volumetric_data['orbital'] = cube['data']
    
    surfaces = itertools.chain(*
        #[({'level': iso_value, 'color': [1, 0.058, 0, 0.55]}, {'level': - iso_value, 'color': [0, 0, 0.8, 0.55]})
        [({'level': -iso_value, 'color': [1, 0.098, 0.04, 0.55]}, {'level': iso_value, 'color': [0.1, 0.1, 0.8, 0.55]})
        for iso_value in args.isovalues
        ]
    )
    
    for index, settings in enumerate(surfaces):
        mol.isosurface.settings[index+1] = settings
    
    mol.isosurface.draw()
    
    # Setup rendering settings.
#     mol.render.engine = 'workbench'
#     mol.render.engine = 'eevee'
    mol.render.engine = 'cycles'
    # Set up cycles for good quality rendering.
    # Prevents early end to rendering (forces us to use the actual number of samples).
    bpy.context.scene.cycles.use_adaptive_sampling = False
    # Quality control, more = better and slower.
    bpy.context.scene.cycles.samples = 256
    # Post-processing to remove noise, works well for coloured backgrounds, useless for transparency.
    bpy.context.scene.cycles.use_denoising = True
    # Ray-tracing options
    bpy.context.scene.cycles.max_bounces = 48
    bpy.context.scene.cycles.transparent_max_bounces = 24
    
    # Change light intensity.
    bpy.data.lights["batoms_light_Default"].node_tree.nodes["Emission"].inputs[1].default_value = 0.3
    
    # Enable to add an outline.
    #bpy.context.scene.render.use_freestyle = True
    
    # Performance options.
    bpy.context.scene.render.threads_mode = 'FIXED'
    bpy.context.scene.render.threads = args.cpus
    
    # Colour tuning.
    # Carbon to 'black'.
    mol["C"].color = (0.095, 0.095, 0.095, 1)
    mol["B"].color = (1.0, 0.396, 0.468, 1)


    
    mol.get_image(viewport = args.orientation, output = args.output)
    
    return 0
    
# If we've been invoked as a program, call main().    
if __name__ == '__main__':
    main()
    #sys.exit(main())