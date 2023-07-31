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

import ase.io
from batoms import Batoms

def main():
    parser = argparse.ArgumentParser(
        prog='Beautiful Atoms Renderer',
        description='Render images with BAtoms')
    
    parser.add_argument("--cube-file", help = "Path to the cube file to read")
    parser.add_argument("--isovalues", help = "List of isovalues to render", nargs = "*", type = float, default = [])
    
    # Both blender and python share the same command line arguments.
    # They are separated by double dash ('--'), everything before is for blender,
    # everything afterwards is for python (except for the first argument, wich is
    # the program name, which is for both).
    if "--" in sys.argv:
        python_argv = sys.argv[sys.argv.index("--") +1:]
    else:
        python_argv = sys.argv
    
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
        [({'level': iso_value, 'color': [1, 1, 0, 0.5]}, {'level': - iso_value, 'color': [0, 0, 0.8, 0.5]})
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
    
    mol.get_image(viewport = [1, 0, 0], output = 'mol.1.png')
    mol.get_image(viewport = [0, 1, 0], output = 'mol.2.png')
    mol.get_image(viewport = [0, 0, 1], output = 'mol.3.png')
    
    return 0
    
# If we've been invoked as a program, call main().    
if __name__ == '__main__':
    main()
    #sys.exit(main())