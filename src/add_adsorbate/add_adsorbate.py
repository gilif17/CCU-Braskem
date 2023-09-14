"""
Script Name: Add Adsorbate to POSCAR
Author: João V. V. Cassiano
Date: 14-09-2023
Version: 1.0.0
Description: 
    Este script adiciona um adsorvato a uma superfície em todas as posições possíveis de um tipo de site especificado.
    O script aceita cinco argumentos de linha de comando: o caminho para o arquivo POSCAR de entrada, o caminho para o arquivo POSCAR de saída, 
    o símbolo do átomo adsorvato, o tipo de site e a altura acima do site para colocar o adsorvato.
Usage:
    python add_adsorbate.py input_poscar_path output_poscar_path adsorbate_atom site height
"""


import argparse
from ase import Atoms
from ase.io import read, write
from ase.build import add_adsorbate
from ase.visualize import view

from scipy.spatial import Delaunay
import numpy as np
    
def add_adsorbate_to_poscar(input_poscar, output_poscar, adsorbate_atom, site, height):
    """
    Add an adsorbate to a surface at all possible positions of a specified site type.

    Parameters:
    - input_poscar (str): Path to the input POSCAR file representing the surface.
    - output_poscar (str): Path to the output POSCAR file to write the new structure.
    - adsorbate_atom (str): Symbol of the adsorbate atom (e.g., 'H').
    - site (str): Site type to add the adsorbate ('top', 'bridge', or 'hole').
    - height (float): Height above the site to place the adsorbate.

    Returns:
    None
    """
    # Read the input POSCAR file to get the slab structure
    slab = read(input_poscar)

    # Create the adsorbate
    adsorbate = Atoms(adsorbate_atom)

    # Identify the surface atoms (e.g., atoms with the highest z-coordinate)
    z_coords = [atom.position[2] for atom in slab]
    max_z = max(z_coords)
    surface_atoms = [i for i, atom in enumerate(slab) if atom.position[2] == max_z]

    # Define a list to store all possible positions for the specified site type
    positions = []

    # Find all possible positions for the specified site type
    if site == 'top':
        for i in surface_atoms:
            positions.append(slab[i].position[:2])
    elif site == 'bridge':
        for i in surface_atoms:
            for j in surface_atoms:
                if i < j:
                    positions.append((slab[i].position[:2] + slab[j].position[:2]) / 2)
    elif site == 'hole':
        # Get the 2D positions of the surface atoms
        surface_positions = [slab[i].position[:2] for i in surface_atoms]

        # Find the Delaunay triangulation of the surface positions
        tri = Delaunay(surface_positions)

        # Find the centroids of the triangles in the triangulation to get the hollow sites
        for simplex in tri.simplices:
            triangle = surface_positions[simplex[0]], surface_positions[simplex[1]], surface_positions[simplex[2]]
            centroid = np.mean(triangle, axis=0)
            positions.append(centroid)
    else:
        raise ValueError("Invalid site specified. Choose from 'top', 'bridge', or 'hole'.")

    # Add the adsorbate at each found position
    for position in positions:
        add_adsorbate(slab, adsorbate, height, position=position)

    # Write the new structure to the output POSCAR file
    write(output_poscar, slab)

    # Visualize the new structure
    view(slab)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add an adsorbate to a surface at all possible positions of a specified site type.')
    
    parser.add_argument('input_poscar', type=str, help='Path to the input POSCAR file representing the surface.')
    parser.add_argument('output_poscar', type=str, help='Path to the output POSCAR file to write the new structure.')
    parser.add_argument('adsorbate_atom', type=str, help='Symbol of the adsorbate atom (e.g., "H").')
    parser.add_argument('site', type=str, choices=['top', 'bridge', 'hole'], help='Site type to add the adsorbate ("top", "bridge", or "hole").')
    parser.add_argument('height', type=float, help='Height above the site to place the adsorbate.')

    args = parser.parse_args()

    add_adsorbate_to_poscar(args.input_poscar, args.output_poscar, args.adsorbate_atom, args.site, args.height)
