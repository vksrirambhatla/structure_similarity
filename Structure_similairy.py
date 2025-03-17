import sys
import argparse
import os
import matplotlib

matplotlib.use('Agg')
from ccdc.io import EntryReader, CrystalWriter, MoleculeWriter
from ccdc.crystal import PackingSimilarity
import numpy as np
import matplotlib.pyplot as plt

def strip_terminal(name, reader):
    allowed_terminal_atoms = ["S", "O", "N"]
    xtal_writer = CrystalWriter(name + "_stripped.cif")

    for entry in reader:
        crystal = entry.crystal
        crystal.assign_bonds()
        molecule = crystal.molecule
        atom_diff = -1
        while atom_diff != 0:
            n_atoms_before = len(molecule.atoms)
            for atom in molecule.atoms:
                if len(atom.bonds) == 1 and atom.atomic_symbol not in allowed_terminal_atoms:
                    molecule.remove_atom(atom)
            n_atoms_after = len(molecule.atoms)
            atom_diff = n_atoms_after - n_atoms_before

        crystal.molecule = molecule
        xtal_writer.write(crystal)

def read_cif_files_from_folder(input_folder):
    cif_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.cif')]
    return cif_files

def main(input_folder, output_folder, matrix_file, n_ps_mols, output_ps_results, conf_threshold, ps_angles, ps_distances, strip,
         n_struct, allow_mol_diff, cluster_mode, pad_length):

    cif_files = read_cif_files_from_folder(input_folder)
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    ps = PackingSimilarity()
    ps.settings.ignore_hydrogen_positions = True
    ps.settings.ignore_bond_types = True
    ps.settings.match_entire_packing_shell = False

    if allow_mol_diff:
        ps.settings.allow_molecular_differences = True
        ps.settings.ignore_hydrogen_counts = True
        ps.settings.ignore_bond_counts = True
    else:
        ps.settings.allow_molecular_differences = False
        ps.settings.ignore_hydrogen_counts = False
        ps.settings.ignore_bond_counts = False

    refcodes = []
    input_name = os.path.basename(input_folder)

    print("--------------------------------------------------------")

    if not matrix_file:
        print("Reading CIF files from:", input_folder)
        
        if strip:
            stripped_cif_files = []
            for cif_file in cif_files:
                structure_reader = EntryReader(cif_file)
                cif_name = os.path.splitext(os.path.basename(cif_file))[0]
                strip_terminal(cif_name, structure_reader)
                stripped_cif_files.append(cif_name + "_stripped.cif")
                structure_reader.close()
            cif_files = stripped_cif_files

        structure_size = len(cif_files)

        matrix = np.zeros((structure_size, structure_size))

        print("Generating matrix of packing similarities")

        overlay_folder = os.path.join(output_folder, "overlays")
        if output_ps_results:
            g = open(os.path.join(output_folder, "packing_similarity_results.txt"), "w")
            g.write("Packing Similarity Analysis for: " + input_folder + "\n")
            if not os.path.exists(overlay_folder):
                os.makedirs(overlay_folder)

        for i, cif_file_i in enumerate(cif_files):
            entry_i = EntryReader(cif_file_i)[0]
            crystal_i = entry_i.crystal
            refcodes.append(str(i+1))

            for j, cif_file_j in enumerate(cif_files):
                if i == j:
                    matrix[i, i] = ps.settings.packing_shell_size
                    continue
                entry_j = EntryReader(cif_file_j)[0]
                crystal_j = entry_j.crystal

                result = ps.compare(crystal_i, crystal_j)

                if result is not None:
                    if result.nmatched_molecules != 1:
                        matrix[i, j] = result.nmatched_molecules
                        matrix[j, i] = result.nmatched_molecules
                    elif result.rmsd < conf_threshold:
                        matrix[i, j] = result.nmatched_molecules
                        matrix[j, i] = result.nmatched_molecules
                    else:
                        matrix[i, j] = 0
                        matrix[j, i] = 0
                    if output_ps_results:
                        g.write(f"{entry_i.identifier} {entry_j.identifier}: {int(matrix[i, j])} molecules\n")
                        overlay_writer = MoleculeWriter(os.path.join(overlay_folder, 
                                                                     f"overlay_{entry_i.identifier}_{entry_j.identifier}.mol2"))
                        mols = result.overlay_molecules()
                        for mol in mols:
                            overlay_writer.write(mol)
                else:
                    matrix[i, j] = 0
                    matrix[j, i] = 0
                    if output_ps_results:
                        g.write(f"{entry_i.identifier} {entry_j.identifier}: {int(matrix[i, j])} molecules\n")

        np.savetxt(os.path.join(output_folder, input_name + '_similarity_matrix.txt'), matrix, delimiter=',')
        print(f"Packing similarity matrix saved to {os.path.join(output_folder, 'similarity_matrix.txt')}")
        if output_ps_results:
            g.close()
    else:
        print(f"Reading input matrix: {matrix_file}")
        matrix = np.loadtxt(matrix_file, delimiter=',')
        structure_size = len(matrix)

    x = np.arange(0, structure_size + 1, 1)
    y = np.arange(0, structure_size + 1, 1)

    # Set figure size larger for more space
    plt.figure(figsize=(10, 10))  # Adjust size as needed

    # Create heat map
    plot = plt.pcolor(x, y, matrix, cmap=plt.get_cmap('rainbow', (n_ps_mols - 1)), vmin=1, vmax=n_ps_mols)
    plt.colorbar(plot, ticks=range(1, n_ps_mols+1))

    # Rotate x-axis labels for better readability
    plt.xticks(np.arange(0.5, len(cif_files)+0.5), [os.path.basename(f) for f in cif_files], rotation=90, fontsize=8, ha='right')
    plt.yticks(np.arange(0.5, len(cif_files)+0.5), [os.path.basename(f) for f in cif_files], fontsize=8)

    # Add axis labels
    plt.xlabel("CIF Files")
    plt.ylabel("CIF Files")

    # Save the heatmap
    plt.savefig(os.path.join(output_folder, input_name + "_heat_map.png"), dpi=300, bbox_inches='tight')
    print(f"Packing similarity heat map saved to {os.path.join(output_folder, input_name + '_heat_map.png')}")
    plt.close()

    print("--------------------------------------------------------")

    sys.exit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('--input_folder', help='Folder containing CIF structures.', default='C:\\vijay\\Solvates')
    parser.add_argument('--output_folder', help='Folder to save the output.', default='C:\\vijay\\Solvates')
    parser.add_argument('-m', '--matrix', type=str, help='NumPy matrix containing existing packing similarity results.', metavar='similarity_matrix.txt')
    parser.add_argument('-ns', '--n_structures', type=int, help='Number of structures to take from input set.', metavar='25')
    parser.add_argument('-nm', '--n_molecules', type=int, default=15, help='Size of molecular packing shell to use for analysis.', metavar="20")
    parser.add_argument('-o', action="store_true", help='Flag for whether to save packing similarity results (text file and mol2 overlays).')
    parser.add_argument('--allow_molecular_differences', action="store_true", help='Flag for whether to allow for molecular differences between structures (e.g. for salts).')
    parser.add_argument('--clustering_type', choices=['complete', 'single', 'average'], default='single', help='Type of clustering to employ')
    parser.add_argument('-s', '--strip', action="store_true", help='Strip all terminal atoms and alkyl chains, up to any hetero atom (O, N, S) or cyclic atom.')
    parser.add_argument('-ct', '--conf_tol', type=float, default=0.5, help='RMSD threshold for considering two conformations to be the same.')
    parser.add_argument('-at', '--angle_tol', type=float, default=25, metavar="25", help="Tolerance for angles (in degrees) used by packing similarity.")
    parser.add_argument('-dt', '--dist_tol', type=float, default=0.3, metavar="0.3", help="Tolerance for distances (in Å) used by packing similarity.")

    args = parser.parse_args()

    main(
        input_folder=args.input_folder,
        output_folder=args.output_folder,
        matrix_file=args.matrix,
        n_ps_mols=args.n_molecules,
        output_ps_results=args.o,
        conf_threshold=args.conf_tol,
        ps_angles=args.angle_tol,
        ps_distances=args.dist_tol,
        strip=args.strip,
        n_struct=args.n_structures,
        allow_mol_diff=args.allow_molecular_differences,
        cluster_mode=args.clustering_type,
        pad_length=3
    )
