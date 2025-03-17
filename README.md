Overview of the script:
Reads CIF files from an input folder or uses a precomputed similarity matrix.
Optionally strips terminal atoms (except S, O, N) from molecular structures to focus on core frameworks.
Computes packing similarity between all pairs of structures using CCDC’s PackingSimilarity tool.
Builds a similarity matrix based on the number of matched molecules or RMSD (root-mean-square deviation).
Visualizes the matrix as a heatmap with labeled axes.
Saves outputs: The similarity matrix as a text file, the heatmap as a PNG, and optionally detailed packing similarity results (text file and MOL2 overlays).
It leverages numpy for matrix operations, matplotlib for plotting, and ccdc for crystallographic functionality, with a command-line interface for configuration.
