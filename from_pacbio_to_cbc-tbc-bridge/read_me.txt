1. python scripts
the script to step-wise construct the bridge connecting cell-barcode (CBC) and viral clonal tracking barcode (TBC). To use the code, follow these steps: 1) download the "scripts" folder and the "00-external" folder into the same directory 2) prepare the additional input data (see below) 3) open the ".py" code in a Python interpreter (Spyder recommended) in the numerical order 4) follow the instructions inside the script

2. 00-external folder
folder containing supporting materials for data processing, including the TBC library IDs, viral backbone sequences, multiplex information, etc

3. required additional input:
matrix data from clonal tracking; fastq files from molecular bridge sequencing by PacBio; H5 file from single cell RNAseq processed by cellranger.