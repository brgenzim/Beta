Description: generation of training and benchmark sets for AlphaFold2 applications
Requirements: 
	- Voronota
	- PSI-BLAST
	- Foldseek 
	- gemmi
Further requirements
	- PDB fasta file
	- SwissProt fasta file	
	- SIFTS pdb chain and uniprot mapping
	- Local version of BIOGRID (in mitab format)
	- cutoff_dates.txt with cutoff dates

Note: Scripts were optimized for HPCs with slurm.



#run 1_pdb_summary_fixer.py 
1. Create a list of PDBs with release date and experimental method. Set paths in the script.
python3 script/1_pdb_summary.py

2. Create foldseek library
foldseek databases PDB pdb tmp

3. Create blast db for PDB (it is recommended to use only protein sequences)
makeblastdb -in pdb_protein.fas -dbtype prot

5. define cutoff dates (supplied)

6. For any PDB code, search for homologous structures and sequences.
python3 script/2_batch_struct.py collects data and use slurm to distribute job (Set paths in the script.)
or run 3_structure_search.py for all pdb codes
python3 script/3_structure_search.py [pdb_code] [pdb_code_fasta] [pdb_blast_db] [psi-blast_binary] [pdb_cif] [foldseek_pdb_lib] [foldseek_binary] [dates] [pdb_summary] [workdir] [output_file]
[pdb_code]: 4 character code
[pdb_code_fasta]: fasta file containing protein sequences for pdb_code
[pdb_blast_db]: generated on step 3
[psi-blast_binary]: path for psi-blast
[pdb_cif]: directory of unzipped version of pdb cif file
[foldseek_pdb_lib]: generated in step 2
[foldseek_binary]: path for foldseek
[dates]: generated in step 5
[pdb_summary]: generated in step 1
[workdir]: working directory
[output file]: output file

7. Collect results. Set paths in the script.
python3 script/4_collect_struct.py

8. Within new PDB chains search for interactions. Set paths in the script.
python3 script/5_batch_voronota.py
or run 6_voronota.py
python3 6_voronota.py [pdb_code] [infile] [outfile]
[pdb_code]: 4 character code
[infile]: cif file for assemblies from PDB (containing most probable oligomerization state)
[outfile]: outfile

9. Collect voronota results. Set paths in the script.
python3 script/7_collect_voronota.py

10. For all SwissProt files, search for homologous sequences. Set paths in the script.
python3 script/8_batch_seq.py collects data and use slurm to distribute job
or run 9_sequence_search.py for all uniprot sequences
python3 script/9_sequence_search.py [uniprot ac] [uniprot fasta] [pdb_blast_db] [psi_blast_binary] [sifts] [dates] [pdb_result] [pdb_summary] [workdir] [output]
[uniprot_ac]: UniProt AC
[uniprot_fasta]: UniProt sequence (single)
[pdb_blast_db]: generated on step 3
[psi_blast_binary]: path for psi-blast
[sifts]: SIFTS file
[dates]: generated in step 5
[pdb_result]: generated in step 7
[pdb_summary]: generated in step 1
[workdir]: working directory
[output]: output file

11. Collect results. Set paths in the script.
python3 script/10_collect_sequences.py

12. Get interactions from BIOGRID. Set paths in the script. 
python3 script/11_sequence_interactions.py


The result should be 4 files from steps 7, 9, 11 and 12. Output files are defined in the first rows (see path for OUT).
