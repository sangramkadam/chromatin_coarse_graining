# chromatin_coarse_graining

File description:
run.in        - Script file for LAMMPS 
conf.dat      - Initial input file with polymer configuration and information about extra chromatin contacts
pos.xyz       - Sample output position data file
com_pos.c     - Analysis code to find center of mass positions of CG beads 
cg_analysis.c - Uses the center of mass file obtained from com_pos.c to measure coarse-graining properties

How to run:
1. Run LAMMPS using command './lmp_serial <run.in'  (Refer LAMMPS documentation for installation).
2. Compile com_pos.c using 'gcc com_pos.c -lm' and run using './a.out'.
3. Compile cg_analysis.c using 'gcc cg_analysis.c -lm' and run using './a.out'.
