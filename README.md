# Predicting scale-dependent chromatin polymer properties from systematic coarse-graining

##### Table of Contents  
[Headers](#headers)  
[Emphasis](#emphasis)  
...snip...    
<a name="headers"/>
## Headers

## System requirements
- LAMMPS is an open source software. For details [click here](https://docs.lammps.org/Install.html).
- The analysis codes were tested on a Linux system (Ubuntu).
- The gcc compiler with math library was used to compile the analysis codes.

## File description:
- run.in        - Script file for LAMMPS 
- conf.dat      - Initial input file with polymer configuration and information about extra chromatin contacts
- pos.xyz       - Demo output position data file from simulation
- com_pos.c     - Analysis code to find center of mass positions of CG beads 
- cg_analysis.c - Uses the center of mass file obtained from com_pos.c to measure coarse-graining properties

## How to run:
- Refer [LAMMPS documentation](https://docs.lammps.org/Install.html) for installation. Run LAMMPS using command 
```
./lmp_serial <run.in
```
- Compile and run com_pos.c using 
```
gcc com_pos.c -lm
./a.out
```
- Compile and run cg_analysis.c using 
```
gcc cg_analysis.c -lm
./a.out
```


Expected output: The simulation generates single position trajectory file. The analysis codes compute average values and distributions of radius of gyration, bond length, and bond angle.


Typical installation for LAMMPS software is around 30 minutes. The simulation run for a single sample trajectory takes around 20 to 30 minutes. Analysis codes for single trajectory take less than a minute.
