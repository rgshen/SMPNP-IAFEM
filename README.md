# bzlu-Group-srg-SMPNP-IAFEM
This program, SMPNP-IAFEM, is developed based on the toolbox Parallel Hierarchical Grid (PHG), whch is currently under active development at State Key Laboratory of Scientific and Engineering Computing of Chinese Academy of Sciences. Please refer to http://lsec.cc.ac.cn/phg/ to get more details about the installation and usage of PHG. Afer PHG has been successfully installed, you need to change the content of 'PHG_MAKEFILE_INC' in Makefile to the installation path of PHG on your computer. Then you can compile and run our program on a parallel computer by running 'smpnp_sphere_run.sh' for the sphere simulation or 'smpnp_gA_run.sh' for the gA channel simulation.

You can modify 'smpnp_xx_run.sh' to adjust some parameters of our program.

-np 32 # the number of cores

-mesh_file ./nano_mesh/nano_l1.mesh # the mesh file of a nanopore

-bias_anode 5.0 # the applied voltage, the unit is V

-bias_surface_charge -1.0e-7 # the surface charge density on the nanopore, the unit is C/cm^2

-analytic_density 1.0e-3 # the bulk concentration density, the unit is M
