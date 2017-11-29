# cryoFIT

- cryoFIT version 0.9 (alpha testing now)

- Confirmed to work well at macOS 10.11.6 and Ubuntu 16.

- Using Centos 7 is OK, but it didn't calculate CC once.

- How to install cryoFIT: Even if user has cryoFIT folder in <user_phenix>/modules, still the user needs to install gromacs_cryo_fit
by python <user_phenix>/modules/cryo_fit/steps/0_install_cryo_fit/install_cryo_fit.py gromacs-4.5.5_cryo_fit_added.zip

- The gromacs-4.5.5_cryo_fit_added.zip can be downloaded by visiting https://github.com/kimdn/cryo_fit_zip and clicking gromacs-4.5.5_cryo_fit_added.zip

- How to install openmpi (maybe needed if user's current mpirun doesn't work): python <user_phenix>/modules/cryo_fit/command_line/install_openmpi.py openmpi-2.1.1.tar.gz

- Forcefield (in .../gromacs-4.5.5.MDfit_serdal/share/gromacs/top) denoted as amber03.ff is actually amber03_gdp_MIA_SEP.ff

- To dos by developers:
Try gromocs forcefield to make ribosome can be run for cryoFIT
