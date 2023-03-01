# KSHELL - Thick-restart block Lanczos method for large-scale shell-model calculations

Noritaka Shimizu, Takahiro Mizusaki, Yutaka Utsuno, Yusuke Tsunoda

Center for Nuclear Study, The University of Tokyo, 7-3-1 Hongo, Bunkyo-ku, Tokyo 113-0033

Japan Institute of Natural Sciences, Senshu University, 3-8-1 Kanda-Jinbocho, Chiyoda-ku, Tokyo 101-8425

Japan Advanced Science Research Center, Japan Atomic Energy Agency, Tokai, Ibaraki 319-1195, Japan

https://doi.org/10.1016/j.cpc.2019.06.011

Code downloaded from https://sites.google.com/alumni.tsukuba.ac.jp/kshell-nuclear/

<details>
<summary>Abstract</summary>
<p>

  We propose a thick-restart block Lanczos method, which is an extension of the thick-restart Lanczos method with the block algorithm, as an eigensolver of the large-scale shell-model calculations. This method has two advantages over the conventional Lanczos method: the precise computations of the near-degenerate eigenvalues, and the efficient computations for obtaining a large number of eigenvalues. These features are quite advantageous to compute highly excited states where the eigenvalue density is rather high. A shell-model code, named KSHELL, equipped with this method was developed for massively parallel computations, and it enables us to reveal nuclear statistical properties which are intensively investigated by recent experimental facilities. We describe the algorithm and performance of the KSHELL code and demonstrate that the present method outperforms the conventional Lanczos method.

  Program summary
  Program Title: KSHELL

  Licensing provisions: GPLv3

  Programming language: Fortran 90

  Nature of problem: The nuclear shell-model calculation is one of the configuration interaction methods in nuclear physics to study nuclear structure. The model space is spanned by the M-scheme basis states. We obtain nuclear wave functions by solving an eigenvalue problem of the shell-model Hamiltonian matrix, which is a sparse, symmetric matrix.

  Solution method: The KSHELL code enables us to solve the eigenvalue problem of the shell-model Hamiltonian matrix utilizing the thick-restart Lanczos or thick-restart block Lanczos methods. Since the number of the matrix elements are too huge to be stored, the elements are generated on the fly at every matrix–vector product. The overhead of the on-the-fly algorithm are reduced by the block Lanczos method.

  Additional comments including restrictions and unusual features: The KSHELL code is equipped with a user-friendly dialog interface to generate a shell script to run a job. The program runs both on a single node and a massively parallel computer. It provides us with energy levels, spin, isospin, magnetic and quadrupole moments, E2/M1 transition probabilities and one-particle spectroscopic factors. Up to tens of billions M-scheme dimension is capable, if enough memory is available.

</p>
</details>


## Prerequisites

<details>
<summary>Click here for prerequisites</summary>
<p>

  * ```Python 3.10``` or newer (kshell_ui.py uses syntax specific to 3.10 and above)
    * `numpy`
    * `matplotlib` (not required but recommended)
    * `kshell-utilities` (not required but recommended)
  * ```gfortran 10.2.0``` or newer (Tested with this version, might work with older versions)
  * ```ifort 19.1.3.304``` (Alternative to gfortran. Tested with this version, might work with other versions.)
  * ```openblas```
  * ```lapack```

  Use `gfortran` Fortran compiler if you plan on running KSHELL on your personal computer and use `ifort` for the Fram supercomputer.
</p>
</details>


## KSHELL on Fram and Betzy

  <details>
  <summary>Click here for KSHELL on Fram and Betzy</summary>
  <p>

  ### Compilation on Fram and Betzy with MPI
  Start by loading the necessary modules which contain the correct additional software to run `KSHELL`. The `intel/2020b` module contains the correct `ifort` version as well as `blas` and `lapack` (double check this), and the module `Python/3.8.6-GCCcore-10.2.0` gives us the correct `Python` version. Load the modules in this order:
  ```
  module load intel/2020b
  module load Python/3.8.6-GCCcore-10.2.0
  ```
  Now, clone this repository to the desired install location. Navigate to the `<install_location>/src/` directory and edit the `Makefile`. We will use the MPI ifort wrapper `mpiifort` to compile `KSHELL`, so make sure that `FC = mpiifort` is un-commented and that all other `FC = ` lines are commented. Comment with `#`. Remember to save the file. Still in the `<install_location>/src/` directory, run the command `make`, and `KSHELL` will be compiled.

  <details>
  <summary>Click here to see the terminal output from the compilation process</summary>
  <p>

    ```
    $ make
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c constant.f90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c model_space.f90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c lib_matrix.F90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c class_stopwatch.F90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c partition.F90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c wavefunction.F90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c rotation_group.f90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c harmonic_oscillator.f90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c operator_jscheme.f90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c operator_mscheme.f90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c bridge_partitions.F90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c sp_matrix_element.f90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c interaction.f90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c bp_io.F90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c lanczos.f90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c bp_expc_val.F90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c bp_block.F90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c block_lanczos.F90
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c kshell.F90
    mpiifort -O3 -qopenmp -no-ipo -DMPI -o kshell.exe kshell.o model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -mkl
    mpiifort -O3 -qopenmp -no-ipo -DMPI  -c transit.F90
    mpiifort -O3 -qopenmp -no-ipo -DMPI -o transit.exe transit.o model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -mkl
    mpiifort -O3 -qopenmp -no-ipo -DMPI -o count_dim.exe count_dim.f90 model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -mkl
    cp kshell.exe transit.exe count_dim.exe ../bin/
    ```

  </p>
  </details>

  `KSHELL` is now compiled! To remove the compiled files and revert back to the starting point, run `make clean` in the `src/` directory.

  ### Queueing job script on Fram and Betzy
  Create a directory in which to store the output from `KSHELL`. In this directory, run `python <install_location>/bin/kshell_ui.py` and follow the instructions on screen. The shell script grenerated by `kshell_ui.py` must begin with certain commands wich will be read by the job queue system, `slurm`. The needed commands will automatically be added to the executable shell script if the keyword `fram` or `betzy` is entered in the first prompt of `kshell_ui.py`. See a section further down in this document for general instructions on how to use `kshell_ui.py`. When the executable shell script has been created, put it in the queue by

  ```
  sbatch executable.sh
  ```

  To see the entire queue, or to filter the queue by username, use

  ```
  squeue
  squeue -u <username>
  ```

  The terminal output from the compute nodes is written to a file, `slurm-*.out`, which is placed in the `KSHELL` output directory you created. Use

  ```
  tail -f slurm-*.out
  ```

  to get a live update on the last 10 lines of terminal output from the compute nodes. If you put in your e-mail address in the executable shell script, you will get an e-mail when the program starts and when it ends (per 2021-12-10, the mailing system is not operative). Following is an example of the commands which must be in the first line of the executable shell script which is generated by `kshell_ui.py`. For running 10 nodes with 32 cores each with an estimated calculation time of 10 minutes on Fram:

  <details>
  <summary>Click here to see the Fram commands</summary>
  <p>

    ```
    #!/bin/bash
    #SBATCH --job-name=Ar28_usda
    #SBATCH --account=<enter account name here (example NN9464K)>
    ## Syntax is d-hh:mm:ss
    #SBATCH --time=0-00:10:00
    #SBATCH --nodes=10
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task=32
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=<your e-mail here>
    module --quiet purge
    module load intel/2020b
    module load Python/3.8.6-GCCcore-10.2.0
    set -o errexit
    set -o nounset
    ```

  </p>
  </details>
    
  For running a job on Betzy with 64 nodes with an estimated time of 1 day, using all 256 (virtual (SMT)) cores per node effectively, the slurm commands look like:
    
  <details>
  <summary>Click here to see the Betzy commands</summary>
  <p>
    ```
    #!/bin/bash
    #SBATCH --job-name=V50_gxpf1a
    #SBATCH --account=<enter account name here (example NN9464K)>
    ## Syntax is d-hh:mm:ss
    #SBATCH --time=0-01:00:00
    #SBATCH --nodes=64
    #SBATCH --ntasks-per-node=8
    #SBATCH --cpus-per-task=16
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=<your e-mail>
    module --quiet purge
    module load intel/2020b
    module load Python/3.8.6-GCCcore-10.2.0
    set -o errexit
    set -o nounset
    export OMP_NUM_THREADS=32
    ```
  </p>
  </details>
  
  The command `export OMP_NUM_THREADS=32` forces 256 virtual cores to be used instead of 128 physical cores per node. SMT is beneficial to use with KSHELL, so use this option for better performance! `--ntasks-per-node=8` specifies 8 MPI ranks per node, and `--cpus-per-task=16` specifies 16 OMP threads per MPI rank (and is extended to 32 by `export OMP_NUM_THREADS=32` which in total per node utilizes 8*32 = 256 virtual cores). The Betzy documentation states that this mix of MPI ranks and OMP threads yields better performance than a pure MPI or pure OMP setup.
    
  Note that the modules must be explicitly loaded in the script file since the modules you load to the login node does not get loaded on the compute nodes. The login node is the computer you control when you SSH to `<username>@fram.sigma2.no` and the compute nodes are other computers which you control via the `slurm` queue system. If you need any other modules loaded, you must add these to the executable shell script. Now, just wait for the program to run its course!

  </p>
  </details>

## KSHELL on your PC
  
  <!-- ### Installation on Ubuntu -->
    
  <details>
  <summary>Installation on Ubuntu</summary>
  <p>

  KSHELL probably works fine on any Linux distro as long as you install the correct versions of Fortran and Python. Following is a recipe for installing and compiling on Ubuntu 20.04.2 LTS.

  #### Fortran compiler
  We start by installing a compatible version of `gfortran`. To get a version newer than 9, we must first add the Ubuntu Toolchain repository:
  ```
  sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
  ```
  Then, install `gfortran` version 10 with:
  ```
  sudo apt install gfortran-10
  ```
  And check that the newly installed Fortran compiler is of version 10.2.0 or above:
  ```
  gfortran-10 --version
  ```
  If the version is incorrect, try installing `gfortran` version 11 instead. Now, install the correct `blas` and `lapack` versions with
  ```
  sudo apt install libopenblas-dev
  ```
  If this specific version does not work or exist for your system, run
  ```
  apt search openblas
  ```
  and try a few different versions to see which one works. You check that the version is correct by compiling `KSHELL` and seeing whether the compile completes or not. Compilation instructions are following.

  #### Python
  For installing the correct version of Python, it is highly recommended to install an environment management system like `miniconda` as to not mess up any other Python dependencies your system has, and to easily download the exact version needed. Start by downloading the latest release of `miniconda` ([alternative downloads here](https://docs.conda.io/en/latest/miniconda.html)):
  ```
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  ```
  Run the installer:
  ```
  bash Miniconda3-latest-Linux-x86_64.sh
  ```
  Accept the terms of service. Choose all default settings except when the installer asks if it should initialize by running conda init. Choose yes. If you have trouble with initializing conda, for example
  ```
  > conda
  conda: command not found
  ```
  cd to `<install_location>/anaconda3/bin` and initialize conda from there (replace `<install_location>` with the path to where you downloaded the file `Miniconda3-latest-Linux-x86_64.sh`). If you for example use `fish` instead of `bash` (you should!), then initialize with
  ```
  ./conda init fish
  ```
  At this point, please close your terminal and open a new one. When the initialization is complete, create an environment named `kshell` with `Python 3.10` with:
  ```
  conda create --name kshell python=3.10
  ```
  Activate the environment with:
  ```
  conda activate kshell
  ```
  Note that any additional Python packages may be installed normally with `pip`. You do not have to type `pip3` or `python3` because conda maps `pip` and `python` to the requested version of Python. The `kshell` environment is only active within your terminal session and does not interfere with any other Python dependencies on your system. This is one of the main reasons why I recommend using an environment manager. Should you at any time in the future need a different version of Python (newer or older), simply create a conda environment with the appropriate Python version. You can see the currently active conda environment in the bottom right or left corner of your terminal. The default environment is called `(base)` and if you have followed these instructions you will see that the active environment is `(kshell)`.

  As an alternative to using the forementioned conda approach, you can download `Python 3.10` with your distro's appropriate packet manager. Use `apt search python` to find the correct name of Python version 3.10 or newer.

  #### Compile KSHELL
  We are now ready to actually install `KSHELL`. Navigate to the directory where you want to install `KSHELL` and clone this repository:
  ```
  git clone https://github.com/GaffaSnobb/kshell.git
  ```
  When the clone is successful use `ls` to see that you have a new directory called `kshell`. Navigate to the `kshell/src/` directory and edit the `Makefile` with your favorite editor. Change `FC = gfortran` to `FC = gfortran-10` (or `-11` if you installed version 11) and make sure that line is un-commented. All other `FC` declarations should be commented or deleted. Save the changes. Still in the `src/` directory, run
  ```
  make
  ```
  to compile. The output should be something like this (mismatch warnings are normal):
  
  <details>
  <summary>Click to see normal terminal output</summary>
  <p>

  ```
  > make
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c constant.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c model_space.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c lib_matrix.F90
  lib_matrix.F90:304:29:

    304 |     call dlarnv(1, iseed, 1, r )
        |                             1
  ......
    312 |     call dlarnv(1, iseed, n, r)
        |                             2
  Warning: Rank mismatch between actual argument at (1) and actual argument at (2) (rank-1 and scalar)
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c class_stopwatch.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c partition.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c wavefunction.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c rotation_group.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c harmonic_oscillator.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c operator_jscheme.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c operator_mscheme.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c bridge_partitions.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c sp_matrix_element.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c interaction.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c bp_io.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c lanczos.f90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c bp_expc_val.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c bp_block.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c block_lanczos.F90
  block_lanczos.F90:548:12:

    548 |             vr(i*nb+1,1), size(vr,1), &
        |            1
  ......
    577 |             -1.d0, vin(i*nb+1, 1), size(vin,1), an, size(an,1), &
        |                                                2
  Warning: Element of assumed-shape or pointer array as actual argument at (1) cannot correspond to actual argument at (2)
  block_lanczos.F90:250:20:

    250 |               1.d0, vi, nc, &
        |                    1
  ......
    577 |             -1.d0, vin(i*nb+1, 1), size(vin,1), an, size(an,1), &
        |                   2
  Warning: Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-2)
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c kshell.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch -o kshell.exe kshell.o model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -llapack -lblas -lm
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch  -c transit.F90
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch -o transit.exe transit.o model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -llapack -lblas -lm
  gfortran-10 -O3 -fopenmp -fallow-argument-mismatch -o count_dim.exe count_dim.f90 model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -llapack -lblas -lm
  cp kshell.exe transit.exe count_dim.exe ../bin/
  ```

  </p>
  </details>
  
  If the output of your terminal is like the expected terminal output listed above then `KSHELL` is compiled correctly and ready to use. See a section further down in this readme for instructions on how to run `KSHELL`. If your terminal reports that the command `gfortran` cannot be found then you need to correctly edit your `Makefile` with `FC = <correct gfortran command>`. It might be `gfortran`, `gfortran-10`, `gfortran-11` or something similar.

  </p>
  </details>

  <!-- ### Installation on macOS -->
    
  <details>
  <summary>Installation on macOS</summary>
  <p>

  #### Homebrew
  `Homebrew` is a packet manager for macOS similar to `apt` for Ubuntu and frankly, every scientist using macOS should have `Homebrew` installed. Install with ([or see detailed install instructions here](https://brew.sh)):
  ```
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  ```

  #### Fortran
  Install the newest Fortran compiler with (per 2021-09-29 version 11.2.0 will be installed):
  ```
  brew install gfortran
  ```
  and check that the version is equal to or greater than 10.2.0 by (version x.y.z with x > 10 should also be fine):
  ```
  gfortran --version
  ```

  #### Python
  For installing the correct version of Python, it is highly recommended to install an environment management system like `miniconda` as to not mess up any other Python dependencies your system has, and to easily download the exact version needed. Start by `cd`ing to your downloads directory with
  ```
  cd ~/Downloads
  ```
  and download the latest release of `miniconda` ([alternative downloads here](https://docs.conda.io/en/latest/miniconda.html)) by:
  ```
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  ```
  Then, run the installer which you just downloaded with:
  ```
  bash Miniconda3-latest-MacOSX-x86_64.sh
  ```
  Follow the on-screen instructions and accept the terms of service. Choose all default settings except when the installer asks if it should initialize by running conda init. For that option select yes. After the installation is complete, close and re-open your terminal session for the conda installation to function properly. Check that conda has been installed by running the command `conda`. If a helpful description of conda shows up in your terminal then the installation is complete. However, if you see
  ```
  > conda
  conda: command not found
  ```
  then conda has not been installed properly. To fix this, `cd` to `<install_location>/anaconda3/bin` and initialize conda from there. Replace `<install_location>` with the path to where you downloaded the file `Miniconda3-latest-Linux-x86_64.sh` which should be `~/Downloads` if you followed these instructions. Then, initialise conda with:
  ```
  ./conda init bash
  ```
  Replace `bash` if you are using a different shell, for example `zsh` or `fish` to mention a few. Close your terminal and open a new session, then run the command `conda` once more to see that the installation is proper. When the initialization is complete, create an environment named `kshell` with `Python 3.10`:
  ```
  conda create --name kshell python=3.10
  ```
  Activate the environment with:
  ```
  conda activate kshell
  ```
  Note that any additional Python packages may be installed normally with `pip`. You do not have to type `pip3` or `python3` because conda maps `pip` and `python` to the requested version of Python. For example, to install `numpy`, simply run:
  ```
  pip install numpy
  ```
  The `kshell` environment is only active within your terminal session and does not interfere with any other Python dependencies on your system. You can see the currently active conda environment in the bottom right or left corner of your terminal. The default environment is called `(base)` and if you have followed these instructions you will see that the active environment is `(kshell)`. This is one of the main reasons why I recommend using an environment manager. You may mess around however much you like inside of any conda environment **other than** `base` without having to worry about a thing. Should you at any time in the future need a different version of Python (newer or older), simply create a conda environment with the appropriate Python version.
  
  You: I need Python 2.7 to run some old code. Lets open the good ole' terminal.
  Apple: Tough luck. 2.7 is not included in newer versions of macOS.
  You: But wait... I have conda!
  ```
  conda create --name <name of environment here> python=2.7
  conda activate <name of environment here>
  python --version
  > Python 2.7.18
  ```
  You: Noice!

  As an alternative to using the forementioned conda approach, you can download `Python 3.10` with `brew`. Use `brew search python` to find the correct name of Python version 3.10 or newer and then use `brew install <correct name>` to install it. Note that you may have to use `python3.10 myfile.py` to actually use version 3.10.

  #### Compile KSHELL
  We are now ready to actually compile (install) `KSHELL`. Navigate to the directory where you want to install `KSHELL`. This guide assumes that you use your home directory (`cd ~/`) but you may use another path if you'd like. Note that you do not need to create a directory called kshell (or whatever), because such a directory will be created during the installation. In your home directory, clone this repo with the command:
  ```
  git clone https://github.com/GaffaSnobb/kshell.git
  ```
  which copies all the needed files to your home inside a directory called `kshell`. `cd` to the `~/kshell/src/` directory and run the command
  ```
  make
  ```
  to compile `KSHELL`. The output should be something like this (mismatch warnings are normal):
  
  <details>
  <summary>Click to see normal terminal output</summary>
  <p>

  ```
  > make
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c constant.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c model_space.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c lib_matrix.F90
  lib_matrix.F90:304:29:

    304 |     call dlarnv(1, iseed, 1, r )
        |                             1
  ......
    312 |     call dlarnv(1, iseed, n, r)
        |                             2
  Warning: Rank mismatch between actual argument at (1) and actual argument at (2) (rank-1 and scalar)
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c class_stopwatch.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c partition.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c wavefunction.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c rotation_group.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c harmonic_oscillator.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c operator_jscheme.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c operator_mscheme.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c bridge_partitions.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c sp_matrix_element.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c interaction.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c bp_io.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c lanczos.f90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c bp_expc_val.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c bp_block.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c block_lanczos.F90
  block_lanczos.F90:548:12:

    548 |             vr(i*nb+1,1), size(vr,1), &
        |            1
  ......
    577 |             -1.d0, vin(i*nb+1, 1), size(vin,1), an, size(an,1), &
        |                                                2
  Warning: Element of assumed-shape or pointer array as actual argument at (1) cannot correspond to actual argument at (2)
  block_lanczos.F90:250:20:

    250 |               1.d0, vi, nc, &
        |                    1
  ......
    577 |             -1.d0, vin(i*nb+1, 1), size(vin,1), an, size(an,1), &
        |                   2
  Warning: Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-2)
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c kshell.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch -o kshell.exe kshell.o model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -llapack -lblas -lm
  gfortran -O3 -fopenmp -fallow-argument-mismatch  -c transit.F90
  gfortran -O3 -fopenmp -fallow-argument-mismatch -o transit.exe transit.o model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -llapack -lblas -lm
  gfortran -O3 -fopenmp -fallow-argument-mismatch -o count_dim.exe count_dim.f90 model_space.o interaction.o harmonic_oscillator.o constant.o rotation_group.o sp_matrix_element.o operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o partition.o  wavefunction.o  bridge_partitions.o bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o -llapack -lblas -lm
  cp kshell.exe transit.exe count_dim.exe ../bin/
  ```

  </p>
  </details>

  If the output of your terminal is like the expected terminal output listed above then `KSHELL` is compiled correctly and is ready to use. See a section further down in this readme for instructions on how to run `KSHELL`. If your terminal reports that the command `gfortran` cannot be found then you need to correctly edit your `~/kshell/src/Makefile` with `FC = <correct gfortran command>`. It might be `gfortran`, `gfortran-10`, `gfortran-11` or something similar. To check, type the command directly into your terminal, hit enter, and see if your terminal can find the command:
  ```
  gfortran
  > gfortran: command not found
  ```
  means that `gfortran` is not a valid command. However:
  ```
  gfortran-12
  > gfortran-12: fatal error: no input files
  > compilation terminated.
  ```
  allthough fatal, means that the command `gfortran-12` is a valid command. Edit the `~/kshell/src/Makefile` with `FC = gfortran-12` in this case.
  
  
  Another possible problem is that your terminal reports this:
  ```
  ld: library not found for -llapack
  collect2: error: ld returned 1 exit status
  make: *** [kshell.exe] Error 1
  ```
  This means that `lapack` is either not installed or not in your `LIBRARY_PATH`. `brew install lapack` should solve this problem. Same for `openblas`. When you now try to compile `KSHELL` again, make sure to clean first:
  ```
  make clean
  make
  ```
  </p>
  </details>

## Usage

  <!-- #### General usage -->

  <details>
  <summary>General usage</summary>
  <p>

  We will here use 20Ne as an example. Create a directory where you want to place the **results** from `KSHELL`. Note that this directory should be at a separate location from where you installed `KSHELL`. For example, create a directory in your home and enter the newly created directory:
  ```
  cd ~/
  mkdir -p kshell_results/ne20
  cd kshell_results/ne20
  ```
  I use the name `ne20` and not `20ne` because some applications, like Python, do not support variable, function, etc. naming with numbers at the start, and I'm trying to be consise and consequent. Note that you now have an install directory: `~/kshell` and a result directory: `~/kshell_results` in your home. Be sure that you understand the difference between these two directories. The former is the location of the actual `KSHELL` program files, and the latter is the location where you wish to place the results from your `KSHELL` calculations. Don't mix these up. You woldn't place your Word documents inside the installation folder of Microsoft Office, would you? Now, inside the **results directory** for 20Ne `~/kshell_results/ne20`, initialise the `KSHELL` setup process by running the command:
  ```
  python ~/kshell/bin/kshell_ui.py
  ```
  This will start the preparations for your 20Ne calculation where you will be asked a series of questions. First, you'll be asked whether you want to use `MPI` (Message Parsing Interface) or not. `MPI` is used for parallelization over multiple nodes (computers) and is mainly applicable for running `KSHELL` on supercomputers. Parallelization over several cores per node is administered by `OpenMP` and is active even though you do not choose `MPI` here. For a regular PC, choose `n`:
  ```
  MPI parallel? Y/N/preset, n nodes (default: N,  TAB to complete) : n
  ```
  You are now asked to choose the model space you wish to use. 20Ne has 10 protons and 10 neutrons which makes the doubly magic 8p 8n core suitable for the inert core. 0d5/2, 1s1/2 and 0d3/2 will then be the model space where the valence nucleons can move about. This is the `USD` model space. Take a look at [this figure](https://periodic-table.org/wp-content/uploads/2019/05/Shell-model-of-nucleus.png) and see if you agree (note the different notation conventions, nlj and (n+1)lj (N = 2n + l)). We choose `usda.snt` for this example.
  ```
  model space and interaction file name (.snt)
  (e.g. w or w.snt,  TAB to complete) : usda.snt
  ```
  Now we specify the nuclide. Here you may enter either the number of valence protons and neutrons or the isotope abbreviation (20ne or ne20, upper or lower case does not matter). 20Ne has 2 valence protons and 2 valence neutrons outside the 8p 8n core, so the input may either be `2, 2` or `ne20`:
  ```
  number of valence protons and neutrons
  (ex.  2, 3 <CR> or 9Be <CR>)    <CR> to quit : ne20
  ```
  We are now prompted for the name of the executable shell script. Press the enter key for the default name:
  ```
  name for script file (default: Ne20_usda ):
  ```
  Choose which angular momentum levels you want to calculate and how many. The default value is to calculate the 10 lowest lying levels. See a section later in this document on how to choose specific angular momenta and parities. For this example we choose the default value (just press enter):
  ```
  J, parity, number of lowest levels
    (ex. 10          for 10 +parity, 10 -parity levels w/o J-proj. (default)
        -5           for lowest five -parity levels,
        0+3, 2+1     for lowest three 0+ levels and one 2+ levels,
        1.5-1, 3.5+3 for lowest one 3/2- levels and three 7/2+ levels
        range        for a range of levels) :
  ```
  We are now asked for truncation information. The model space is small and the number of nucleos is low, so we dont need to truncate this system. The default is no truncation. 20Ne in the `USD` model space only allows positive parity levels, so we are only asked for truncation of the positive parity levels. See a section later in this document for truncation details. Choose the default value of no truncation:
  ```
  truncation for "+" parity state in  Ne20_usda_p.ptn
  truncation scheme ?
        0 : No truncation (default)
        1 : particle-hole truncation for orbit(s)
        2 : hw truncation
        3 : Both (1) and (2)

  ```
  At this point we are asked whether we want to edit any other parameters, like the proton and neutron effective charges, the gyroscopic spin factor and the number of Lanczos iterations. Change these to your needs (tab complete is supported). In this demo, we'll leave them to the default values:
  ```
  Modify parameters?
  Example: maxiter = 300 for parameter change or <CR> for no more modification.
  Available paramters are:
  ['max_lanc_vec', 'maxiter', 'n_restart_vec', 'hw_type', 'mode_lv_hdd', 'n_block', 'eff_charge', 'gl', 'gs', 'beta_cm', 'fn_int', 'is_obtd', 'is_ry_sum', 'is_calc_tbme', 'sq', 'quench', 'is_tbtd']


 --- set parameters ---
  beta_cm = 0.0
  eff_charge = 1.5, 0.5,
  gl = 1.0, 0.0,
  gs = 5.585, -3.826,
  hw_type = 2
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  n_block = 0
  n_restart_vec = 10

  :
  ```
  The transition probabilities are calculated by default, but they can be omitted. For this example we want to calculate the transition probabilities, so please select `y`:
  ```
  compute transition probabilities (E2/M1/E1) for
      Ne20_usda ? Y/N (default: Y) : y
  ```
  Now you may repeat the process and add parameters for another nuclide (per 2022-08-29 only one nuclide is supported at a time so the program will crash if you try to add an additional nuclide). Press enter to skip this step. For the last step you are asked if you want to split the commands into separate shell scripts. This is handy for running very large calculations on supercomputers, but not for running smaller calculations on single PCs. We'll choose `n`: 
  ```
  Split shell files? y/n (default: n): n
  Setup complete. Exiting...
  ```
  At this point the preparations before running the actual calculations are complete. Your data directory `~/kshell_results/ne20` should now contain these files:

  ```
  Ne20_usda.sh
  Ne20_usda_p.ptn
  collect_logs.py
  count_dim.py
  kshell.exe
  save_input_ui.txt
  transit.exe
  usda.snt
  ```
  See file descriptions later in this document if you want to know what they all do. Run your 20Ne `KSHELL` calculation by:
  ```
  ./Ne20_usda.sh
  ```
  If the program runs successfully, you will see:
  ```
  start running log_Ne20_usda_m0p.txt ...
  start running log_Ne20_usda_tr_m0p_m0p.txt ...
  Finish computing Ne20_usda.
  ```
  Congrats! You have just performed your first `KSHELL` calculation! To use these results, please se the *General usage* section later in this document.

  </p>
  </details>

  <!-- #### How to choose spin and parity states -->

  <details>
  <summary>How to choose spin and parity levels</summary>
  <p>
  
  `kshell_ui.py` asks you to choose what spin and parity levels you want to calculate:
  ```
  J, parity, number of lowest levels
    (ex. 100          for 100 +parity, 100 -parity levels w/o J-proj. (default)
        -5           for lowest five -parity levels,
        0+3, 2+1     for lowest three 0+ levels and one 2+ levels,
        1.5-1, 3.5+3 for lowest one 3/2- levels and three 7/2+ levels)
        range        for a range of levels) :
  ```
  * Entering an integer `N` will ask `KSHELL` to produce the `N` lowest lying energy levels, regardless of spin and parity. Example: Inputting `1337` will produce the 1337 lowest lying energy levels.
  * Prepending a plus sign (`+`) or a minus sign (`-`) to the integer will specify which parity you want to calculate the levels for. Note that your chosen nuclide and model space might only be able to produce either positive or negative parity levels. Example: `+1337` will produce the 1337 lowest lying positive parity levels.
  * You can request the `N` lowest lying levels of a specific spin and parity. Example: `0+3` will produce the three lowest lying levels with spin 0 and positive parity.
  * You can request several different specific spin and parity levels. Example: `1.5-1, 3.5+3` will produce the lowest lying level of spin 3/2 and negative parity, as well as the three lowest lying levels of spin 7/2 and positive parity.
  * The `range` functionality lets you easily select `N` levels for a range of different angular momenta. Any invalid choice will be filtered away, like choosing integer angular momenta for a nucleus of half integer angular momenta, or choosing a parity which the interaction does not support. In the following example we select 10 positive parity levels for angular momenta 0, 1, 2, and 3:
    
  ```
range
Start spin: 0
End spin (included): 3
Number of states per spin: 10
Parity (+, -, both): +
Chosen states: ['0+10', '1+10', '2+10', '3+10']
  ```

  </p>
  </details>

  <!-- #### How to calculate the dimensionality -->

  <details>
  <summary>How to calculate the dimensionality</summary>
  <p>

  After answering all the questions from `kshell_ui.py` it might be reasonable to check the dimensionality of the configuration to see if your computer will actually manage to solve the calculations. At this point, the results directory will look something like this:
  ```
  Ne20_usda.sh
  Ne20_usda_p.ptn
  collect_logs.py
  count_dim.py
  kshell.exe
  save_input_ui.txt
  transit.exe
  usda.snt
  ```
  The `.snt` file contains the two-body matrix elements (TBME) in the current model space (here `usda`). The `.ptn` contains the possible different proton and neutron combinations. Count the dimensionality by:
  ```
  python count_dim.py usda.snt Ne20_usda_p.ptn
  ```
  or by (you will be prompted for available `.snt` and `ptn` files)
  ```
  python count_dim.py
  ```
  which generates the output
  ```
        2*M        M-scheme dim.          J-scheme dim.
  dim.    16                    4                    4   4.00x10^ 0  4.00x10^ 0
  dim.    14                   16                   12   1.60x10^ 1  1.20x10^ 1
  dim.    12                   52                   36   5.20x10^ 1  3.60x10^ 1
  dim.    10                  116                   64   1.16x10^ 2  6.40x10^ 1
  dim.     8                  225                  109   2.25x10^ 2  1.09x10^ 2
  dim.     6                  354                  129   3.54x10^ 2  1.29x10^ 2
  dim.     4                  497                  143   4.97x10^ 2  1.43x10^ 2
  dim.     2                  594                   97   5.94x10^ 2  9.70x10^ 1
  dim.     0                  640                   46   6.40x10^ 2  4.60x10^ 1
  ```
  The M- and J-scheme dimensionalities are both very small in this configuration and the calculations will take only a few seconds to run on a normal laptop. The J-scheme dimensionality tells us how many levels of the different spins are available. From the above table we read that this configuration has 46 possible spin 0 states, 97 spin 1 states, 143 spin 2 states, and so on. We can also read from the table that this configuration has 640 possible M = 0 states (projection of J on the z-axis), 594 M = 1 states, and so on. The two last columns displays the M- and J-scheme dimensionalities in scientific notation.

  We now look at a much larger configuration, namely V50 with the `GXPF` model space:
  ```
  python count_dim.py gxpf1a.snt V50_gxpf1a_p.ptn
  ```
  gives:
  ```
        2*M        M-scheme dim.          J-scheme dim.
  dim.    44                    4                    4   4.00x10^ 0  4.00x10^ 0
  dim.    42                   46                   42   4.60x10^ 1  4.20x10^ 1
  dim.    40                  263                  217   2.63x10^ 2  2.17x10^ 2
  dim.    38                 1069                  806   1.07x10^ 3  8.06x10^ 2
  dim.    36                 3489                 2420   3.49x10^ 3  2.42x10^ 3
  dim.    34                 9737                 6248   9.74x10^ 3  6.25x10^ 3
  dim.    32                23975                14238   2.40x10^ 4  1.42x10^ 4
  dim.    30                53304                29329   5.33x10^ 4  2.93x10^ 4
  dim.    28               108622                55318   1.09x10^ 5  5.53x10^ 4
  dim.    26               205136                96514   2.05x10^ 5  9.65x10^ 4
  dim.    24               362005               156869   3.62x10^ 5  1.57x10^ 5
  dim.    22               600850               238845   6.01x10^ 5  2.39x10^ 5
  dim.    20               942669               341819   9.43x10^ 5  3.42x10^ 5
  dim.    18              1403670               461001   1.40x10^ 6  4.61x10^ 5
  dim.    16              1990227               586557   1.99x10^ 6  5.87x10^ 5
  dim.    14              2694122               703895   2.69x10^ 6  7.04x10^ 5
  dim.    12              3489341               795219   3.49x10^ 6  7.95x10^ 5
  dim.    10              4331494               842153   4.33x10^ 6  8.42x10^ 5
  dim.     8              5160580               829086   5.16x10^ 6  8.29x10^ 5
  dim.     6              5907365               746785   5.91x10^ 6  7.47x10^ 5
  dim.     4              6502475               595110   6.50x10^ 6  5.95x10^ 5
  dim.     2              6886407               383932   6.89x10^ 6  3.84x10^ 5
  dim.     0              7019100               132693   7.02x10^ 6  1.33x10^ 5
  ```
  The `GXPF` model space uses the 0f7/2, 1p3/2, 0f5/2 and 1p1/2 orbitals for the valence nucleons. V50 has 3 valence protons and 7 valence neutrons free to move about in the model space. Compared to 20Ne in the `USD` model space, V50 has both more valence nucleons and more states for them to be in, thus the larger M- and J-scheme dimensionalities. The V50 `GXPF` configuration might be possible to run on a multicore laptop for a small number of requested states. Running the configuration for the 100 lowest lying states for spins 0 to 14 takes approximately 1-2 hours on the Fram supercomputer using 32 nodes.

  </p>
  </details>

  <!-- #### How to truncate -->

  <details>
  <summary>How to truncate</summary>
  <p>
  
  #### Particle-hole truncation
  `kshell_ui.py` asks you if you want to truncate the model space. For large configurations (many valence nucleons and many shells for them to occupy) truncation might be necessary for `KSHELL` to actually complete the calculations. We use V50 in the `GXPF` model space as an example. This configuration has a dimensionality of (see above section on how to calculate the dimensionality):
  ```
        2*M        M-scheme dim.          J-scheme dim.
  dim.    44                    4                    4   4.00x10^ 0  4.00x10^ 0
  dim.    42                   46                   42   4.60x10^ 1  4.20x10^ 1
  dim.    40                  263                  217   2.63x10^ 2  2.17x10^ 2
  dim.    38                 1069                  806   1.07x10^ 3  8.06x10^ 2
  dim.    36                 3489                 2420   3.49x10^ 3  2.42x10^ 3
  dim.    34                 9737                 6248   9.74x10^ 3  6.25x10^ 3
  dim.    32                23975                14238   2.40x10^ 4  1.42x10^ 4
  dim.    30                53304                29329   5.33x10^ 4  2.93x10^ 4
  dim.    28               108622                55318   1.09x10^ 5  5.53x10^ 4
  dim.    26               205136                96514   2.05x10^ 5  9.65x10^ 4
  dim.    24               362005               156869   3.62x10^ 5  1.57x10^ 5
  dim.    22               600850               238845   6.01x10^ 5  2.39x10^ 5
  dim.    20               942669               341819   9.43x10^ 5  3.42x10^ 5
  dim.    18              1403670               461001   1.40x10^ 6  4.61x10^ 5
  dim.    16              1990227               586557   1.99x10^ 6  5.87x10^ 5
  dim.    14              2694122               703895   2.69x10^ 6  7.04x10^ 5
  dim.    12              3489341               795219   3.49x10^ 6  7.95x10^ 5
  dim.    10              4331494               842153   4.33x10^ 6  8.42x10^ 5
  dim.     8              5160580               829086   5.16x10^ 6  8.29x10^ 5
  dim.     6              5907365               746785   5.91x10^ 6  7.47x10^ 5
  dim.     4              6502475               595110   6.50x10^ 6  5.95x10^ 5
  dim.     2              6886407               383932   6.89x10^ 6  3.84x10^ 5
  dim.     0              7019100               132693   7.02x10^ 6  1.33x10^ 5
  ```
  which is too large to run on a regular computer for any decent amount of requested states. Lets see how the dimensionality changes with truncation. When `kshell_ui.py` asks for truncation, enter `1` to apply particle-hole truncation:

  ```
  truncation for "+" parity state in  V50_gxpf1a_p.ptn
  truncation scheme ?
        0 : No truncation (default)
        1 : particle-hole truncation for orbit(s)
        2 : hw truncation
        3 : Both (1) and (2)

  1
  ```
  which outputs:
  ```
    #    n,  l,  j, tz,    spe
    1    0   3   7  -1    -8.624     p_0f7/2
    2    1   1   3  -1    -5.679     p_1p3/2
    3    0   3   5  -1    -1.383     p_0f5/2
    4    1   1   1  -1    -4.137     p_1p1/2
    5    0   3   7   1    -8.624     n_0f7/2
    6    1   1   3   1    -5.679     n_1p3/2
    7    0   3   5   1    -1.383     n_0f5/2
    8    1   1   1   1    -4.137     n_1p1/2
  specify # of orbit(s) and min., max. occupation numbers for restriction

  # of orbit(s) for restriction?  (<CR> to quit):
  ```
  Here we see the valence orbitals 0f7/2, 1p3/2, 0f5/2 and 1p1/2, for both protons and neutrons. The `l` column denotes the angular momentum of the orbital, `j` the total angular momentum of the orbital, and `tz` the isospin. Let us now restrict the number of protons and neutrons allowed in the 0f7/2 orbital. In the above table we can see that the 0f7/2 orbitals are labeled 1 (protons) and 5 (neutrons). Set the maximum number of protons and neutrons to 2 in those orbitals by:
  ```
  # of orbit(s) for restriction?  (<CR> to quit): 1,5
  min., max. restricted occupation numbersfor the orbit(s) (or max only) : 2
  ```
  We now check the dimensionality of the truncated configuration:
  ```
        2*M        M-scheme dim.          J-scheme dim.
  dim.    36                    5                    5   5.00x10^ 0  5.00x10^ 0
  dim.    34                   58                   53   5.80x10^ 1  5.30x10^ 1
  dim.    32                  303                  245   3.03x10^ 2  2.45x10^ 2
  dim.    30                 1148                  845   1.15x10^ 3  8.45x10^ 2
  dim.    28                 3474                 2326   3.47x10^ 3  2.33x10^ 3
  dim.    26                 8930                 5456   8.93x10^ 3  5.46x10^ 3
  dim.    24                20129                11199   2.01x10^ 4  1.12x10^ 4
  dim.    22                40732                20603   4.07x10^ 4  2.06x10^ 4
  dim.    20                75106                34374   7.51x10^ 4  3.44x10^ 4
  dim.    18               127691                52585   1.28x10^ 5  5.26x10^ 4
  dim.    16               201896                74205   2.02x10^ 5  7.42x10^ 4
  dim.    14               298865                96969   2.99x10^ 5  9.70x10^ 4
  dim.    12               416333               117468   4.16x10^ 5  1.17x10^ 5
  dim.    10               547983               131650   5.48x10^ 5  1.32x10^ 5
  dim.     8               683573               135590   6.84x10^ 5  1.36x10^ 5
  dim.     6               810023               126450   8.10x10^ 5  1.26x10^ 5
  dim.     4               913390               103367   9.13x10^ 5  1.03x10^ 5
  dim.     2               981186                67796   9.81x10^ 5  6.78x10^ 4
  dim.     0              1004814                23628   1.00x10^ 6  2.36x10^ 4
  ```
  where we see that the dimensionality has been reduced by up to an order of magnitude for some spins.
    
  #### hw truncation
  Some interactions (model spaces), like sdpf-mu, span over several major shells. In this case we can use hw (hbar omega) truncation to limit the number of particles which are allowed to cross the major shell gap. Lets use 44Sc as an example. 44Sc with the sdpf-mu interaction has a dimensionality so large that we are not even able to calculate the dimensionality, let alone perform the calculations. Here we need to use hw truncation to drastically reduce the size. Choose option 2 when you are prompted for truncation (or option 3 if you plan on using particle-hole truncation in addition to hw):
  ```
  truncation for "+" parity state in  Sc44_sdpf-mu_p.ptn
  truncation scheme ?
       0 : No truncation (default)
       1 : particle-hole truncation for orbit(s)
       2 : hw truncation
       3 : Both (1) and (2)

  2
  (min. and) max hw for excitation : 3
  lowest hw, maxhw  60 63
  generating partition file ............ done.
  ```
  In the case of 44Sc with sdpf-mu you will be prompted for truncation on the negative parity states too. This example uses the same truncation for both + and -. In this example, 3 particles are allowed to cross the major shell gaps which results in a dimensionality of
  ```
        2*M        M-scheme dim.          J-scheme dim.
  dim.    42                    8                    8   8.00x10^ 0  8.00x10^ 0
  dim.    40                   84                   76   8.40x10^ 1  7.60x10^ 1
  dim.    38                  513                  429   5.13x10^ 2  4.29x10^ 2
  dim.    36                 2250                 1737   2.25x10^ 3  1.74x10^ 3
  dim.    34                 7950                 5700   7.95x10^ 3  5.70x10^ 3
  dim.    32                23800                15850   2.38x10^ 4  1.58x10^ 4
  dim.    30                62464                38664   6.25x10^ 4  3.87x10^ 4
  dim.    28               146820                84356   1.47x10^ 5  8.44x10^ 4
  dim.    26               313940               167120   3.14x10^ 5  1.67x10^ 5
  dim.    24               617562               303622   6.18x10^ 5  3.04x10^ 5
  dim.    22              1127352               509790   1.13x10^ 6  5.10x10^ 5
  dim.    20              1922531               795179   1.92x10^ 6  7.95x10^ 5
  dim.    18              3079113              1156582   3.08x10^ 6  1.16x10^ 6
  dim.    16              4651003              1571890   4.65x10^ 6  1.57x10^ 6
  dim.    14              6648334              1997331   6.65x10^ 6  2.00x10^ 6
  dim.    12              9018026              2369692   9.02x10^ 6  2.37x10^ 6
  dim.    10             11633108              2615082   1.16x10^ 7  2.62x10^ 6
  dim.     8             14296260              2663152   1.43x10^ 7  2.66x10^ 6
  dim.     6             16760154              2463894   1.68x10^ 7  2.46x10^ 6
  dim.     4             18762983              2002829   1.88x10^ 7  2.00x10^ 6
  dim.     2             20072284              1309301   2.01x10^ 7  1.31x10^ 6
  dim.     0             20527802               455518   2.05x10^ 7  4.56x10^ 5
  ```
  and
  ```
        2*M        M-scheme dim.          J-scheme dim.
  dim.    50                    6                    6   6.00x10^ 0  6.00x10^ 0
  dim.    48                   95                   89   9.50x10^ 1  8.90x10^ 1
  dim.    46                  735                  640   7.35x10^ 2  6.40x10^ 2
  dim.    44                 3972                 3237   3.97x10^ 3  3.24x10^ 3
  dim.    42                16782                12810   1.68x10^ 4  1.28x10^ 4
  dim.    40                59228                42446   5.92x10^ 4  4.24x10^ 4
  dim.    38               181116               121888   1.81x10^ 5  1.22x10^ 5
  dim.    36               492378               311262   4.92x10^ 5  3.11x10^ 5
  dim.    34              1210949               718571   1.21x10^ 6  7.19x10^ 5
  dim.    32              2729673              1518724   2.73x10^ 6  1.52x10^ 6
  dim.    30              5695210              2965537   5.70x10^ 6  2.97x10^ 6
  dim.    28             11083379              5388169   1.11x10^ 7  5.39x10^ 6
  dim.    26             20241387              9158008   2.02x10^ 7  9.16x10^ 6
  dim.    24             34862609             14621222   3.49x10^ 7  1.46x10^ 7
  dim.    22             56856340             21993731   5.69x10^ 7  2.20x10^ 7
  dim.    20             88092886             31236546   8.81x10^ 7  3.12x10^ 7
  dim.    18            130029311             41936425   1.30x10^ 8  4.19x10^ 7
  dim.    16            183263256             53233945   1.83x10^ 8  5.32x10^ 7
  dim.    14            247098324             63835068   2.47x10^ 8  6.38x10^ 7
  dim.    12            319234048             72135724   3.19x10^ 8  7.21x10^ 7
  dim.    10            395690620             76456572   3.96x10^ 8  7.65x10^ 7
  dim.     8            471046277             75355657   4.71x10^ 8  7.54x10^ 7
  dim.     6            539002617             67956340   5.39x10^ 8  6.80x10^ 7
  dim.     4            593209277             54206660   5.93x10^ 8  5.42x10^ 7
  dim.     2            628208483             34999206   6.28x10^ 8  3.50x10^ 7
  dim.     0            640309604             12101121   6.40x10^ 8  1.21x10^ 7
  ```
  </p>
  </details>

  <!-- #### How to use the output from KSHELL -->

  <details>
  <summary>How to use the output from KSHELL</summary>
  <p>

  After running `KSHELL`, your work directory will look similar to this:
  ```
  Ne20_usda.sh
  Ne20_usda_m0p.wav
  Ne20_usda_p.ptn
  count_dim.py
  kshell.exe
  log_Ne20_usda_m0p.txt
  log_Ne20_usda_tr_m0p_m0p.txt
  save_input_ui.txt
  transit.exe
  usda.snt
  ```
  All the level data are located in `log_Ne20_usda_m0p.txt` and all the transition data are located in `log_Ne20_usda_tr_m0p_m0p.txt`.

  #### Load and view data from KSHELL

  The log files are easily read with the `kshell-utilities` package. See the docstrings in the [kshell-utilities repository](https://github.com/GaffaSnobb/kshell-utilities) for extended documentation. Install the package with `pip`:
  ```
  pip install kshell-utilities
  ```
  Create a blank Python file with your favourite editor. Lets name it `ne20.py` and lets place it in the results directory of the 20Ne calculation which is `~/kshell_results/ne20` according to this guide. However, the use of `~/` as a shortcut to your home directory is not standard in Python and is discouraged to be used, so if you wish to specify the path to your home directory, use the actual path. For macOS: `/Users/<your username>/kshell_results/ne20`. For most Linux distros: `/home/<your username>/kshell_results/ne20`. We use the `loadtxt` function to read the results from `KSHELL`:
  ``` python
  import kshell_utilities as ksutil

  def main():
    ne20 = ksutil.loadtxt(path=".")

  if __name__ == "__main__":
    main()
  ```
  The use of a name guard (`if __name__ == "__main__":`) is required because `kshell-utilities` uses Python's `multiprocessing` module which requires this to function properly. Note that `path` is a period (`.`). This simply means that the `KSHELL` results are located in the same directory as the Python file `ne20.py`. If we do not place `ne20.py` in the same directory as the `KSHELL` results, then we need to specify either the relative path to the `KSHELL` results from `ne20.py` or the absolute path of the `KSHELL` results which is (macOS) `/Users/<your username>/kshell_results/ne20`. Back to the `loadtxt` function. `ne20` is an instance containing several useful attributes. To see the available attributes:
  ``` python
  > print(ne20.help)
  ['debug',
  'fname_ptn',
  'fname_summary',
  'gamma_strength_function_average_plot',
  'gsf',
  'help',
  'level_density_plot',
  'level_plot',
  'levels',
  'model_space',
  'negative_spin_counts',
  'neutron_partition',
  'nucleus',
  'parameters',
  'path',
  'proton_partition',
  'transitions_BE1',
  'transitions_BE2',
  'transitions_BM1',
  'truncation']
  ```
  To see the energy, 2\*angular momentum and parity of each level:
  ``` python
  > print(ne20.levels)
  [[-40.467   0.      1.   ]
   [-38.771   4.      1.   ]
   [-36.376   8.      1.   ]
   [-33.919   0.      1.   ]
   [-32.882   4.      1.   ]
   [-32.107  12.      1.   ]
   ...
   [-25.978  12.      1.   ]
   [-25.904  10.      1.   ]
   [-25.834   8.      1.   ]
   [-25.829   2.      1.   ]]
  ```
  Slice the array to get only selected values, if needed (`ne20.levels[:, 0]` for only the energies). To see 2\*spin_initial, parity_initial, Ex_initial, 2\*spin_final, parity_final, Ex_final, E_gamma, B(.., i->f), B(.., f<-i)] for the M1 transitions:
  ``` python
  > print(ne20.transitions_BM1)
  [[4.0000e+00 1.0000e+00 1.6960e+00 ... 7.5850e+00 5.8890e+00 0.0000e+00]
  [4.0000e+00 1.0000e+00 1.6960e+00 ... 9.9770e+00 8.2810e+00 4.8200e-01]
  [4.0000e+00 1.0000e+00 7.5850e+00 ... 9.9770e+00 2.3920e+00 1.1040e+00]
  ...
  [4.0000e+00 1.0000e+00 1.3971e+01 ... 1.4638e+01 6.6700e-01 6.0000e-03]
  [0.0000e+00 1.0000e+00 1.4126e+01 ... 1.4638e+01 5.1200e-01 2.0000e-02]
  [2.0000e+00 1.0000e+00 1.4336e+01 ... 1.4638e+01 3.0200e-01 0.0000e+00]]
  ```

  #### Visualise data from KSHELL 
  ##### Create a level density plot

  You can easily create a level density plot by
  ``` python
  ne20.level_density_plot()
  ```
  An alternative way is:
  ``` python
  ground_state_energy: float = ne20.levels[0, 0]
  ksutil.level_density(
      levels = ne20.levels[:, 0] - ground_state_energy,
      bin_width = 0.2,
      plot = True
  )
  ```
  Note that scaling the excitation energies by the ground state energy is required with this method. If you want greater control of `matplotlib` plotting parameters, use this method:
  ``` python
  import matplotlib.pyplot as plt

  ground_state_energy: float = ne20.levels[0, 0]
  bins, density = ksutil.level_density(
      levels = ne20.levels[:, 0] - ground_state_energy,
      bin_width = 0.2,
      plot = False,
  )
  plt.step(bins, density)
  plt.show()
  ```
  The `bin_width` is in the same energy units as your results, which for `KSHELL` is MeV. The two latter ways of generating the plot does not require that the data comes from `KSHELL`. Use any energy level data normalised to the ground state energy. The plot will look like this:
  
  <details>
  <summary>Click to see level density plot</summary>
  <p>

  ![level_density_plot](https://github.com/GaffaSnobb/kshell-utilities/blob/main/doc/level_density_plot_ne20.png)

  </p>
  </details>

  ##### Create a level plot / level scheme

  To generate a level plot:
  ``` python
  ne20.level_plot()
  ```
  or
  ``` python
  import matplotlib.pyplot as plt

  fig, ax = plt.subplots()
  ksutil.level_plot(
      levels = ne20.levels,
      ax = ax
  )
  plt.show()
  ```

  <details>
  <summary>Click to see level plot</summary>
  <p>

  ![level_plot](https://github.com/GaffaSnobb/kshell-utilities/blob/main/doc/level_plot_ne20.png)

  </p>
  </details>

  Both ways of generating the level plot supports selecting what total angular momenta to include in the plot, and how many levels per angular momentum. 
  ``` python
  ne20.level_plot(
      include_n_levels = 3,
      filter_spins = [0, 3, 5]
  )
  ```

  <details>
  <summary>Click to see filtered level plot</summary>
  <p>

  ![filtered_level_plot](https://github.com/GaffaSnobb/kshell-utilities/blob/main/doc/level_plot_filtered_ne20.png)

  </p>
  </details>

  ##### Create a gamma strength function plot
  
  The gamma strengh function (averaged over total angular momenta and parities) can easily be calculated in several ways. The quickest way is
  ``` python
    ne20.gsf()
  ```
  which is an alias for the following function call:
  ``` python
    ne20.gamma_strength_function_average_plot(
        bin_width = 0.2,
        Ex_max = 5,
        Ex_min = 20,
        multipole_type = "M1",
        plot = True,
        save_plot = False
    )
  ```
  The default parameters are applied if no function arguments are supplied. If you want to have greater control over the plotting procedure, then this solution is better:
  ``` python
    import matplotlib.pyplot as plt
    
    bins, gsf = ne20.gamma_strength_function_average_plot(
        bin_width = 0.2,
        Ex_max = 50,
        Ex_min = 5,
        multipole_type = "M1",
        plot = False,
        save_plot = False
    )
    plt.plot(bins, gsf)
    plt.show()
  ```
  since you yourself have control over the `matplotlib` calls. Note that `Ex_max` is set to way higher energy than you get from the `KSHELL` calculations. Typical max energy from a `KSHELL` calculation is in the range `[8, 12]`MeV. The default upper limit is set large as to include all levels of any `KSHELL` calculation. The final way of doing it is:
  ``` python
  import matplotlib.pyplot as plt

  bins, gsf = ksutil.gamma_strength_function_average(
      levels = ne20.levels,
      transitions = ne20.transitions_BM1,
      bin_width = 0.2,
      Ex_min = 5,
      Ex_max = 20,
      multipole_type = "M1"
  )
  plt.plot(bins, gsf)
  plt.show()
  ```
  where the difference is that you supply the `levels` and `transitions` arrays. I'd not recommend this final solution unless you have level and transition data from some other place than `KSHELL`. The parameters `bin_width`, `Ex_max` and `Ex_min` are in the same unit as the input energy levels, which from `KSHELL` is in MeV. `bin_width` is the width of the bins when the level density is calculated. `Ex_min` and `Ex_max` are the lower and upper limits for the excitation energy of the initial state of the transitions.

  <details>
  <summary>Click to see gamma strength function plot</summary>
  <p>

  ![gsf_plot](https://github.com/GaffaSnobb/kshell-utilities/blob/main/doc/gsf_ne20.png)

  </p>
  </details>

  </p>
  </details>

  <details>
  <summary>A more advanced example of using output from KSHELL</summary>
  <p>

  ##### Acquire some beefy 44Sc results

  For these more advanced examples, we need beefier files than the previous 20Ne example. Lets use a scandium-44 calculation I performed for my master's thesis. Start by creating a new directory for the 44Sc results in your results directory:
  ``` bash
  cd ~/kshell_results
  mkdir sc44
  cd sc44
  ```
  Then, copy the three files `000_Sc44_GCLSTsdpfsdgix5pn_tr_j0p_j2p.sh`, `save_input_ui.txt`, and `summary_Sc44_GCLSTsdpfsdgix5pn_000.tgz` from [here](https://github.com/GaffaSnobb/master-tasks/tree/main/Sc44/sdpf-sdg/200_levels/3hw) to the `sc44` directory you just created. The `.tgz` file has a download button, but the `.sh` and `.txt` files you have to copy-paste. While in the directory `~/kshell_results/sc44`, create these files with your favourite editor, for example VSCode, by:
  ``` bash
  code save_input_ui.txt
  code 000_Sc44_GCLSTsdpfsdgix5pn_tr_j0p_j2p.sh
  ```
  and copy-paste [the contents for the .sh file](https://github.com/GaffaSnobb/master-tasks/blob/main/Sc44/sdpf-sdg/200_levels/3hw/000_Sc44_GCLSTsdpfsdgix5pn_tr_j0p_j2p.sh) and [the contents for the .txt file](https://github.com/GaffaSnobb/master-tasks/blob/main/Sc44/sdpf-sdg/200_levels/3hw/save_input_ui.txt) to their respective files which you just created, and be sure to save the files.


  The `.tgz` file contains the 44Sc results from `KSHELL`, but the file must be un-compressed before it can be used by `kshell-utilities`. Still in the `~/kshell_results/sc44` directory, run the command
  ``` bash
  tar -xzvf summary_Sc44_GCLSTsdpfsdgix5pn_000.tgz
  ```
  to un-compress the file. You now have another file, `summary_Sc44_GCLSTsdpfsdgix5pn_000.txt`, in the same directory! Great! We need to download one more file which you can find [here](https://github.com/GaffaSnobb/master-tasks/blob/main/Sc44/Sc44_gsf.txt) (it has a download button). This file contains the experimental gamma strength function of 44Sc. Please place it in the same directory, namely `~/kshell_results/sc44`.

  ##### Load the 44Sc data into kshell-utilities

  While in the directory `~/kshell_results/sc44`, create a Python script named `sc44.py` and read the newly un-compressed summary file by:
  
  ```python
  import kshell_utilities as ksutil

  def main():
    sc44 = ksutil.loadtxt(
      path = "summary_Sc44_GCLSTsdpfsdgix5pn_000.txt"
    )

  if __name__ == "__main__":
      main()
  ```
  This summary file is quite large and will take 10-30 seconds to load. Your terminal should look like this when the process is done:

  ```bash
  > python sc44.py
  Thread 0 loading Energy values...
  Thread 1 loading B(E1) values...
  Thread 2 loading B(M1) values...
  Thread 3 loading B(E2) values...
  Thread 0 finished loading Energy values in 0.03 s
  Thread 2 finished loading B(M1) values in 6.43 s
  Thread 1 finished loading B(E1) values in 6.61 s
  Thread 3 finished loading B(E2) values in 10.51 s
  ```
  
  Note that your `~/kshell_results/sc44` directory now has a new directory called `tmp`. This new directory contains the `KSHELL` data from the summary file stored as binary `numpy` arrays. If you run `sc44.py` again, you will se that the output is different and that the program uses 1-2 seconds instead of 10-30 seconds to run:
  
  ```bash
  > python sc44.py
  Summary data loaded from .npy! Use loadtxt parameter load_and_save_to_file = 'overwrite' to re-read data from the summary file.
  ```

  Instead of re-reading the data from the summary text file, `kshell-utilities` now loads the binary `numpy` arrays which is much faster. You may at any time delete the `tmp` directory without losing any data. The only downside is that the next time you run the program it will use some time reading the summary text file again. The reason to include the `save_input_ui.txt` and `000_Sc44_GCLSTsdpfsdgix5pn_tr_j0p_j2p.sh` files is because they contain specific information about the calculation parameters of the 44Sc calculations, like the number of levels per angular momentum-parity pair, the truncation, the interaction used, etc. `kshell-utilities` uses this information to generate unique identifiers for the contents of the `tmp` directory in case the `tmp` directory should contain data from several different 44Sc `KSHELL` calculations. Not strictly necessary for this example, but this is the intended way to use `kshell-utilities`.


  ##### Take a look at the gamma strength function of 44Sc
  Lets take a look at a properly calculated gamma strength function. For reference, this 44Sc calculation took a few days of calculation time on [Betzy, Norway's most powerful supercomputer](https://documentation.sigma2.no/hpc_machines/betzy.html). Extend your Python script to include the following:

  ```python
  import kshell_utilities as ksutil
  import numpy as np
  import matplotlib.pyplot as plt
  ksutil.latex_plot()
  ksutil.flags["debug"] = True

  BIN_WIDTH = 0.2
  EX_MIN = 5
  EX_MAX = 9.699    # S(n).

  def main():
    fig, ax = plt.subplots()
    N, Ex, gsf_experimental, gsf_std = np.loadtxt("Sc44_gsf.txt", skiprows=2, unpack=True)
    ax.errorbar(Ex, gsf_experimental, yerr=gsf_std, fmt=".", capsize=1, elinewidth=0.5, label="Exp", color="black")
    
    sc44 = ksutil.loadtxt(
      path = "summary_Sc44_GCLSTsdpfsdgix5pn_000.txt",
    )
    bins, gsf_M1 = sc44.gsf(
      bin_width = BIN_WIDTH,
      Ex_min = EX_MIN,
      Ex_max = EX_MAX,
      multipole_type = "M1",
      plot = False
    )
    bins, gsf_E1 = sc44.gsf(
      bin_width = BIN_WIDTH,
      Ex_min = EX_MIN,
      Ex_max = EX_MAX,
      multipole_type = "E1",
      plot = False
    )
    ax.step(bins, (gsf_M1 + gsf_E1), label=r"SM $E1 + M1$", color="grey")
    ax.step(bins, gsf_M1, label=r"SM $M1$", color="red")
    ax.step(bins, gsf_E1, label=r"SM $E1$", color="blue")

    ax.set_yscale('log')
    ax.set_xlabel(r"E$_{\gamma}$ [MeV]")
    ax.set_ylabel(r"GSF [MeV$^{-3}$]")
    ax.legend()
    plt.show()
  ```
  The function call `ksutil.latex_plot()` makes your plots look nicer by making it "Latex style", whatever that means. Well, it actually means changing a few fonts and sizes, and you can see [the exact code here](https://github.com/GaffaSnobb/kshell-utilities/blob/ab0d7f9b261692a412d50508c6c66349f7208862/kshell_utilities/parameters.py#L11). It looks much prettier than default `matplotlib` and it fits right into your thesis. The line `ksutil.flags["debug"] = True` makes `kshell-utilities` be more verbose and can help you resolve issues. If you ever get tired of the terminal output you can set it to `False`.
  
  
  
  Run `sc44.py` again now and let it think for a few seconds. You should see a bunch of debug information like so:
  
  <details>
  <summary>Click here to see a bunch of debug information</summary>
  <p>

  ```bash
  > python sc44.py
  Summary data loaded from .npy! Use loadtxt parameter load_and_save_to_file = 'overwrite' to re-read data from the summary file.
  loadtxt_time = 0.1195173840096686 s
  --------------------------------
  transit_gsf_time = 0.770901508978568 s
  level_density_gsf_time = 0.0019428109808359295 s
  gsf_time = 0.0072257140127476305 s
  avg_gsf_time = 8.191799861378968e-05 s
  total_gsf_time = 0.7898409570043441 s
  multipole_type = 'M1'
  Skips: Transit: Energy range: 698614
  Skips: Transit: Number of levels: 0
  Skips: Transit: Parity: 0
  Skips: Level density: Energy range: 2320
  Skips: Level density: Number of levels: 0
  Skips: Level density: Parity: 0
  transit_total_skips = 698614
  n_transitions = 898504
  n_transitions_included = 199890
  level_density_total_skips = 2320
  n_levels = 3600
  n_levels_included = 1280
  --------------------------------
  --------------------------------
  transit_gsf_time = 0.44835006099310704 s
  level_density_gsf_time = 0.0018729210132732987 s
  gsf_time = 0.007104302989318967 s
  avg_gsf_time = 7.715600077062845e-05 s
  total_gsf_time = 0.4653512270015199 s
  multipole_type = 'E1'
  Skips: Transit: Energy range: 879173
  Skips: Transit: Number of levels: 0
  Skips: Transit: Parity: 0
  Skips: Level density: Energy range: 2320
  Skips: Level density: Number of levels: 0
  Skips: Level density: Parity: 0
  transit_total_skips = 879173
  n_transitions = 958400
  n_transitions_included = 79227
  level_density_total_skips = 2320
  n_levels = 3600
  n_levels_included = 1280
  --------------------------------
  ```
  </p>
  </details>

  and you'll se a nice GSF plot with both experimental values and `KSHELL` calculations. But, hold on... There is something strange about this plot...


  BREAKING NEWS: Ola Nordmann (43) was SHOCKED when he discovered why there is such a big difference between the experimental data and the calculated GSF of 44Sc. [Read the full story here!](https://github.com/GaffaSnobb/master-tasks/blob/main/doc/masters_thesis_final.pdf)


  Note that when you run `sc44.py` again it is much faster than the first run. If you peek inside the `tmp` directory you'll see that there are now additional files there. The GSF has been stored as binary `numpy` arrays and it does not have to be re-calculated during subsequent runs of the program. This means that you can make all your millions of tiny plot adjustments without waiting for a long time to show the changes. Neat, eh? Note also that if you change any of the parameters of the GSF, like `bin_width`, `Ex_min` and `Ex_max`, then `kshell-utilities` will understand that this is a different calculation from your previous one and it will perform new calculations and save these as binary `numpy` arrays too. These saved `.npy` GSF files only take up a few hundred bytes so don't worry about storing many different calculations (the saved `.npy` files of the transition calculations from `KSHELL` however can take several hundred megabytes but these are only generated once per `KSHELL` calculation).

  ```
  > python sc44.py
  Summary data loaded from .npy! Use loadtxt parameter load_and_save_to_file = 'overwrite' to re-read data from the summary file.
  loadtxt_time = 0.1136299180216156 s
  Sc44 M1 GSF data loaded from .npy!
  Sc44 E1 GSF data loaded from .npy!
  ```

  ##### Level density as a function of energy, angular momentum, and parity

  Let's look at some other fancy stuff, shall we? When I made the following functionality my intentions were to study how the total angular momentum distribution looked like with regards to energy. The result however turned out to be a level density heatmap where the level density is plotted as a function of total angular momentum, energy, and parity. Still a nice result, but the name of the function is a bit off. Add this to your code and run it:

  ```python
  sc44.angular_momentum_distribution_plot(
    bin_width = 1,
    E_min = EX_MIN,
    E_max = 15,
    filter_parity = "+",
    save_plot = False,
    # j_list = [0, 2, 4, 7]
  )
  ```

  You can specify a selection of total angular momenta with the `j_list` parameter. Note that `E_min` and `E_max` do not mean the exact same thing as `Ex_min` and `Ex_max`. 

  ##### What effect does the number of levels have?
  Have you ever lay awake at night, wondering about what the hell would happen if you changed the number of levels per angular momentum and per parity in your `KSHELL` calculations? Me too! Lets stop wondering:
  
  ```python
  n_levels = [60, 100, 200]
  colors = ["cyan", "dodgerblue", "blue"]
  fig_0, ax_0 = plt.subplots()
  fig_1, ax_1 = plt.subplots()

  for levels, color in zip(n_levels, colors):
    bins_gsf_E1, gsf_E1 = sc44.gsf(
      bin_width = BIN_WIDTH,
      Ex_min = EX_MIN,
      Ex_max = EX_MAX,
      multipole_type = "E1",
      include_n_levels = levels,
      plot = False
    )
    ax_0.step(bins_gsf_E1, gsf_E1, label=f"{levels} levels per " + r"$j^{\pi}$", color=color)

    bins_nld, nld = sc44.nld(
      bin_width = BIN_WIDTH,
      include_n_levels = levels,
      plot = False
    )
    ax_1.step(bins_nld, nld, label=f"{levels} levels per " + r"$j^{\pi}$", color=color)
  
  ax_0.set_yscale('log')
  ax_0.set_xlabel(r"$E_{\gamma}$ [MeV]")
  ax_0.set_ylabel(r"GSF [MeV$^{-3}$]")
  ax_0.legend()

  ax_1.set_xlabel(r"$E$ [MeV]")
  ax_1.set_ylabel(r"NLD [MeV$^{-1}$]")
  ax_1.legend()
  plt.show()
  ```

  ##### A small generalised Brink-Axel test
  This one is a real treat! Do you wonder if the gBA holds for your `KSHELL` calculations? This might shed some light on the matter:

  ```python
  fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6.4, 4.8*2))
  j_list = [0, 1, 2, 3, 4, 5, 6, 7, 8]
  n_j = len(j_list)

  bins_M1_all_j, gsf_M1_all_j = sc44.gsf(
    bin_width = BIN_WIDTH,
    Ex_min = EX_MIN,
    Ex_max = EX_MAX,
    multipole_type = "M1",
    plot = False,
  )
  bins_E1_all_j, gsf_E1_all_j = sc44.gsf(
    bin_width = BIN_WIDTH,
    Ex_min = EX_MIN,
    Ex_max = EX_MAX,
    multipole_type = "E1",
    plot = False,
  )
  ax[0].plot(bins_M1_all_j, gsf_M1_all_j, color="black", label=r"All $j_i$")
  ax[1].plot(bins_E1_all_j, gsf_E1_all_j, color="black", label=r"All $j_i$")

  for i in range(n_j):
    bins_M1_one_j, gsf_M1_one_j = sc44.gsf(
      bin_width = BIN_WIDTH,
      Ex_min = EX_MIN,
      Ex_max = EX_MAX,
      multipole_type = "M1",
      partial_or_total = "partial",
      filter_spins = [j_list[i]],
      plot = False,
    )
    bins_E1_one_j, gsf_E1_one_j = sc44.gsf(
      bin_width = BIN_WIDTH,
      Ex_min = EX_MIN,
      Ex_max = EX_MAX,
      multipole_type = "E1",
      partial_or_total = "partial",
      filter_spins = [j_list[i]],
      plot = False,
    )
    ax[0].plot(bins_M1_one_j, gsf_M1_one_j, color="black", alpha=0.2)
    ax[1].plot(bins_E1_one_j, gsf_E1_one_j, color="black", alpha=0.2)

  ax[0].set_title(r"$^{44}$Sc, $M1$")
  ax[0].set_yscale("log")
  ax[0].set_ylabel(r"GSF [MeV$^{-3}$]")
  ax[0].plot([0], [0], color="black", alpha=0.2, label=r"Single $j_i$")  # Dummy for legend.
  ax[0].legend(loc="lower left")

  ax[1].set_title(r"$^{44}$Sc, $E1$")
  ax[1].set_yscale("log")
  ax[1].set_xlabel(r"$E_{\gamma}$ [MeV]")
  ax[1].set_ylabel(r"GSF [MeV$^{-3}$]")
  ax[1].plot([0], [0], color="black", alpha=0.2, label=r"Single $j$")  # Dummy for legend.
  ax[1].legend(loc="lower left")
  plt.show()
  ```
  For this one it is really nice that `kshell-utilities` saves the GSFs as `.npy` because you need like 18 of them to generate the plots. Run it once more and BAM! The plots show up instantly.

  ##### The Porter-Thomas distribution

  We can't mention gBA without mentioning the Porter-Thomas distribution. The following code will plot a histogram of B values (reduced transition probabilities) from selections of Ei values (thanks to Jørgen Midtbø for creating the figure from which the following is heavily inspired):

  ```python
  sc44.porter_thomas_Ei_plot(
    Ei_range_min = EX_MIN,
    Ei_range_max = EX_MAX,
    Ei_values = np.linspace(EX_MIN, EX_MAX, 3),
    Ei_bin_width = 0.2,
    BXL_bin_width = 0.1,
    multipole_type = "M1",
  )
  ```

  [The docstring of this function](https://github.com/GaffaSnobb/kshell-utilities/blob/ab0d7f9b261692a412d50508c6c66349f7208862/kshell_utilities/kshell_utilities.py#L743) explains in detail what all the parameters are. A similar plot but analysed for total angular momentum instead of excitation energy can be created by:

  ```python
  sc44.porter_thomas_j_plot(
    Ex_max = EX_MAX,
    Ex_min = EX_MIN,
    j_lists = [[0, 1, 2], [3, 4, 5], [6, 7, 8]],
  )
  ```

  and the parameters are described in [the docstring](https://github.com/GaffaSnobb/kshell-utilities/blob/ab0d7f9b261692a412d50508c6c66349f7208862/kshell_utilities/kshell_utilities.py#L1008).


  
  If you wonder what all of this stuff might mean, check out [my masters thesis](https://github.com/GaffaSnobb/master-tasks/blob/main/doc/masters_thesis_final.pdf).
  </p>
  </details>
  
  <!-- #### KSHELL file descriptions -->
    
  <details>
  <summary>Compare data from different KSHELL runs</summary>
  <p>

    It might be very useful to compare data from different KSHELL calculations. Maybe you want to compare the energy levels of 50V and 51V? We've got your back! Note that the `ksutil.latex_plot()` call is optional since it works only to make the plot prettier. If you want to use the exact same data as the example code, you can find it here: [50V](https://github.com/GaffaSnobb/master-tasks/tree/main/V50/sdpf-mu/200_levels/1hw), [51V](https://github.com/GaffaSnobb/master-tasks/tree/main/V51/sdpf-mu/200_levels)
    ```python
    import kshell_utilities as ksutil
    ksutil.latex_plot()

    def main():
        V50 = ksutil.loadtxt("summary_V50_sdpf-mu_000.txt")
        V51 = ksutil.loadtxt("summary_V51_sdpf-mu_000.txt")

        comp = ksutil.Compare(kshell_outputs=[V50, V51])
        comp.plot_levels()

    if __name__ == "__main__":
        main()
    ```
    
  <details>
  <summary>Click to see the level comparison plot produced by the above code</summary>
  <p>

  ![level_comparison_plot](https://github.com/GaffaSnobb/kshell-utilities/blob/main/doc/level_comparison_plot.png)

  </p>
  </details>
    
  </p>
  </details>

  <details>
  <summary>KSHELL file descriptions</summary>
  <p>

  #### .sh
  The `.sh` file(s) is (are) generated by `kshell_ui` and contain the run commands for `KSHELL`. This is the file you run to start `KSHELL`.
  #### .wav
  The `.wav` files are generated after running the `KSHELL` executable. They contain the eigenvectors of the Hamiltonian matrix and are used to compute the transition probabilities.
  #### .snt
  The `.snt` files contain the parameters for each of the interactions. For example `usda.snt`, `gxpf1a.snt` etc. They are located in `<install_directory>/snt` and after completing the `kshell_ui` setup, the chosen interaction file is copied to the run directory.
  #### .ptn
  The `.ptn` files are generated by `kshell_ui` and contain the possible proton and neutron configurations of the chosen model space and nucleus with the chosen truncation.
  #### .exe
  The `.exe` files are the compiled executable program files. These are generated by the compilation process and are located in `<install_directory>/src`. They are copied to the run directory after the `kshell_ui` setup.
  #### .input
  The `.input` files contain run parameters for the `.exe` files. They are deleted after a successful calculation. If you see such a file in your run directory after the program has terminated, then something went wrong during the calculation.
  #### log\_\*.txt
  There are log files for level information and separate log files for transition information. The log files contain all level and transition information from `KSHELL` in addition to debug parameters like RAM usage, time usage, and much more. If specific angular momenta were chosen during the `kshell_ui` setup, then there will be one level log file for each of the angular momentum choices; if you also chose to calculate transition probabilities then there will be one transition log file for each unique initial angular momentum and initial parity to final angular momentum and final parity pair. If you chose just a number of levels without specifying any angular momenta, then there will be only one level log file and one transition log file.
  #### summary\_\*.txt
  Older versions of `KSHELL` will generate a summary file after all calculations have been executed. The summary file contains all the nuclear data from the log files without any of the `KSHELL` debug information. Newer versions of `KSHELL` will however not generate such a summary file. `kshell_utilities` will generate a summary file if you point `kshell_utilities.loadtxt` to the directory of the log files. The summary file is not strictly needed for `kshell_utilities` to work since `kshell_utilities` can read the data from the log files. The summary is generated for legacy reasons.

  </p>
  </details>

          
## Pitfalls

<details>
<summary>Click here for pitfalls</summary>
<p>
  
  #### Crashes on Betzy
  ``` bash
  srun: error: b5272: task 10: Broken pipe
  [mpiexec@b1373.betzy.sigma2.no] wait_proxies_to_terminate (../../../../../src/pm/i_hydra/mpiexec/intel/i_mpiexec.c:527): downstream from host b1373 exited with status 141
  ```
  Try to increase or decrease `n_block`. This error occurred for me, crashing the job after just a few seconds, when calculating 200 1- states for 68Zn with gs8 (sdpf-sdg) with `n_block = 0`. Setting `n_block = 8` solved the problem.

  #### error [dcg]: invalid j or m
  KSHELL might raise this error, meaning that the projection `m` is larger than the angular momentum of the state, `j`. This error probably occurs in combination with using block Lanczos (`n_block = 8` for example). Setting `n_block = 0` should resolve this problem, though at an increase in computation time.
  
  #### Small M- / J-scheme dimensionalities on many cores
  If you are performing a calculation of relatively small dimensionality, be sure to not use too many CPU cores. This is not applicable to normal desktop / laptop computers, but to supercomputer with thousands of cores. Best case, the program crashes. Worst case, the program does nothing for the enitre duration of the allocated time. The program might run fine, but not using all the allocated resources and thus wasting CPU hours. As an example, Sc45 in sdpf-sdg with a max M-scheme dimensionality of 1.4e6 does not run well on 64 nodes on Betzy and just uses all of the allocated time doing nothing. Reducing the number of nodes to 4 solved the calculations in under 5 minutes. Dimensionality above 1e7 should work fine with any number of nodes on betzy.
  
  2021-09-29 UPDATE: `kshell_ui.py` now checks if the number of requested states exceeds the maximum possible number of states for the given model space and configuration and adjusts accordingly. This error should not be a problem anymore for single PC compilation. We still do experience this issue when compiled with `-DMPI`, but running KSHELL a with a small number of possible configurations on several nodes is nonsenical; reduce the number of nodes.

  KSHELL version 2 has undefined behavior if you request more states than the configuration and model space allows. As an example, take 28Ar in the USDA model space. By running the `count_dim.py` script we get
  ```
  python <path>/count_dim.py usda.snt Ar28_usda_p.ptn
        2*M        M-scheme dim.          J-scheme dim.
  dim.    16                    4                    4   4.00x10^ 0  4.00x10^ 0
  dim.    14                   16                   12   1.60x10^ 1  1.20x10^ 1
  dim.    12                   52                   36   5.20x10^ 1  3.60x10^ 1
  dim.    10                  116                   64   1.16x10^ 2  6.40x10^ 1
  dim.     8                  225                  109   2.25x10^ 2  1.09x10^ 2
  dim.     6                  354                  129   3.54x10^ 2  1.29x10^ 2
  dim.     4                  497                  143   4.97x10^ 2  1.43x10^ 2
  dim.     2                  594                   97   5.94x10^ 2  9.70x10^ 1
  dim.     0                  640                   46   6.40x10^ 2  4.60x10^ 1
  ```
  The `J-scheme dim.` column indicates how many different states of the spin indicated in the `2*M` column that can be calculated in this model space with this configuration of protons and neutrons. 28Ar in USDA has 10 valence protons and 2 valence neutrons, and from `count_dim.py` we see that this model space and configuration allows 46 0+ states, 97 1+ states, 143 2+ states, and so on. Take the 0+ states as an example. If you request more than 46 0+ states, say 100, the best case scenario is that KSHELL gives you 46 0+ states and 54 invalid / undefined states. Worst case scenario is that KSHELL gives no output. The current best solution is to request exactly 46 0+ states if you want them all.

</p>
</details>

## Notes from before
Mostly outdated info.

<details>
<summary>Click here for notes from before</summary>
<p>

  ### Additions by jorgenem

  I have added some Python scripts in the bin/ folder, namely `shellmodelutilities.py` and `spin_selection.py`. The latter is a small tool to ease setup of calculations, while the first is a comprehensive library of tools to calculate level density (NLD) and gamma-ray strength function (gSF) from shell model files. 

  The folder example_nld_gsf/ contains an example of just that, using the `shellmodelutilities` library. There is also an example summary file on Ne20 with the USDa interaction, to demonstrate the use of the script. The calculated NLD and gSF is not very interesting, however, but I cannot put a large file on Github. If you like, you can download a more interesting calculation summary file from the supplemental material to our PRC on M1 systematics ([arXiv:1807.04036 [nucl-th]](https://arxiv.org/abs/1807.04036)) from this link: https://doi.org/10.5281/zenodo.1493220

  ### Technical notes (NB: THESE CHANGES WERE OVERWRITTEN IN THE VERSION 2 UPDATE OF KSHELL (2021-04-29))
  * I have modified the `transit.f90` file slightly so it prints transition strengths with more decimal precision, to facilitate the gSF calculations. I have updated `collect_logs.py` accordingly. 
  * I have modified `collect_logs.py` to ensure it does not double-count transitions. 
  * I have added some lines to kshell_ui.py so that it does an automatic backup of all the text files from the run into a folder called `KSHELL_runs` under the home path. This is mainly useful when running on a supercomputer, where the calculation is typically run on a scratch disk where files are deleted after some weeks.

</p>
</details>

### Notes to self
MPI compile wrapper mpiifort
intel/2020b og Python/3.8.6-GCCcore-10.2.0
100 lowest states for spins 0 to 14 took 39 minutes on Fram with 32 nodes
