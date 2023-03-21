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

## Usage
  
  ### All usage instructions will be / have been moved to the Wiki on this repository! Please see the Wiki for detailed instructions.
          
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
  
#### MPI issues
  Remember that preprocessor lines cannot be indented. I had a problem with program crashes due to `MPI` errors and this was simply because three preprocessor lines were indented and therefore didnt compile properly.
  
#### WARNING: norm of initital vector too small *** decrease mpi_cnk *** / failed to read wf
  There is a parameter in `constant.f90` called `mpi_cnk` but I have had no luck in removing this error message by reducing the numeric value of this parameter. I have solved the issue by re-calculating the `.wav` file in question, and at one point I had to re-run the calculation of the `.wav` with `n_block = 0`.

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
