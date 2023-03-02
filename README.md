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
