import sys
import kshell_utilities as ksutil
import numpy as np

def main():
    if sys.argv[1] == "long":
        Ne22_v2 = ksutil.loadtxt(path="Ne22_v2/", load_and_save_to_file=False)
        Ne22_v4 = ksutil.loadtxt(path="Ne22_v4/", load_and_save_to_file=False)
        V50_v2 = ksutil.loadtxt(path="V50_v2/", load_and_save_to_file=False)
        V50_v4 = ksutil.loadtxt(path="V50_v4/", load_and_save_to_file=False)
        Ne22_dev = ksutil.loadtxt(path="Ne22_dev/", load_and_save_to_file=False)
        V50_dev = ksutil.loadtxt(path="V50_dev/", load_and_save_to_file=False)
    
    elif sys.argv[1] == "short":
        Ne22_v2 = ksutil.loadtxt(path="Ne22_v2/", load_and_save_to_file=True)
        Ne22_v4 = ksutil.loadtxt(path="Ne22_v4/", load_and_save_to_file=True)
        V50_v2 = ksutil.loadtxt(path="V50_v2/", load_and_save_to_file=True)
        V50_v4 = ksutil.loadtxt(path="V50_v4/", load_and_save_to_file=True)
        Ne22_dev = ksutil.loadtxt(path="Ne22_dev/", load_and_save_to_file=False)
        V50_dev = ksutil.loadtxt(path="V50_dev/", load_and_save_to_file=False)

    elif sys.argv[1] == "supershort":
        Ne22_v2 = ksutil.loadtxt(path="Ne22_v2/", load_and_save_to_file=True)
        Ne22_v4 = ksutil.loadtxt(path="Ne22_v4/", load_and_save_to_file=True)
        V50_v2 = ksutil.loadtxt(path="V50_v2/", load_and_save_to_file=True)
        V50_v4 = ksutil.loadtxt(path="V50_v4/", load_and_save_to_file=True)
        Ne22_dev = ksutil.loadtxt(path="Ne22_dev/", load_and_save_to_file=True)
        V50_dev = ksutil.loadtxt(path="V50_dev/", load_and_save_to_file=True)

    assert np.all(Ne22_v2.levels == Ne22_v4.levels)
    assert np.all(Ne22_v2.transitions_BE1 == Ne22_v4.transitions_BE1)
    assert np.all(Ne22_v2.transitions_BM1 == Ne22_v4.transitions_BM1)
    assert np.all(Ne22_v2.transitions_BE2 == Ne22_v4.transitions_BE2)

    assert np.all(Ne22_v2.levels == Ne22_dev.levels)
    assert np.all(Ne22_v2.transitions_BE1 == Ne22_dev.transitions_BE1)
    assert np.all(Ne22_v2.transitions_BM1 == Ne22_dev.transitions_BM1)
    assert np.all(Ne22_v2.transitions_BE2 == Ne22_dev.transitions_BE2)

    assert np.all(V50_v2.levels == V50_v4.levels)
    assert np.all(V50_v2.transitions_BE1 == V50_v4.transitions_BE1)
    assert np.all(V50_v2.transitions_BM1 == V50_v4.transitions_BM1)
    assert np.all(V50_v2.transitions_BE2 == V50_v4.transitions_BE2)

    assert np.all(V50_v2.levels == V50_dev.levels)
    assert np.all(V50_v2.transitions_BE1 == V50_dev.transitions_BE1)
    assert np.all(V50_v2.transitions_BM1 == V50_dev.transitions_BM1)
    assert np.all(V50_v2.transitions_BE2 == V50_dev.transitions_BE2)

    print("LEVELS")
    print("------")
    print(f"Ne22_v2: {ksutil.get_timing_data(path='Ne22_v2/log_Ne22_usda_m0p.txt')} s")
    print(f"Ne22_v4: {ksutil.get_timing_data(path='Ne22_v4/log_Ne22_usda_m0p.txt')} s")
    print(f"Ne22_dev: {ksutil.get_timing_data(path='Ne22_dev/log_Ne22_usda_m0p.txt')} s")

    print(f"V50_v2: {ksutil.get_timing_data(path='V50_v2/log_V50_gxpf1a_m0p.txt')} s")
    print(f"V50_v4: {ksutil.get_timing_data(path='V50_v4/log_V50_gxpf1a_m0p.txt')} s")
    print(f"V50_dev: {ksutil.get_timing_data(path='V50_dev/log_V50_gxpf1a_m0p.txt')} s")
    print(f"\nV51_dev: {ksutil.get_timing_data(path='V51_dev/log_V51_gxpf1a_m1n.txt')} s")

    print("TRANSIT")
    print("-------")
    print(f"Ne22_v2: {ksutil.get_timing_data(path='Ne22_v2/log_Ne22_usda_tr_m0p_m0p.txt')} s")
    print(f"Ne22_v4: {ksutil.get_timing_data(path='Ne22_v4/log_Ne22_usda_tr_m0p_m0p.txt')} s")
    print(f"Ne22_dev: {ksutil.get_timing_data(path='Ne22_dev/log_Ne22_usda_tr_m0p_m0p.txt')} s")

    print(f"V50_v2: {ksutil.get_timing_data(path='V50_v2/log_V50_gxpf1a_tr_m0p_m0p.txt')} s")
    print(f"V50_v4: {ksutil.get_timing_data(path='V50_v4/log_V50_gxpf1a_tr_m0p_m0p.txt')} s")
    print(f"V50_dev: {ksutil.get_timing_data(path='V50_dev/log_V50_gxpf1a_tr_m0p_m0p.txt')} s")
    print(f"\nV51_dev: {ksutil.get_timing_data(path='V51_dev/log_V51_gxpf1a_tr_m1n_m1n.txt')} s")

if __name__ == "__main__":
    main()