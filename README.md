# HybridizationWaveforms

The current version only do hybridization of extrapolated waveforms.

To run the code, go to line 713 of Hybridization.py, enter the start time of matching window (the end time of matching window is sat as -5000 currently), and in line 714, enter the directory $data_dir of extrapolated waveform data (must also contain Horizon.h5 file), and the directory $out_dir in where you want to store the outputs. Then run: sbatch Hybrid.sh
