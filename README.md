# Optimal Pseudo-random Noise for DOA Estimation

From this folder in MATLAB, run the following commands to reproduce the paper figures.

| MATLAB command | Experiment |
| --- | --- |
| `main_snapshots` | Fig. 1. RMSE versus the number of snapshots at −2 dB SNR. |
| `main_SNR` | Fig. 2. RMSE versus SNR with 50 snapshots. |
| `main_measured_rmse` | Fig. 4. DOA estimation RMSE on measured acoustic data. |

To reproduce Fig. 4, open `main_measured_rmse.m` in MATLAB and click **Run**. The recording (`ACDataset_NWPU_Aug2026.mat`) is included in this folder. It was collected in the anechoic chamber at Northwestern Polytechnical University in August 2026 and comes from the [NWPU Anechoic Chamber Dataset on GitHub](https://github.com/tian1yao/NWPUAnechoicChamberDataset-Aug2026).

If you find this code useful in your research, please consider citing our paper once it is published.
