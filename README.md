# Robust-PSZ-algorithm
The LMI-based framework implementation of a robust hybrid method for PSZ.

key word: sound field control; robust control; SDP; robust least square

## Usage Notes

### ⚠️ Important Notes
* **Data Exclusion:** Please note that raw data files (`.mat`) and MATLAB figure files (`.fig`) are **not included** in this repository to keep the package lightweight.
* **Dependencies:** The visualization scripts require the **brewermap** package to render colormaps correctly. 
    * Download link: [brewermap (MATLAB Central File Exchange)](https://ww2.mathworks.cn/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps)
* **Functionality:** A brief description of each script's primary purpose is provided in the file tree comments below.

## System Requirements

* **Memory (RAM):** At least **16 GB** of RAM is required. Some modules utilize parallel computing, which can be memory-intensive. Please ensure your system meets this threshold to avoid out-of-memory errors during large-scale simulations.
* **Software Dependencies:** 
    * MATLAB R2021b or later.
    * **Parallel Computing Toolbox:** Required for accelerated data processing and cross-validation scripts.
    * **CVX Toolbox:** With a professional solver (e.g., MOSEK or SDPT3) installed.

## Data Availability

- **Simulated Data:** All simulated Room Impulse Responses (RIRs) and Acoustic Transfer Functions (ATFs) can be reproduced locally by executing the provided scripts in the `step1_data_generator.m` pipeline.
- **Measured Data:** The real-world car cabin measurement dataset is available upon request. Due to commercial confidentiality and partnership agreements, this dataset is **strictly restricted to academic research purposes only**. It may not be redistributed or used for any commercial applications. 
  - To obtain access to the cloud storage link, please contact the author via email.

## Project Structure

```text
Robust-PSZ-algorithm/
│
├── data/                           # Acoustic Datasets (Local storage)
│   ├── arrayGeometry/              # Array configurations (Loudspeakers/Microphones)
│   │   ├── array.mat               # Raw coordinates of transducers
│   │   ├── array_layout.mat        # Pre-defined topological parameters
│   │   └── arrayGeometry.fig       # Visualization of the array setup
│   ├── Cabin_Measurements/         # Real-world ATF dataset (60 sets from car cabin)
│   ├── MonitorGrid/                # High-density grids for sound field visualization
│   └── SimulateRIR/                # Simulated Room Impulse Responses (RIRs)
│       ├── position/               # RIRs with microphone position perturbations
│       ├── snr/                    # RIRs with varying Signal-to-Noise Ratios
│       └── temperature/            # RIRs with temperature-induced sound speed mismatch
│
├── lib/                            # External toolboxes (e.g., CVX, RIR-Generator)
├── logs/                           # Runtime logs and intermediate variables
├── results/                        # Numerical simulation outputs and figures
├── real_measured_results/          # Results validated using car cabin data
│
├── src/                            # Source Code
│   ├── algorithm/                  # Core algorithm implementations
│   │   ├── vast/                   # VAST (Variable Span Trade-off) method based on the public code
│   │   ├── ACC.m                   # Acoustic Contrast Control (Baseline)
│   │   ├── ACC_PM.m                # Hybrid ACC-PM (Non-robust baseline)
│   │   ├── PM.m                    # Pressure Matching (Baseline)
│   │   ├── POTDC_RACC.m            # Robust ACC via POTDC (Iterative SDP)
│   │   ├── RACC_PM.m               # Proposed Robust ACC-PM (Core framework)
│   │   ├── RACC_PM_GLS.m           # Proposed RACC-PM (Global Large Scale / Monolithic)
│   │   ├── RACC_PM_Sub.m           # Proposed RACC-PM (Decomposed / Efficient version)
│   │   ├── RPM.m                   # Robust Pressure Matching (SOCP-based)
│   │   └── wcACC.m                 # Worst-case Robust ACC (Diagonal loading)
│   │
│   ├── evaluations/                # Performance evaluation metrics
│   │   ├── calculate_AC.m          # Compute Acoustic Contrast (AC)
│   │   ├── calculate_AE.m          # Compute Array Effort (AE)
│   │   ├── calculate_NSRE.m        # Compute Normalized Squared Reproduction Error (NSRE)
│   │   ├── calculate_planarity.m   # Compute sound field Planarity
│   │   ├── evaluate_performance.m  # Main performance evaluation wrapper
│   │   └── evaluate_performance_V2.m # Updated performance evaluation wrapper
│   │
│   ├── simulations/                # Acoustic environment and scenario configurations
│   ├── visualization/              # Sound field maps, phase plots, and performance curves
│   |── utils/                      # General utility functions and helper modules
│   |    ├── compute_atf.m          # Compute Acoustic Transfer Functions (ATFs) from RIRs
│   |    ├── configure_freq_parameters.m # Initialize frequency-domain simulation parameters
│   |    ├── design_filters.m       # High-level wrapper for loudspeaker filter design
│   |    ├── evaluate_performance.m # Performance evaluation and metric calculation
│   |    ├── get_bound_paras.m      # Calculate uncertainty bounds for robust optimization
│   |    ├── get_data_filename.m    # Utility for automated data file naming/management
│   |    ├── get_real_measurement_bound.m # Extract uncertainty bounds from measured data (discard)
│   |    ├── log_message.m          # Logging utility for tracking simulation progress
│   |    ├── MaxEigenvector.m       # Math utility: Principal eigenvector extraction
│   |    ├── msal_token_cache_outlook.json # Token cache for the Outlook notification system
│   |    ├── pagemtimes_tmp.m       # Page-wise matrix multiplication (optimization)
│   |    ├── pagenorm_tmp.m         # Page-wise matrix norm calculation
│   |    ├── plane_wave_generator.m # Target sound field (ideal plane wave) generation
│   |    ├── precompute_steering_matrix.m # Precompute steering matrices for planarity metrics
│   |    ├── send_graphmail.m       # MATLAB interface for sending emails via MS Graph API
│   |    ├── send_notification.py   # Python backend for the automated notification system
│   |    └── temp2speed.m           # Convert temperature to sound speed (for robustness analysis)
│   └── debug/                      # Internal debugging scripts and variable validation
│
│   % --- Main Experimental Pipeline ---
├── step1_data_generator.m          # Step 1: Generate/Load ATFs and pre-process data
├── step2_run_cross_validation_V2.m # Step 2: Execute robustness and cross-validation tests
├── step3_results_show.m            # Step 3: Summarize and plot simulation results
├── step4_Pareto_Front.m            # Step 4: Analyze the trade-off between AC and NSRE, sensitivity analysis for weighting parameter
├── step_between_4&5_data_processing.m # Data formatting for real-world validation
├── step5_Real_Measured.m           # Step 5: Validate algorithms using measured cabin data
│
│ % --- Analysis & Utility Scripts ---
├── Revisit_Robustness.m            # In-depth analysis of robustness mechanisms
├── get_design_data.m               # Script to extract experimental configurations
└── README.md                       # Project documentation
```

## Citation

If you find this code or dataset useful for your research, please cite our paper:
> L. Zhou, Y. Zhu, C. Huang, Y. Wang, L. Shi, L. Gan, and H. Liu, "A Unified LMI-Based Framework for Robust Personal Sound Zone Control," *IEEE Transactions on Audio, Speech, and Language Processing*, 2025. (Under Review)

**BibTeX:**
```bibtex
@article{zhou2026robust,
  title={A Unified LMI-Based Framework for Robust Personal Sound Zone Control},
  author={Zhou, Lei and Zhu, Yaqi and Huang, Chen and Wang, Yuewen and Shi, Liming and Gan, Lu and Liu, Hongqing},
  journal={IEEE Transactions on Audio, Speech, and Language Processing},
  year={2026},
  publisher={IEEE}
}
```

## Contact

For any questions, please contact: zhouleicqupt2016@outlook.com.