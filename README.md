# Robust Personal Sound Zone Control (RACC-PM)

MATLAB research code and representative listening demonstrations for the
manuscript **“A Unified LMI-Based Framework for Robust Personal Sound Zone
Control.”** The work is currently under major revision for IEEE/ACM
Transactions on Audio, Speech, and Language Processing; it must not yet be
cited or described as a published TASLP article.

## Scope

The repository implements a worst-case robust framework that combines
acoustic contrast control and pressure matching under norm-bounded acoustic
transfer-function (ATF) uncertainty. It includes:

- the maintained decomposed RACC-PM semidefinite program;
- the raw/global counterpart used for numerical-equivalence and runtime tests;
- ACC, ACC-Reg, PM, ACC-PM, WCRACC, WCRPM, and POTDC-RACC baselines;
- the NoCT-WCRACC/Full-WCRACC cross-term ablation;
- the SICER-VAST and M-ACC comparison implementations used in the revision;
- AC, NSRE, array-effort, and planarity evaluation functions; and
- scripts and listening examples added during the major revision.

The paper reports a two-zone, single-program validation. The code should not
be interpreted as a jointly optimized arbitrary multi-program controller.

## Repository layout

```text
.
├── src/
│   ├── algorithm/       Core and baseline filter-design methods
│   ├── evaluations/     AC, NSRE, array effort, and planarity
│   ├── simulations/     Array, monitor-grid, and RIR/ATF generation
│   ├── utils/           Shared numerical and data-loading utilities
│   └── visualization/   Sound-field and metric plotting utilities
├── lib/rir_generator/   RIR generator source and Windows MEX binary
├── RQ_response/         Experiments added for the major revision
│   └── listening_demo/  Three 10-algorithm BZ/DZ listening comparisons
├── step1_data_generator.m
└── generate_plane_wave_target.m
```

Large generated MAT files, intermediate solver checkpoints, measured cabin
data, logs, and editable MATLAB figures are intentionally excluded.

## Algorithm names

The internal function/field names retained for compatibility map to the paper
as follows:

| Paper name | Main implementation |
|---|---|
| ACC | `ACC_Unregularized.m` |
| ACC-Reg | `ACC.m` |
| ACC-PM | `ACC_PM.m` |
| WCRACC | `wcACC.m` |
| WCRPM | `RPM.m` |
| NoCT-WCRACC | `NoCT_WCRACC.m` |
| Full-WCRACC | `Full_WCRACC.m` |
| RACC-PM | `RACC_PM_Sub.m` |
| Raw/global RACC-PM | `RACC_PM_GlobalMatched.m` |

`RACC_PM_Sub.m` is the maintained implementation used for the full-scale
RACC-PM experiments. `RACC_PM_GlobalMatched.m` is retained for the matched
global-versus-decomposed benchmark.

## Requirements

The revision experiments were run with MATLAB R2024b. The core code is
expected to work with recent MATLAB releases, subject to the following
dependencies:

- CVX with MOSEK for the SDP/SOCP formulations;
- Signal Processing Toolbox;
- Parallel Computing Toolbox for the large sweeps;
- Statistics and Machine Learning Toolbox for selected analyses; and
- Audio Toolbox plus external PESQ/PEAQ implementations only for the
  perceptual experiment.

The supplied `rir_generator.mexw64` is Windows-specific. Source files are
included in `lib/rir_generator/` for rebuilding on other platforms.

## Reproducing the simulation workflow

Run MATLAB from the repository root.

1. Add the source and RIR generator to the MATLAB path:

   ```matlab
   addpath(genpath('src'));
   addpath(genpath('lib'));
   addpath('src/evaluations', '-begin');
   ```

2. Generate the array geometry and temperature-dependent RIR/ATF database:

   ```matlab
   step1_data_generator
   generate_plane_wave_target
   ```

   The full database is large and can take substantial time and memory.

3. Run the matched full-band benchmark used in the revision:

   ```matlab
   run('RQ_response/run_all_algorithms_fair_evaluation.m')
   ```

The major-revision experiments, their roles, and additional data requirements
are listed in [`RQ_response/README.md`](RQ_response/README.md).

Before running any pipeline, confirm that MATLAB resolves the maintained
evaluator first:

```matlab
which evaluate_performance -all
```

The expected first result is `src/evaluations/evaluate_performance.m`.

## Data availability

- Simulated RIRs/ATFs can be regenerated locally with the supplied geometry
  and RIR generator.
- The measured cabin ATFs are not distributed because they originate from an
  industrial collaboration and are subject to confidentiality restrictions.
  The corresponding evaluation script is included to document the protocol.
- EBU SQAM research audio is not included and was not used in the public demo
  package.

## Listening demonstrations

`RQ_response/listening_demo/` contains 66 stereo WAV files for speech,
synthetic music, and a 1-kHz tone. Each item provides the desired BZ plane-wave
reference, the desired silent DZ reference, and the BZ/DZ outputs of ten
algorithms. The two channels are pressure signals at a horizontal 18-cm
spatial proxy pair in the selected zone; they are not HRTF-rendered binaural
signals and do not constitute a formal listening test. Within each program
item, both zones, both references, and every algorithm share one common
digital safety gain.

See the demo README and manifest for provenance, algorithm mapping, and
licensing information.

### Interactive web player

The repository also includes a browser-based A/B player that preserves the
playback position when switching algorithms:

**[Open the interactive listening demonstrations](https://minikatty.github.io/Robust-PSZ-algorithm/)**

The page source is stored in [`docs/index.html`](docs/index.html). The link
becomes active after GitHub Pages is enabled for the `main` branch and
`/docs` folder. Until then, the page can be previewed locally from the
repository checkout.

## Citation

Until an archival version is available, please cite the manuscript as under
review rather than as a published TASLP paper:

```bibtex
@misc{zhou2026raccpm,
  author = {Lei Zhou and Yaqi Zhu and Chen Huang and Yuewen Wang and
            Liming Shi and Lu Gan and Hongqing Liu},
  title  = {A Unified LMI-Based Framework for Robust Personal Sound Zone Control},
  year   = {2026},
  note   = {Manuscript under major revision}
}
```

## License

The source code and documentation in this repository are released under the
[MIT License](LICENSE). Third-party materials retain their original licenses;
in particular, the LibriSpeech-derived demo files are distributed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) and require
attribution as documented in
[`RQ_response/listening_demo/README.md`](RQ_response/listening_demo/README.md).

## Contact

Lei Zhou — `zhouleicqupt2016 [at] outlook [dot] com`
