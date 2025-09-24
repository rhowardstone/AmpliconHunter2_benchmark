# [AmpliconHunter2](https://github.com/rhowardstone/AmpliconHunter2)_benchmark
Reproducibility scripts for the manuscript of [AmpliconHunter2](https://github.com/rhowardstone/AmpliconHunter2)

### Setup & Installation

Please run the following to install [AmpliconHunter1.1](https://github.com/rhowardstone/AmpliconHunter), our test data, and prepare test directories:

```bash
git clone https://github.com/rhowardstone/AmpliconHunter2_benchmark.git
cd AmpliconHunter2_benchmark
bash setup.sh
```
(Note: If you do not have gdown installed, you will be prompted to download the compressed set of 204,800 bacterial genomes manually from: [Link](https://drive.google.com/file/d/1Nt7MjwfL3pIa5Axa3I2z3xoO2Ait_ito/view?usp=drive_link))

### Reproduction of results:
```bash
bash run_benchmarks.sh results
python generate_publication_figures.py results results-plots
```

The following shorthand was applied in testing:
 - **AHv1** is [*AmpliconHunter1.1*](https://github.com/rhowardstone/AmpliconHunter) in the manuscript,

- **AHv2** is *AmpliconHunter2.&alpha;*,

- **AHv3** is *AmpliconHunter2.&beta;*,

- **AHv4** is *AmpliconHunter2.&gamma;*, and

- **AHv5** is [*AmpliconHunter2*](https://github.com/rhowardstone/AmpliconHunter2).
