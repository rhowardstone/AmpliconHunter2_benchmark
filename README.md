# AmpliconHunter2_benchmark
Reproducibility scripts for the manuscript of AmpliconHunter2

### Setup & Installation

Please run the following to install AmpliconHunter1.1, our test data, and prepare test directories:

```bash
git clone https://github.com/rhowardstone/AmpliconHunter2_benchmark.git
cd AmpliconHunter2_benchmark
bash setup.sh
```
(Note: If you do not have gdown installed, you will be prompted to download the compressed set of 204,800 bacterial genomes manually from: [Link](https://drive.google.com/file/d/1Nt7MjwfL3pIa5Axa3I2z3xoO2Ait_ito/view?usp=drive_link))

### Reproduction of results:
```bash
bash run_benchmarks.sh results
```

The following shorthand was applied in testing:
 - **AHv1** is *AmpliconHunter1.1* in the manuscript,

- **AHv2** is *AmpliconHunter2.\alpha*,

- **AHv3** is *AmpliconHunter2.\beta*, and

- **AHv4** is *AmpliconHunter2*.

