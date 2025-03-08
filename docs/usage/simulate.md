# Quicat Simulate Command Documentation

The `simulate` command generates synthetic genetic barcoding data for DNA and single-cell (sc) applications. The output includes FASTQ files, ground truth CSVs, and optionally, sequence data for downstream analysis.

# **Outputs**

- **Ground Truth CSV**: Contains barcode-cell relationships.
- **FASTQ Files**: `_R1.fastq` and `_R2.fastq` files (if `art` is enabled).
- **FASTA Sequences**: If `art` is disabled, sequences are saved for later processing.

# **Usage**

To run a simulation, use:

```sh
quicat simulate --config <path_to_config.yaml>
```

# **Input Files**

## YAML File Example

This configuration file specifies all parameters for the simulation.

```yaml
general:
  output_path: "" # Path to save the ground truth CSV file
  samples: [] # Prefix for output FASTQ files
  simulation_type: "sc" # Options: 'dna', 'sc'
  n_threads: 1 # Number of threads to use

common:
  art: False # Generate FASTQ directly or return sequences in FASTA
  flank_left: "AGT" # Sequence prepended to each barcode
  flank_right: "TCA" # Sequence appended to each barcode
  barcode_length: 20 # Barcode length
  min_hamming_dist: 4 # Minimum Hamming distance between barcodes
  sequence_length: 100 # Total sequence length including flanks
  num_pcr_chimeras: 5 # Number of PCR chimeras per barcode
  max_pcr_hamming_dist: 2 # Max allowed Hamming distance for PCR chimeras
  pcr_error_rate: 0.1 # Mutation probability during PCR
  pcr_error_fraction: 0.1 # Fraction of reads/UMIs replaced by PCR chimeras

distribution_params:
  distribution_type: "uniform" # Options: 'uniform', 'random', 'normal', 'powerlaw'
  counts: 1 # Used for 'uniform' distribution
  max_count: 100 # Max count for 'random' distribution
  mean: 50 # Mean for 'normal' distribution
  std_dev: 10 # Standard deviation for 'normal' distribution
  exponent: 1.5 # Exponent for 'powerlaw' distribution
  min_frequency: 0.001 # Minimum barcode frequency

dna:
  num_barcodes: 10 # Number of unique barcodes to generate

sc:
  num_cells: 10 # Number of cells to simulate
  num_clones: 2 # Unique clones to generate
  cell_bc_length: 16 # Cell barcode length
  umi_length: 10 # UMI sequence length
  umis_per_cell: 5 # Number of UMIs per cell
```

# **YAML Parameters: In-Depth Explanation**

## **General Parameters**

| Parameter             | Description                              | Default Value |
| --------------------- | ---------------------------------------- | ------------- |
| **`output_path`**     | Directory for storing output files       | `""`          |
| **`samples`**         | Prefix for FASTQ files                   | `[]`          |
| **`simulation_type`** | Simulation type: `dna` or `sc`           | `sc`          |
| **`n_threads`**       | Number of threads for parallel execution | `1`           |

## **Common Simulation Parameters**

| Parameter                        | Description                                           | Default Value   |
| -------------------------------- | ----------------------------------------------------- | --------------- |
| **`art`**                        | Generate FASTQ (`true`) or return sequences (`false`) | `false`         |
| **`flank_left` / `flank_right`** | Flanking sequences added to barcodes                  | `"AGT" / "TCA"` |
| **`barcode_length`**             | Length of the barcode sequences                       | `20`            |
| **`min_hamming_dist`**           | Ensures sufficient barcode diversity                  | `4`             |
| **`sequence_length`**            | Length of extended sequences                          | `100`           |
| **`num_pcr_chimeras`**           | Simulates PCR chimeras                                | `5`             |
| **`pcr_error_rate`**             | Probability of mutations in PCR                       | `0.1`           |

## **Barcode Distribution Parameters**

| Parameter               | Description                                    | Default Value |
| ----------------------- | ---------------------------------------------- | ------------- |
| **`distribution_type`** | Distribution model (`uniform`, `normal`, etc.) | `"uniform"`   |
| **`counts`**            | Count for uniform distribution                 | `1`           |
| **`max_count`**         | Max barcode count (random distribution)        |               |
