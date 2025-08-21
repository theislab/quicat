# quiCAT

`quiCAT` is a versatile tool divided into two main modules—a CLI and an API.

### **CLI**

The CLI allows users to:

- **Extract genetic barcoding data** from raw sequencing files.
- **Generate synthetic datasets** for testing and benchmarking purposes.

---

### CLI Commands

The CLI exposes two primary commands, both of which take a `yaml` configuration file as input. Refer to the specific documentation for details on how to use each command.

#### **1. `extract`**

Extract genetic barcodes from sequencing data. [extract documentation](extract.md)

#### **2. `simulate`**

Generate synthetic data for testing and evaluation purposes. [simulate documentation](simulate.md)

---

### **API**

The API enables seamless integration of genetic barcoding data into analysis workflows. It includes:

- Handy reading functions to import barcoding data into an [AnnData object](https://anndata.readthedocs.io/en/stable/), a widely used format for single-cell data.
- Compatibility with tools like [scanpy](https://scanpy.readthedocs.io/en/stable/) and the broader [scverse](https://scverse.org) ecosystem, facilitating downstream analyses and visualization.
- Some tools and plotting functions to extend `scanpy`
