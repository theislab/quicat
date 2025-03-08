import polars as pl
import scanpy as sc
import scipy.sparse as sp

METADATA_DNA = ["sample", "replicate", "condition"]
METADATA_SC = ["cell", "sample", "replicate", "condition"]


def read_dna(csv_path: str) -> sc.AnnData:
    """
    This function reads the bulkDNA CSV file from the extract pipeline into a Polars DataFrame, extracts the specified
    columns as metadata, and creates a sparse matrix from the remaining columns.
    The resulting AnnData object is prepared with the metadata and count data.

    Parameters:
        csv_path (str): Path to the CSV file.
        columns (list of str): List of column names to extract as cell metadata.

    Returns:
        anndata.AnnData: The constructed AnnData object with the sample metadata and count matrix.
    """

    df = pl.read_csv(csv_path)
    metadata_columns_present = [col for col in (METADATA_DNA) if col in df.columns]

    cell_metadata = df.select(metadata_columns_present)
    barcodes = [col for col in df.columns if col not in metadata_columns_present]
    matrix = sp.csr_matrix(df.select(barcodes).to_numpy()).astype(float)

    adata = sc.AnnData(X=matrix, var={"var_names": barcodes}, obs=cell_metadata.to_pandas())
    adata.obs_names = adata.obs["sample"].values

    for col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].astype("category")
    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None, var_type="barcodes")
    return adata


def read_sc(csv_path: str) -> sc.AnnData:
    """
    This function reads the single-cell or spatial transcriptomics CSV file from the extract pipeline into a Polars DataFrame, extracts the specified
    columns as metadata, and creates a sparse matrix from the remaining columns.
    The resulting AnnData object is prepared with the metadata and expression data.

    Parameters:
        csv_path (str): Path to the CSV file.
        columns (list of str): List of column names to extract as cell metadata.

    Returns:
        anndata.AnnData: The constructed AnnData object with the cell/spot metadata and expression matrix.
    """

    df = pl.read_csv(csv_path)
    metadata_columns_present = [col for col in (METADATA_SC) if col in df.columns]

    cell_metadata = df.select(metadata_columns_present)
    barcodes = [col for col in df.columns if col not in metadata_columns_present]
    matrix = sp.csr_matrix(df.select(barcodes).to_numpy()).astype(float)

    adata = sc.AnnData(X=matrix, var={"var_names": barcodes}, obs=cell_metadata.to_pandas())
    adata.obs_names = adata.obs["cell"].values
    adata.obs_names_make_unique()
    del adata.obs["cell"]

    for col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].astype("category")
    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None, var_type="barcodes")
    return adata
