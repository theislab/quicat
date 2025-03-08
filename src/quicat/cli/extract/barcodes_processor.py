import gc
from typing import Optional

import polars as pl

from .sequences_collapser import SequencesCollapser


class BarcodesProcessor:
    """
    Processes barcode data for sequencing technologies.

    Attributes:
        sequencing_technology (str): The sequencing technology used ('DNA', 'scRNA-seq', etc.).
        filter_barcodes_raw_numbers (int): The raw number threshold for filtering low abundance barcodes.
        filter_barcodes_relative_abundance (float): The relative abundance threshold for filtering barcodes.
        distance_threshold (int): The distance threshold for collapsing similar barcodes.
        n_threads (int): Number of threads to use for parallel processing.
        barcode_ratio (float): Ratio used to collapse similar barcodes.
    """

    def __init__(
        self,
        sequencing_technology: str,
        filter_barcodes_raw_numbers: int,
        filter_barcodes_relative_abundance: float,
        distance_threshold: int,
        n_threads: int,
        barcode_ratio: int,
    ) -> None:
        """
        Initializes BarcodesProcessor with the provided parameters.

        Parameters:
            sequencing_technology (str): The sequencing technology used ('DNA', 'scRNA-seq', etc.).
            filter_barcodes_raw_numbers (int): The raw number threshold for filtering low abundance barcodes.
            filter_barcodes_relative_abundance (float): The relative abundance threshold for filtering barcodes.
            distance_threshold (int): The distance threshold for collapsing similar barcodes.
            n_threads (int): Number of threads to use for parallel processing.
            barcode_ratio (int): Ratio used to collapse similar barcodes.
        """
        self.sequencing_technology = sequencing_technology
        self.filter_barcodes_raw_numbers = filter_barcodes_raw_numbers
        self.filter_barcodes_relative_abundance = filter_barcodes_relative_abundance
        # TODO: single cell only retain top barcode
        self.distance_threshold = distance_threshold
        self.n_threads = n_threads
        self.barcode_ratio = barcode_ratio

    def filter_low_abundant_barcodes(self, data: pl.DataFrame, groupby: Optional[str] = None) -> pl.DataFrame:
        """
        Filters low abundance barcodes from the DataFrame.

        Parameters:
            data (pl.DataFrame): The input DataFrame containing barcode information.

        Returns:
            pl.DataFrame: The filtered DataFrame.
        """

        if self.filter_barcodes_raw_numbers:
            data = data.filter(pl.col("counts") >= self.filter_barcodes_raw_numbers)

        elif self.filter_barcodes_relative_abundance:
            if groupby == "sample":
                data = (
                    data.with_columns(
                        (pl.col("counts") / pl.col("counts").sum().over(groupby) * 100).alias("relative_abundance")
                    )
                    .filter(pl.col("relative_abundance") >= self.filter_barcodes_relative_abundance)
                    .drop("relative_abundance")
                )
            elif groupby == "replicate":
                data = data.with_columns(
                    pl.col("counts").sum().over(["replicate", "barcode"]).alias("total_counts_summed_over_replicate")
                )
                data = (
                    data.with_columns(
                        (
                            pl.col("total_counts_summed_over_replicate") / pl.col("counts").sum().over(groupby) * 100
                        ).alias("relative_abundance")
                    )
                    .filter(pl.col("relative_abundance") >= self.filter_barcodes_relative_abundance)
                    .drop(["total_counts_summed_over_replicate", "relative_abundance"])
                )
        return data

    def _build_dataframe_DNA(self, df: pl.DataFrame) -> pl.DataFrame:
        """
        Builds a DataFrame for DNA sequencing barcodes, filtering and collapsing similar barcodes.

        Parameters:
            data (List[Dict[str, Any]]): The input data containing barcode information.

        Returns:
            pl.DataFrame: A pivoted DataFrame of filtered and collapsed barcode counts.
        """
        df = df[[s.name for s in df if s.null_count() != df.height]]
        grouping_columns, groupby = ["sample"], "sample"

        if "replicate" in df.columns:
            grouping_columns.append("replicate")
            groupby = "replicate"
        if "condition" in df.columns:
            grouping_columns.append("condition")
        if self.filter_barcodes_raw_numbers or self.filter_barcodes_relative_abundance:
            df = self.filter_low_abundant_barcodes(df, groupby=groupby)
        df = df.cast({pl.Int64: pl.Int32})

        # Collapse barcodes
        if self.distance_threshold and self.distance_threshold > 0:
            df = df.sort("counts", descending=True)
            unique_barcodes = df.unique(subset=["barcode"], maintain_order=True)["barcode"].to_list()
            collapser = SequencesCollapser(
                unique_barcodes, error_threshold=self.distance_threshold, n_threads=self.n_threads
            )
            barcode_to_component = collapser.cluster()

            df = df.with_columns(
                pl.Series(name="cluster", values=[barcode_to_component[barcode] for barcode in df["barcode"]])
            )
            if not self.barcode_ratio:
                self.barcode_ratio = 1
            if self.barcode_ratio < 1:
                self.barcode_ratio = 1

            max_counts = df.group_by("cluster").agg(
                [
                    pl.col("counts").max().alias("max_count"),
                    pl.col("barcode").filter(pl.col("counts") == pl.col("counts").max()).first().alias("max_barcode"),
                ]
            )
            df = df.join(max_counts, on="cluster")
            df = df.with_columns(
                pl.when(pl.col("max_count") >= (pl.col("counts") * self.barcode_ratio))
                .then(pl.col("max_barcode"))
                .otherwise(pl.col("barcode"))
                .alias("barcode")
            )
            df = df.group_by([*grouping_columns, "barcode"]).agg(pl.col("counts").sum().alias("counts"))

        df_pivoted = (
            df.pivot(values="counts", index=grouping_columns, on="barcode", aggregate_function="sum")
            .fill_null(0)
            .sort("sample")
        )
        gc.collect()
        return df_pivoted

    def _build_dataframe_sc(self, df: pl.DataFrame) -> pl.DataFrame:
        """
        Builds a DataFrame for single-cell or spatial sequencing barcodes, filtering and collapsing similar barcodes.

        Parameters:
            data (List[Dict[str, Any]]): The input data containing barcode information.

        Returns:
            pl.DataFrame: A pivoted DataFrame of filtered and collapsed barcode counts.
        """
        # df = df.explode(["barcode", "counts"])

        df = df[[s.name for s in df if s.null_count() != df.height]]
        grouping_columns = ["sample", "cell"]
        if "replicate" in df.columns:
            grouping_columns.append("replicate")
        if "condition" in df.columns:
            grouping_columns.append("condition")

        if self.filter_barcodes_raw_numbers or self.filter_barcodes_relative_abundance:
            df = self.filter_low_abundant_barcodes(df, groupby="sample")

        df = df.cast({pl.Int64: pl.Int32})

        # Collapse barcodes
        if self.distance_threshold and self.distance_threshold > 0:
            df = df.sort("counts", descending=True)
            unique_barcodes = df.unique(subset=["barcode"], maintain_order=True)["barcode"].to_list()

            collapser = SequencesCollapser(
                unique_barcodes, error_threshold=self.distance_threshold, n_threads=self.n_threads
            )
            barcode_to_component = collapser.cluster()

            df = df.with_columns(
                pl.Series(name="cluster", values=[barcode_to_component[barcode] for barcode in df["barcode"]])
            )
            if not self.barcode_ratio:
                self.barcode_ratio = 1
            if self.barcode_ratio < 1:
                self.barcode_ratio = 1
            max_counts = df.group_by("cluster").agg(
                [
                    pl.col("counts").max().alias("max_count"),
                    pl.col("barcode").filter(pl.col("counts") == pl.col("counts").max()).first().alias("max_barcode"),
                ]
            )
            df = df.join(max_counts, on="cluster")
            df = df.with_columns(
                pl.when(pl.col("max_count") >= (pl.col("counts") * self.barcode_ratio))
                .then(pl.col("max_barcode"))
                .otherwise(pl.col("barcode"))
                .alias("barcode")
            )
            df = (
                df.group_by(*grouping_columns, "barcode")
                .agg(pl.col("counts").sum().alias("counts"))
                .sort(grouping_columns, descending=[False] * len(grouping_columns))
            )
        df = df.with_columns(
            (pl.col("counts") / pl.col("counts").sum().over([grouping_columns]) * 100).alias("relative_abundance")
        )

        df_pivoted = df.pivot(
            values="counts", index=grouping_columns, on="barcode", aggregate_function="sum"
        ).fill_null(0)
        gc.collect()
        return df_pivoted

    def build_dataframe(self, data: pl.DataFrame) -> pl.DataFrame:
        """
        Builds a DataFrame based on the sequencing technology.

        Parameters:
            data (List[Dict[str, Any]]): The input data containing barcode information.

        Returns:
            pl.DataFrame: A pivoted DataFrame of filtered and collapsed barcode counts.
        """
        if self.sequencing_technology == "DNA":
            return self._build_dataframe_DNA(data)
        else:
            return self._build_dataframe_sc(data)
