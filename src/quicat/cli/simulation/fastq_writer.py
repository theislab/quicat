from pandas import DataFrame

from .simulator_utils import generate_quality_string


class FastqWriter:
    def __init__(self, output_prefix: str) -> None:
        """
        Initialize the FastqWriter.

        Parameters:
            output_prefix (str): The prefix for the output FASTQ file(s).
        """
        self.output_prefix = output_prefix

    def write_fastq(self, df: DataFrame, is_dna: bool = True) -> None:
        """
        Writes FASTQ file(s) based on the input DataFrame.

        Parameters:
            df (DataFrame): DataFrame containing the sequences to write.
                            For DNA, it should have a 'read' column.
                            For 10x, it should have 'read_1' and 'read_2' columns.
            is_dna (bool): Whether the input data is DNA. If False, it's assumed to be from 10x. Default is True.
        """
        if is_dna:
            self._write_single_fastq(df)
        else:
            self._write_paired_fastq(df)

    def _write_single_fastq(self, df: DataFrame) -> None:
        """
        Writes a single FASTQ file for DNA reads.

        Parameters:
            df (DataFrame): DataFrame containing the 'read' column.
        """
        with open(f"{self.output_prefix}_R1.fastq", "w") as fq:
            for idx, row in df.iterrows():
                read = row["read"]
                quality = generate_quality_string(len(read))
                fq.write(f"@read_{idx + 1}\n{read}\n+\n{quality}\n")

    def _write_paired_fastq(self, df: DataFrame) -> None:
        """
        Writes paired-end FASTQ files for 10x reads.

        Parameters:
            df (DataFrame): DataFrame containing the 'read_1' and 'read_2' columns.
        """
        with open(f"{self.output_prefix}_R1.fastq", "w") as fq1, open(f"{self.output_prefix}_R2.fastq", "w") as fq2:
            for idx, row in df.iterrows():
                read1 = row["read_1"]
                read2 = row["read_2"]
                quality1 = generate_quality_string(len(read1))
                quality2 = generate_quality_string(len(read2))
                fq1.write(f"@read_{idx + 1}\n{read1}\n+\n{quality1}\n")
                fq2.write(f"@read_{idx + 1}\n{read2}\n+\n{quality2}\n")
