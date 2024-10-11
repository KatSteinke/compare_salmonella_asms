"""A collection of input functions for the workflow"""

from typing import Literal

import pandas as pd

def get_reads_from_overview(sample: str, sample_overview: pd.DataFrame,
                            read_type: Literal["r1", "r2", "extra"]) -> str:
    """Get the path to reads of the specified type for the given sample.

    Arguments:
        sample:             the isolate number for which reads should be retrieved
        sample_overview:    the file of filenames mapping sample number to read files
        read_type:          the read type (r1, r2 or extra - the last for long reads)

    Returns:
        The path to the requested read file for the isolate.

    Raises:
        KeyError:   if the requested isolate is not present in the file or does not have the
                    requested kind of reads (e.g. forward read for ONT data)
        ValueError: if an invalid read type is requested
    """
    valid_types = ['r1', 'r2', 'extra']
    if read_type not in valid_types:
        raise ValueError(f"{read_type} is not a valid read type. "
                         f"Valid read types are {valid_types}.")
    if sample not in sample_overview["sample"].values:
        missing_sample_error = f"Isolate {sample} not found in read overview."
        raise KeyError(missing_sample_error)
    sample_read = sample_overview.loc[sample_overview["sample"] == sample, read_type].squeeze()
    if pd.isna(sample_read):
        missing_read_error = f"No {read_type} reads found for {sample}."
        raise KeyError(missing_read_error)

    return sample_read