from turtle import right
import polars as pl
import polars_distance as pld

from umi_tools import network

from cite_seq_count.constants import (
    SUBSET_COLUMN,
    BARCODE_COLUMN,
    FEATURE_NAME_COLUMN,
    R2_COLUMN,
    UMI_COLUMN,
    COUNT_COLUMN,
    UNMAPPED_NAME,
)


def correct_barcodes_pl(
    barcodes_df: pl.DataFrame,
    barcode_subset_df: pl.DataFrame,
    hamming_distance: int,
) -> tuple[pl.LazyFrame, int, dict]:
    """Corrects barcodes using a subset based on join_asof from polars.
    Uses both forward and backward strategy to dinf the closest barcode

    Args:
        barcodes_df (pl.DataFrame): All barcodes with their respective counts
        barcode_subset_df (pl.DataFrame): Barcode reference used to correct
        hamming_distance (int): Max hamming distance allowed
        mapped_barcodes (dict): Dict of mapped barcodes

    Returns:
        tuple[pl.DataFrame, int]: The corrected version of the input barcodes_df, number of corrected barcodes
    """
    print("Correcting barcodes")
    corrected_barcodes_pl = pl.DataFrame(
        schema={
            BARCODE_COLUMN: pl.String,
            "count": pl.UInt32,
            SUBSET_COLUMN: pl.String,
        }
    )
    current_barcodes = pl.DataFrame(
        schema={
            BARCODE_COLUMN: pl.String,
            "count": pl.UInt32,
            SUBSET_COLUMN: pl.String,
        }
    )
    current_barcodes_to_correct = barcodes_df.shape[0]
    last_iteration_barcodes_to_correct = 0
    unknown_barcodes = barcodes_df.filter(
        ~pl.col(BARCODE_COLUMN).is_in(barcode_subset_df[SUBSET_COLUMN])
    )
    n_iterations = 0
    while (
        current_barcodes_to_correct > 0
        and current_barcodes_to_correct != last_iteration_barcodes_to_correct
    ):
        methods = ["backward", "forward"]
        for method in methods:
            current_barcodes = (
                (
                    unknown_barcodes.filter(
                        (
                            ~pl.col(BARCODE_COLUMN).is_in(
                                corrected_barcodes_pl[BARCODE_COLUMN]
                            )
                        )
                    )
                    .sort(BARCODE_COLUMN)
                    .join_asof(
                        barcode_subset_df.sort(SUBSET_COLUMN),
                        left_on=BARCODE_COLUMN,
                        right_on=SUBSET_COLUMN,
                        strategy=method,  # type: ignore
                    )
                )
                .filter(~pl.col(SUBSET_COLUMN).is_null())
                .with_columns(
                    pld.col(BARCODE_COLUMN)
                    .dist_str.hamming(pl.col(SUBSET_COLUMN))
                    .cast(pl.UInt32)
                    .alias("hamming_distance")
                )
                .filter(pl.col("hamming_distance") <= hamming_distance)
                .drop("hamming_distance")
            )
            corrected_barcodes_pl = pl.concat(
                [
                    corrected_barcodes_pl,
                    current_barcodes,
                ]
            )
        current_barcodes_to_correct = current_barcodes.shape[0]
        barcode_subset_df = barcode_subset_df.filter(
            ~pl.col(SUBSET_COLUMN).is_in(corrected_barcodes_pl[SUBSET_COLUMN])
        )
        n_iterations += 1
    print(f"Corrected barcodes in {n_iterations} iterations")
    print(f"Number of uncorrected barcodes: {current_barcodes.shape[0]}")
    mapped_barcodes = dict(
        corrected_barcodes_pl.select(BARCODE_COLUMN, SUBSET_COLUMN).iter_rows()  # type: ignore
    )
    final_corrected = (
        barcodes_df.with_columns(
            pl.col(BARCODE_COLUMN).map_dict(mapped_barcodes, default=pl.first())
        )
        .group_by(BARCODE_COLUMN)
        .agg(pl.sum("count"))
    )
    print("Barcodes corrected")
    n_corrected_barcodes = corrected_barcodes_pl.shape[0]

    return final_corrected.lazy(), n_corrected_barcodes, mapped_barcodes


def summarise_unmapped_df(main_df: pl.LazyFrame, unmapped_r2_df: pl.LazyFrame):
    """Merge main df and unmapped df to get a summary of the unmapped reads

    Args:
        main_df (pl.DataFrame): _description_
        unmapped_r2_df (pl.DataFrame): _description_
    """
    unmapped_r2_df = (
        unmapped_r2_df.filter(pl.col(FEATURE_NAME_COLUMN) == UNMAPPED_NAME)
        .join(main_df, on=R2_COLUMN, how="left")
        .with_columns(
            pl.when(pl.col(FEATURE_NAME_COLUMN).is_null())
            .then(pl.col(R2_COLUMN))
            .otherwise(pl.col(FEATURE_NAME_COLUMN))
            .alias(FEATURE_NAME_COLUMN)
        )
    )
    unmapped_df = (
        unmapped_r2_df.group_by(R2_COLUMN)
        .agg(pl.sum(COUNT_COLUMN))
        .sort(COUNT_COLUMN, descending=True)
        .head(1000)
    )

    return unmapped_df


def generate_umi_counts(read_counts: pl.LazyFrame) -> pl.LazyFrame:
    """Generate umi counts from read counts

    Args:
        read_counts (pl.DataFrame): Read counts

    Returns:
        pl.DataFrame: Umi counts
    """
    umi_counts = read_counts.group_by([BARCODE_COLUMN, FEATURE_NAME_COLUMN]).agg(
        pl.count()
    )
    return umi_counts


def update_main_df(main_df: pl.LazyFrame, mapped_barcodes: dict) -> pl.LazyFrame:
    """Update the main data df with the corrected barcodes

    Args:
        main_df (pl.DataFrame): Data of all reads
        mapped_barcodes (dict): Mapped barcodes from correction

    Returns:
        pl.DataFrame: Data of all reads with barcodes corrected
    """
    main_df = (
        main_df.with_columns(
            pl.col(BARCODE_COLUMN).map_dict(mapped_barcodes, default=pl.first())
        )
        .group_by([BARCODE_COLUMN, UMI_COLUMN, R2_COLUMN])
        .agg(pl.sum(COUNT_COLUMN))
    )
    return main_df


def correct_umis_df(
    main_df: pl.LazyFrame,
    mapped_r2_df: pl.LazyFrame,
    umi_distance: int,
    barcode_subset: pl.LazyFrame,
    cluster_method: str = "directional",
    max_umis: int = 20000,
) -> tuple[pl.LazyFrame, int, list]:
    """Take main df and mapped reads to correct umis.

    Args:
        main_df (pl.LazyFrame): Main df mapping r2, barcode and umi
        mapped_r2_df (pl.LazyFrame): mapped reads
        umi_distance (int): Hamming distance for umi correction
        cluster_method (str, optional): Cluster methods used by the umi clusterer. Defaults to "directional".
        max_umis (int, optional): Threshold for clustered cells. Defaults to 20000.

    Returns:
        tuple[pl.LazyFrame, int, list]: read counts, number of umis corrected, clustered cells lf
    """
    merged_lf = mapped_r2_df.join(main_df, on=R2_COLUMN).join(
        barcode_subset, left_on=BARCODE_COLUMN, right_on=SUBSET_COLUMN, how="inner"
    )
    if umi_distance > 0:
        print("UMI correction")
        temp = (
            merged_lf.with_columns(umi=pl.col(UMI_COLUMN).cast(pl.Binary))
            .drop(R2_COLUMN)
            .group_by([FEATURE_NAME_COLUMN, BARCODE_COLUMN])
            .agg(pl.struct(pl.col(UMI_COLUMN), pl.col(COUNT_COLUMN)))
        ).collect()

        clustered_cells = (
            temp.filter(pl.col(UMI_COLUMN).list.len() > max_umis)
            .select(BARCODE_COLUMN)
            .unique()
            .get_column(BARCODE_COLUMN)
            .to_list()
        )

        umi_mapping_lf = find_umis_to_correct(
            temp=temp,
            cluster_method=cluster_method,
            max_umis=max_umis,
            umi_distance=umi_distance,
        )
        read_counts = update_umis_and_create_read_counts(
            merged_lf=merged_lf, umi_mapping_lf=umi_mapping_lf
        )

        n_corrected_umis = umi_mapping_lf.collect().shape[0]
    else:
        read_counts = (
            merged_lf.group_by([BARCODE_COLUMN, FEATURE_NAME_COLUMN, UMI_COLUMN])
            .agg(pl.sum(COUNT_COLUMN))
            .drop(R2_COLUMN)
        )
        clustered_cells = []
        n_corrected_umis = 0

    return read_counts, n_corrected_umis, clustered_cells


def find_umis_to_correct(
    temp: pl.DataFrame, cluster_method: str, max_umis: int, umi_distance: int
) -> pl.LazyFrame:
    """Iterate through all umis that might need correcting and return a mapping lf

    Args:
        temp (pl.DataFrame): Filtered aggregated UMI per barcode per features
        cluster_method (str): What cluster method to use
        max_umis (int): Threshold for clustered cells
        umi_distance (int): Hamming distance for umi correction

    Returns:
        pl.LazyFrame: Mapping to correct umis
    """
    mapping_list = []
    umi_clusterer = network.UMIClusterer(cluster_method=cluster_method)
    for feature_name, barcode, umis in temp.filter(
        (pl.col(UMI_COLUMN).list.len() > 1) & (pl.col(UMI_COLUMN).list.len() < max_umis)
    ).iter_rows():
        corrected_umis = correct_umis(
            umis,
            umi_distance=umi_distance,
            umi_clusterer=umi_clusterer,
        )
        if len(corrected_umis) == 0:
            continue
        for umi_set in corrected_umis:
            for index, umi in enumerate(umi_set):
                if index != 0:
                    mapping_list.append([feature_name, barcode, umi, umi_set[0]])
    return get_umi_mapping_lf(mapping_list)


def update_umis_and_create_read_counts(
    merged_lf: pl.LazyFrame, umi_mapping_lf: pl.LazyFrame
) -> pl.LazyFrame:
    """Update corrected umis and format to a read count table

    Args:
        merged_lf (pl.LazyFrame): Main df and R2 df joined
        umi_mapping_lf (pl.LazyFrame): UMIs to update

    Returns:
        pl.LazyFrame: Final read counts
    """
    return (
        merged_lf.join(
            umi_mapping_lf,
            on=[FEATURE_NAME_COLUMN, BARCODE_COLUMN, UMI_COLUMN],
            how="left",
        )
        .with_columns(
            pl.when(pl.col("replace").is_null())
            .then(pl.col(UMI_COLUMN))
            .otherwise(pl.col("replace"))
            .alias(UMI_COLUMN)
        )
        .drop("replace")
        .group_by([BARCODE_COLUMN, FEATURE_NAME_COLUMN, UMI_COLUMN])
        .agg(pl.sum(COUNT_COLUMN))
    )


def get_umi_mapping_lf(mapping_list: list) -> pl.LazyFrame:
    """Convert a list of UMI per barcode per features to correct to a lazy frame

    Args:
        mapping_list (list): List of UMIs to correct

    Returns:
        pl.LazyFrame: LazyFrame to correct UMIs
    """
    mapping_df = (
        pl.LazyFrame(
            mapping_list,
            schema={
                FEATURE_NAME_COLUMN: pl.String,
                BARCODE_COLUMN: pl.String,
                UMI_COLUMN: pl.Binary,
                "replace": pl.Binary,
            },
        )
        .with_columns(
            pl.col(UMI_COLUMN).cast(pl.String),
            pl.col("replace").cast(pl.String),
        )
        .drop("orig")
    )
    return mapping_df


def correct_umis(umis_list, umi_distance, umi_clusterer) -> list[list]:
    """Find corrected umis from a pl.struct of UMI and counts

    Args:
        umis_list (_type_): pl.struct of UMI and their counts
        umi_distance (_type_): Hamming distance for a correction
        umi_clusterer (_type_): umi_cluster object

    Returns:
        list[list]: List of list. First member is the
    """
    umis = dict([(i["umi"], i["count"]) for i in umis_list])

    res = umi_clusterer(umis, umi_distance)
    corrected = [corrected_umis for corrected_umis in res if len(corrected_umis) > 1]
    return corrected


# def generate_mtx_counts(
#     main_df: pl.DataFrame,
#     barcode_subset: pl.DataFrame,
#     mapped_r2_df: pl.DataFrame,
#     data_type: str,
# ) -> pl.DataFrame:
#     if data_type == "read":
#         return (
#             main_df.join(barcode_subset, on=BARCODE_COLUMN)
#             .join(mapped_r2_df, on=R2_COLUMN)
#             .group_by([BARCODE_COLUMN, FEATURE_NAME_COLUMN])
#             .agg(pl.sum(COUNT_COLUMN))
#         )
#     else:
#         return (
#             main_df.join(barcode_subset, on=BARCODE_COLUMN)
#             .join(mapped_r2_df, on=R2_COLUMN)
#             .group_by([BARCODE_COLUMN, FEATURE_NAME_COLUMN])
#             .agg(pl.count())
#         )


def find_closest_match(
    df: pl.LazyFrame,
    source_column: str,
    target_df: pl.LazyFrame,
    Levenshtein_distance=1,
):
    return df.join(
        target_df, left_on=source_column, right_on=target_df.columns[0], how="cross"
    ).with_columns(
        pld.col(source_column)
        .dist_str.levenshtein(target_df.columns[0])
        .alias("distance")
    )
