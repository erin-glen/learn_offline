# helper functions used in multiple parts of the gp tool
import arcpy
import pandas as pd
from lookups import nlcdParentRollupCategories, nlcdCategories
import numpy as np

arcpy.env.overwriteOutput = True


def feature_class_to_pandas_data_frame(feature_class, field_list):
    """
	Load data into a Pandas Data Frame for subsequent analysis.
	https://joelmccune.com/arcgis-to-pandas-data-frame/
	:param feature_class: Input ArcGIS Feature Class.
	:param field_list: Fields for input.
	:return: Pandas DataFrame object.
	"""
    return pd.DataFrame(
        arcpy.da.FeatureClassToNumPyArray(
            in_table=feature_class,
            field_names=field_list,
            skip_nulls=False,
            null_value=-99999,
        )
    )


def rollupToParentClass(df, columnsToGather, groupBy=None):
    """
	Rollup values to the standard/parent categories using sum
	:param df: dataframe
	:param columnsToGather: list of columns to gather
	:param groupBy: optional list of columns to use in the groupBy in addition to the NLCD parent classes
	:return: dataframe with values summed using the standard/parent categories
	"""
    df2 = df.copy()  # create a copy to avoid manipulating original
    # use the lookup dictionary to add the standard categories for each row
    df2["NLCD1_parentclass"] = df2["NLCD1_class"].map(nlcdParentRollupCategories)
    df2["NLCD2_parentclass"] = df2["NLCD2_class"].map(nlcdParentRollupCategories)
    groupColumns = ["NLCD1_parentclass", "NLCD2_parentclass"]
    if groupBy is not None:
        groupColumns += groupBy
    gb = df2.groupby(groupColumns)[columnsToGather].sum().reset_index()
    return gb


def showPivot(df, values, index_category, column_category):
    """
	Pivot dataframe like an excel pivot table for display using the value sums
	:param df: dataframe
	:param values: column to sum
	:param index_category: the column to use as index (vertical)
	:param column_category: the column to use as columns (horizontal)
	:return: pivoted dataframe
	"""
    cell_count_pivot = pd.pivot_table(df, values=values, index=[index_category], columns=[column_category], aggfunc=np.sum)
    return cell_count_pivot


def tabulateAreaByStratification(stratificationRaster, valueRaster, outputName, pixelSize=30, colNameArea="Hectares"):
    """
	Tabulate the forest age type values for each of the stratification raster land use change classes
	:param stratificationRaster: the stratification raster (clipped to the AOI)
	:param valueRaster: path to the second raster to tabulate area
	:param outputName: name for the column that stores the valueRaster values
	:param pixelSize: pixel size for tabulating error (default is 30 meters)
	:param colNameArea: string for the column name that contains the area (units are hectares)
	:return: pandas dataframe with the area in HA for all of the value raster types for each stratification type
	"""
    # use the tabulate area spatial analyst tool to get the cross tabulated area for combinations
    crossTab = arcpy.sa.TabulateArea( stratificationRaster, "Value", valueRaster, "Value",
                                      "in_memory/crossTab", pixelSize,)
    # convert the results from the tabulation area to a pandas dataframe (all columns)
    crossTab_df = feature_class_to_pandas_data_frame(crossTab, "*")
    # convert wide data into long by gethering all of the forest type columns
    df_melt = pd.melt(crossTab_df, id_vars=["OBJECTID", "VALUE"], value_name="Area")
    # extract just the raster value from the variable column ("VALUE_123" -> 123)
    df_melt["Raster_VALUE"] = df_melt["variable"].apply(
        lambda x: int(x.split("_")[1])
    )  # grab just the number

    # unpack the NLCD codes from the stratification raster
    df_melt["NLCD1_value"] = df_melt["VALUE"].apply(lambda x: int(x / 100))
    df_melt["NLCD2_value"] = df_melt["VALUE"].apply(lambda x: x % 100)

    # Use the NLCD categories dictionary to get the NLCD categories
    df_melt["NLCD1_class"] = df_melt["NLCD1_value"].map(nlcdCategories)
    df_melt["NLCD2_class"] = df_melt["NLCD2_value"].map(nlcdCategories)

    # remove any rows with no area
    df_melt_filter = df_melt[df_melt["Area"] > 0]  # removes any rows without any area
    df_melt_filter[colNameArea] = df_melt_filter.Area / 10000
    df_melt_filter.rename(columns={"Raster_VALUE": outputName, "VALUE": "StratificationValue"}, inplace=True)

    if arcpy.Exists("in_memory/crossTab"):
        arcpy.Delete_management("in_memory/crossTab")

    return df_melt_filter[["StratificationValue", "NLCD1_class", "NLCD2_class", outputName, colNameArea]]


def ZonalSumByStratification(stratificationRaster, valueRaster, columnName, cellSize=30):
    """
    Use zonal stats by table to sum up the values in a raster by each landuse stratifiction class
    :param stratificationRaster: the landuse stratification raster (clipped to AOI)
    :param valueRaster: the raster to summarize
    :param columnName: the column name to store the totals
    :param cellSize: the processing cell size (default is 30 meters)
    :return: pandas dataframe with the stratification classes and the total summed up
    """
    # zonal stats to sum up the average canopy values per stratification class
    zonalStatsSum = arcpy.sa.ZonalStatisticsAsTable(
        stratificationRaster,
        "Value",
        valueRaster,
        "in_memory/zonalout1",
        statistics_type="SUM",
    )

    # convert to pandas dataframe
    zs_df = feature_class_to_pandas_data_frame(zonalStatsSum, ["VALUE", "COUNT", "SUM"])

    # dataframe clean up
    zs_df.columns = ["StratificationValue", "CellCount", columnName]
    zs_df["NLCD1_value"] = zs_df["StratificationValue"].apply(lambda x: int(x / 100))
    zs_df["NLCD2_value"] = zs_df["StratificationValue"].apply(lambda x: x % 100)
    zs_df["NLCD1_class"] = zs_df["NLCD1_value"].map(nlcdCategories)
    zs_df["NLCD2_class"] = zs_df["NLCD2_value"].map(nlcdCategories)

    # area in meters squared is the cell count times cell size (30 x 30 meters). Divide by 10k to get hectares
    zs_df["Hectares"] = (zs_df["CellCount"].astype('int64') * cellSize ** 2) / 10000

    if arcpy.Exists("in_memory/zonalout1"):
        arcpy.Delete_management("in_memory/zonalout1")

    return zs_df[["StratificationValue", "NLCD1_class", "NLCD2_class", "Hectares", "CellCount", columnName]]


def df2jsonstr(df_dict):
    """
    convert a dictionary of pandas dataframe to a json like string
    :param df_dict: dictionary of pandas data frames
    :return: json like string
    """
    return str({k: v.to_json(orient="records") for (k,v) in df_dict.iteritems()})
