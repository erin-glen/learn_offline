# helper functions used in multiple parts of the gp tool
import arcpy
import os
import pandas as pd
from datetime import datetime
arcpy.env.overwriteOutput = True
from lookups import disturbanceLookup
from funcs import tabulateAreaByStratification, ZonalSumByStratification

"""
This is a test
"""

# pandas options
pd.options.mode.chained_assignment = None  # suppressed chained assignment warnings

# check out the spatial analyst extension (required)
if arcpy.CheckExtension("Spatial") == "Available":
    arcpy.CheckOutExtension("Spatial")
else:
    # raise a custom exception
    raise Exception


if __name__ == "__main__":
    wd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print(wd)

    dataFolder = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters"  # CHANGE this to the data folder location
    alternateDataFolder = r"U:\eglen\Projects\LEARN Tools\Data\AlternateData"

    year1 = str(input("Year 1: "))
    year2 = str(input("Year 2: "))
    cellsize = input("Cell size: ")
    aoi_name = input("Name of AOI? ")
    tree_canopy = input("Tree Canopy source? ")

    # hardcoded AOI for development - Montgomery County, Maryland
    # aoi = os.path.join(wd, "data", "AOI", "MontgomeryMD.shp")
    aoi_path = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters\AOI"
    aoi = os.path.join(aoi_path,(aoi_name + ".shp"))


    if "NLCD" in tree_canopy:
        treecanopy_path = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters\TreeCanopy\NLCD"
    else:
        treecanopy_path = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters\TreeCanopy"

    #####ALTER INPUT HERE
    inputConfig = dict(
        nlcd_1=os.path.join(dataFolder, "LandCover", ("NLCD_" + year1 + "_Land_Cover_l48_20210604.tif")),
        nlcd_2=os.path.join(dataFolder, "LandCover", ("NLCD_" + year2 + "_Land_Cover_l48_20210604.tif")),
        forestAgeRaster=os.path.join(dataFolder, "ForestType", "forest_raster_07232020.tif"),
        ####NLCD tree canopy paths > comment out if not in use
        treecanopy_1=os.path.join(treecanopy_path, "nlcd_tcc_conus_" + year1 + "_v2021-4.tif"),
        treecanopy_2=os.path.join(treecanopy_path,  "nlcd_tcc_conus_" + year2 + "_v2021-4.tif"),
        ####HighRes tree canopy paths > comment out if not in use
        #treecanopy_1=os.path.join(treecanopy_path, "nlcd_tcc_conus_" + year1 + "_v2021-4"),
        #treecanopy_2=os.path.join(treecanopy_path, "nlcd_tcc_conus_" + year2 + "_v2021-4"),
        carbon_ag_bg_us=os.path.join(dataFolder, "Carbon", "carbon_ag_bg_us.tif"),
        carbon_sd_dd_lt=os.path.join(dataFolder, "Carbon", "carbon_sd_dd_lt.tif"),
        carbon_so=os.path.join(dataFolder, "Carbon", "carbon_so.tif"),
        disturbanceRasters=[
            #os.path.join(dataFolder, "Disturbances", "disturbance_0104.tif"),
            #os.path.join(dataFolder, "Disturbances", "disturbance_0406.tif"),
            #os.path.join(dataFolder, "Disturbances", "disturbance_0608.tif"),
            #os.path.join(dataFolder, "Disturbances", "disturbance_0811.tif"),
            #os.path.join(dataFolder, "Disturbances", "disturbance_1113.tif"),
            #os.path.join(dataFolder, "Disturbances", "disturbance_1316.tif"),
            os.path.join(dataFolder, "Disturbances", "disturbance_1619.tif")
        ]
    )

    startTime = datetime.now()
    inputConfig["aoi"] = aoi # add the AOI to the inputConfig dictionary

    #define the output directory
    parentOutputDirectory = "U:/eglen/Projects/LEARN Tools/Data/Outputs/"
    dateTime = datetime.now()
    dateFormat = dateTime.strftime("%m_%d")
    #outputFolderName = input("Name of output directory: ")
    outputFolderName = dateFormat + "_" + year1 + "_" + year2 + "_" + aoi_name
    outputPath = os.path.join(parentOutputDirectory, outputFolderName)
    outputDirectory = os.mkdir(outputPath)
    text_doc = os.path.join(outputPath,"doc")
    with open(text_doc,'w') as doc:
        doc.write("Year 1: " + year1 + "\n")
        doc.write("Year 2: " + year2 + "\n")
        doc.write("Cellsize: " + cellsize + "\n")
        doc.write("Date: " + str(datetime.now()) + "\n" + "\n")
        doc.write(str(inputConfig.values()))

def landuseStratificationRaster(nlcdRaster1, nlcdRaster2, aoi):
    """
    Create the two year NLCD stratification raster by offsetting one year by 100 and then adding the second year.
    This allow all landuse change to be tracked between two years. Rasters must follow the NLCD classification schema.
    Example: a pixel has a value of 43 (mixed forest) in the first raster and a value 81 (Pasture/hay) in the second
    representing the transition from mixed forest to pasture, the stratification raster will have a pixel value of 4381.
    :param nlcdRaster1: the path to the NLCD raster to use for the initial state
    :param nlcdRaster2: the path the second NLCD raster to use as the after stat
    :param aoi: the area of interest
    :return: the landuse change raster where the pixel values combines both year classifcation values
    """
    # https://gis.stackexchange.com/questions/94469/is-clip-data-management-or-extract-by-mask-spatial-analyst-more-efficient
    arcpy.Clip_management(nlcdRaster1, "#", "in_memory/nlcd_before", aoi, "", "ClippingGeometry")
    arcpy.Clip_management(nlcdRaster2, "#", "in_memory/nlcd_after", aoi, "", "ClippingGeometry")
    stratRast = arcpy.Raster("in_memory/nlcd_before") * 100 + arcpy.Raster("in_memory/nlcd_after")
    return stratRast


def main(aoi, nlcd_1, nlcd_2, forestAgeRaster, treecanopy_1, treecanopy_2, carbon_ag_bg_us, carbon_sd_dd_lt,
         carbon_so, disturbanceRasters):
    """
    Landuse change stratification summaries for forest age, treecanopy, carbon & disturbance
    :param aoi: area of interest
    :param nlcd_1: path to NLCD landuse for year 1
    :param nlcd_2: path to NLCD landuse for year 2
    :param forestAgeRaster: path to the forest age region type raster
    :param treecanopy_1: path to the tree canopy raster for year 1
    :param treecanopy_2: path to the tree canopy raster for year 2
    :param carbon_ag_bg_us: path to the above ground, below ground, and understory carbon raster
    :param carbon_sd_dd_lt: path to the standing dead, dead down, and litter raster
    :param carbon_so: path to the soil organic carbon raster
    :param disturbanceRasters: list of all the disturbance rasters that encompass the years selected
    :return: pandas data frames for tree canopy, forestAge, disturbance, carbon
    """
    # print out the paths to all the inputs
    arcpy.AddMessage(
        "INPUTS: {} {} {} {} {} {} {} {} {} {}".format(aoi, nlcd_1, nlcd_2, forestAgeRaster, treecanopy_1, treecanopy_2,
                                                       carbon_ag_bg_us, carbon_sd_dd_lt, carbon_so, disturbanceRasters))

    # arcpy environment settings
    arcpy.env.snapRaster = nlcd_1
    arcpy.env.cellSize = cellsize
    arcpy.env.overwriteOutput = True
    arcpy.env.extent = aoi

    # stratification raster
    arcpy.AddMessage("STEP 1: Creating land use stratification raster for all classes of landuse")
    stratRast = landuseStratificationRaster(nlcd_1, nlcd_2, aoi)

    # Tree Canopy - sum up all the pixels values for the average and difference (year 1 - year 2)
    arcpy.AddMessage("STEP 2: Summing up the tree canopy average & difference by stratification class")

    # mask the tree canopy to the area of interest
    # mask the tree canopy to the area of interest
    if treecanopy_1 is not None and treecanopy_2 is not None:

        # calculate the average tree canopy per pixel
        with arcpy.EnvManager(mask=aoi, cellSize=(int(cellsize))):
            treeCanopyAvg = (
                                    arcpy.Raster(treecanopy_1) + arcpy.Raster(treecanopy_2)
                            ) / 2  # average

        # tree canopy loss this only happens on the pixels where there is loss
        # logic would be something like Con(tree2 < tree1, tree1-tree2, 0)
        # find the different where there is tree canopy loss
        with arcpy.EnvManager(mask=aoi, cellSize=(cellsize)):
            treeCanopyDiff = arcpy.sa.Con(
                arcpy.Raster(treecanopy_2) < arcpy.Raster(treecanopy_1),
                arcpy.Raster(treecanopy_1) - arcpy.Raster(treecanopy_2),
                0,
            )

        # Optional hardcoded tree canopy loss: if using this, comment out the above tree canopy loss calculation
        # treeCanopyDiff = r"U:\eglen\Projects\LEARN Tools\Data\AlternateData\Cincinnati\TreeCanopy_Loss.tif"

        # zonal sum for each stratification class
        tc_avg_df = ZonalSumByStratification(
            stratRast, treeCanopyAvg, "TreeCanopy_HA", cellSize=30
        )
        tc_diff_df = ZonalSumByStratification(
            stratRast, treeCanopyDiff, "TreeCanopyLoss_HA", cellSize=30
        )

        # drop columns to avoid duplicate cols
        tc_diff_df.drop(["Hectares", "CellCount"], axis=1, inplace=True)
        # merge the two dataframes
        treeCover = tc_avg_df.merge(
            tc_diff_df,
            how="outer",
            on=["StratificationValue", "NLCD1_class", "NLCD2_class"],
        )

        """
            tree canopy notes...
            cell size is now set to automatically update for TC based on tree canopy selection (1,30)
            for non aggregated 1 m data: 0.0001 
            if running at 1 meter, take * 900 out of formula, so just * 0.0001)
            for aggregated data: 0.0001 (this one you can run at 30 meter cell size)
        """

        if "NLCD" in tree_canopy:
            # ie 50% is 50/100 * 900 square meters  * 0.0001 HA = 0.045 canopy in HA
            treeCover["TreeCanopy_HA"] = (
                    treeCover["TreeCanopy_HA"] / 100 * 900 * 0.0001
            )  # % canopy * pixel size * sq m to ha
            treeCover["TreeCanopyLoss_HA"] = (
                    treeCover["TreeCanopyLoss_HA"] / 100 * 900 * 0.0001
            )  # % canopy * pixel size * sq m to ha
        else:
            treeCover["TreeCanopy_HA"] = (
                    treeCover["TreeCanopy_HA"] * 0.0001 * 900
            )  # % canopy * pixel size * sq m to ha
            treeCover["TreeCanopyLoss_HA"] = (
                    treeCover["TreeCanopyLoss_HA"] * 0.0001 * 900
            )  # % canopy * pixel size * sq m to ha

        # drop columns to avoid duplicates names when merging"StratificationValue", "NLCD1_class", "NLCD2_class"
        treeCover.drop(["Hectares", "CellCount"], axis=1, inplace=True)

    else:
        arcpy.AddMessage("Skipping Tree Canopy - no data.")
        treeCover = pd.DataFrame(columns=["StratificationValue", "NLCD1_class", "NLCD2_class",
                                          "TreeCanopy_HA", "TreeCanopyLoss_HA"])

    # Disturbance - tabulate the area
    arcpy.AddMessage("STEP 3: Cross tabulating disturbance area by stratification class")
    arcpy.AddMessage("Number of disturbance rasters: {}".format(len(disturbanceRasters)))
    if len(disturbanceRasters) == 1:
        disturbRast = disturbanceRasters[0]
    else:
        disturbRast = arcpy.sa.CellStatistics(disturbanceRasters, "MAXIMUM", ignore_nodata="DATA")
    disturbance = tabulateAreaByStratification(stratRast, disturbRast, outputName="Disturbance")
    disturbance["DisturbanceClass"] = disturbance["Disturbance"].map(disturbanceLookup)

    # convert disturbances to wide dataframe and sum classes that are the same together to get total area
    disturbance_wide = disturbance.pivot_table(index=['StratificationValue', 'NLCD1_class', 'NLCD2_class'],
                                               columns='DisturbanceClass', values='Hectares',
                                               aggfunc='sum').reset_index()

    if not any(x in ["fire_HA", "harvest_HA", "insect_damage_HA"] for x in disturbance_wide.columns):
        disturbance_wide = pd.DataFrame(columns=['StratificationValue', 'NLCD1_class', 'NLCD2_class'])

    # Carbon - zonal sum
    arcpy.AddMessage("STEP 4: Zonal statistics sum for carbon rasters by stratification class")
    carbon_ag_bg_us_df = ZonalSumByStratification(stratRast, carbon_ag_bg_us, "carbon_ag_bg_us")
    carbon_sd_dd_lt_df = ZonalSumByStratification(stratRast, carbon_sd_dd_lt, "carbon_sd_dd_lt")
    carbon_so_df = ZonalSumByStratification(stratRast, carbon_so, "carbon_so")

    # drop columns to avoid duplicates names (hectares + areas are the same since the zones are identical)
    carbon_sd_dd_lt_df.drop(["Hectares", "CellCount"], axis=1, inplace=True)
    carbon_so_df.drop(["Hectares", "CellCount"], axis=1, inplace=True)

    # merge the three dataframes into one by joining on the two NLCD classes
    carbon = carbon_ag_bg_us_df.merge(
        carbon_sd_dd_lt_df, how='outer', on=["StratificationValue", "NLCD1_class", "NLCD2_class"]
    ).merge(carbon_so_df, how='outer', on=["StratificationValue", "NLCD1_class", "NLCD2_class"])

    # unit conversions
    # pixels are in metric tons Carbon per hectare
    # metric tons C / hectare * 1 hectare / 10000 square meters * 900 square meters  = metric tons carbon
    carbon["carbon_ag_bg_us"] = carbon["carbon_ag_bg_us"] / 10000 * 900
    carbon["carbon_sd_dd_lt"] = carbon["carbon_sd_dd_lt"] / 10000 * 900
    carbon["carbon_so"] = carbon["carbon_so"] / 10000 * 900

    # merge disturbance area, tree cover totals, carbon totals by stratification class
    groupByLanduseChangeDF = carbon.merge(treeCover, how='outer',
                                          on=["StratificationValue", "NLCD1_class", "NLCD2_class"]).merge(
        disturbance_wide, how='outer', on=["StratificationValue", "NLCD1_class", "NLCD2_class"]
    )

    # Forest Age Type - tabulate the area
    arcpy.AddMessage("STEP 5: Tabulating total area for the forest age types by stratification class")
    forestAge = tabulateAreaByStratification(stratRast, forestAgeRaster, outputName="ForestAgeTypeRegion")

    # merge forestAge area total + disturbance areas by forestAge
    arcpy.AddMessage("STEP 6: Tabulating disturbance area for the forest age types by stratification class for fires")

    # calculate disturbances (insects, fires, harvest) by forestAgeRegionType and landuse change type
    # pull the disturbance codes from the dictionary and loop through them
    disturbance_cats = list(set(disturbanceLookup.values()))
    for d in disturbance_cats:
        pixValues = [k for k, v in disturbanceLookup.items() if v == d]
        tempDisturbanceDF = tabulateAreaByStratification(
            arcpy.sa.Con(arcpy.sa.InList(disturbRast, pixValues), stratRast, ""),
            forestAgeRaster,
            outputName="ForestAgeTypeRegion",
            colNameArea=d)
        forestAge = forestAge.merge(tempDisturbanceDF, how='outer', on=["StratificationValue", "NLCD1_class",
                                                                        "NLCD2_class", "ForestAgeTypeRegion"])

        #fill NA with zero to avoid calculation errors
        forestAge['fire_HA'] = forestAge['fire_HA'].fillna(0)
        forestAge['insect_damage_HA'] = forestAge['insect_damage_HA'].fillna(0)
        forestAge['harvest_HA'] = forestAge['harvest_HA'].fillna(0)

        forestAge['undisturbed_HA'] = forestAge['Hectares'] - forestAge['fire_HA'] - forestAge['harvest_HA'] - \
                                      forestAge['insect_damage_HA']

        # code to return all attribites from forest table
        forest_lookup_csv = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters\ForestType\forest_raster_09172020.csv"
        forest_table = pd.read_csv(forest_lookup_csv)
        #forestAge = forestAge.merge(forest_table, left_on="ForestAgeTypeRegion", right_on="Value")

        # code to return only factors from forest table (and convert to CO2)
        col_list = ['ForestAgeTypeRegion', 'Nonforest to Forest Removal Factor',
                    'Forests Remaining Forest Removal Factor', 'Fire Emissions Factor',
                    'Insect Emissions Factor', 'Harvest Emissions Factor']
        forest_table = pd.read_csv(forest_lookup_csv, usecols=col_list)
        print(forest_table)
        print(forestAge)
        #forestAge = pd.merge(forestAge, forest_table)

        forestAge['Annual_Removals_Undisturbed_C02'] = (
                (forestAge['undisturbed_HA'] * forestAge['Forests Remaining Forest Removal Factor']) * (44 / 12))
        # todo add a field for nonforest to forest hectares
        forestAge['Annual_Removals_N_to_F_C02'] = (
                (forestAge['undisturbed_HA'] * forestAge['Nonforest to Forest Removal Factor']) * (44 / 12))
        forestAge['Annual_Emissions_Fire_CO2'] = (
                (forestAge['fire_HA'] * forestAge['Fire Emissions Factor']) * (44 / 12)
                / (int(year2) - int(year1)))
        forestAge['Annual_Emissions_Harvest_CO2'] = (
                (forestAge['harvest_HA'] * forestAge['Harvest Emissions Factor']) * (44 / 12)
                / (int(year2) - int(year1)))
        forestAge['Annual_Emissions_Insect_CO2'] = (
                (forestAge['insect_damage_HA'] * forestAge['Insect Emissions Factor']) * (44 / 12)
                / (int(year2) - int(year1)))

        # todo add function to calculate emissions using dictionary
        # todo add calculations to dataframe

        # return groupByLanduseChangeDF.sort_values(by=["Hectares"], ascending=False), forestAge.sort_values(
        #     by=["Hectares"],
        #     ascending=False)


    return groupByLanduseChangeDF.sort_values(by=["Hectares"], ascending=False), forestAge.sort_values(by=["Hectares"],
                                                                                                       ascending=False)

# stratification areas/sums for each of the value rasters
landuse_result, forestType_result = main(**inputConfig)

# build the json like string that contains all the results
print(
    '{{"stratificationByLanduse": {}, "stratificationByForestAgeRegionType": {}}}'.format(
        landuse_result.to_json(orient="records"),
        forestType_result.to_json(orient="records")
    )
)

# save the dataframes as csvs
strat_csv = os.path.join(outputPath, "stratificationByLanduse.csv")
strat_forest_csv = os.path.join(outputPath, "stratificationByForestAgeRegionType.csv")

print(strat_csv)
print(strat_forest_csv)

landuse_result.to_csv(strat_csv, index=False)
forestType_result.to_csv(strat_forest_csv, index=False)

print("Total processing time: {}".format(datetime.now() - startTime))
