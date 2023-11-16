# helper functions used in multiple parts of the gp tool
import arcpy
import os
import pandas as pd
from datetime import datetime

arcpy.env.overwriteOutput = True
import pprint
from lookups import disturbanceLookup, nlcdParentRollupCategories
from funcs import tabulateAreaByStratification, ZonalSumByStratification, calculate_category, save_results, \
    calculate_FRF, fillNA, calculate_plantable, calculate_treeCanopy, calculate_disturbances, zonal_sum_carbon

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
    aoi = os.path.join(aoi_path, (aoi_name + ".shp"))

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
        treecanopy_2=os.path.join(treecanopy_path, "nlcd_tcc_conus_" + year2 + "_v2021-4.tif"),
        plantableAreas="None",
        ####HighRes tree canopy paths > comment out if not in use
        # treecanopy_1=os.path.join(treecanopy_path, "nlcd_tcc_conus_" + year1 + "_v2021-4"),
        # treecanopy_2=os.path.join(treecanopy_path, "nlcd_tcc_conus_" + year2 + "_v2021-4"),
        carbon_ag_bg_us=os.path.join(dataFolder, "Carbon", "carbon_ag_bg_us.tif"),
        carbon_sd_dd_lt=os.path.join(dataFolder, "Carbon", "carbon_sd_dd_lt.tif"),
        carbon_so=os.path.join(dataFolder, "Carbon", "carbon_so.tif"),
        disturbanceRasters=[
            # os.path.join(dataFolder, "Disturbances", "disturbance_0104.tif"),
            # os.path.join(dataFolder, "Disturbances", "disturbance_0406.tif"),
            # os.path.join(dataFolder, "Disturbances", "disturbance_0608.tif"),
            # os.path.join(dataFolder, "Disturbances", "disturbance_0811.tif"),
            # os.path.join(dataFolder, "Disturbances", "disturbance_1113.tif"),
            # os.path.join(dataFolder, "Disturbances", "disturbance_1316.tif"),
            os.path.join(dataFolder, "Disturbances", "disturbance_1619.tif")
        ]
    )

    startTime = datetime.now()
    inputConfig["aoi"] = aoi  # add the AOI to the inputConfig dictionary

    # define the output directory
    parentOutputDirectory = "U:/eglen/Projects/LEARN Tools/Data/Outputs/"
    dateTime = datetime.now()
    dateFormat = dateTime.strftime("%m_%d")
    # outputFolderName = input("Name of output directory: ")
    outputFolderName = dateFormat + "_" + year1 + "_" + year2 + "_" + aoi_name
    outputPath = os.path.join(parentOutputDirectory, outputFolderName)

    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

    text_doc = os.path.join(outputPath, "doc")
    with open(text_doc, 'w') as doc:
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
         carbon_so, disturbanceRasters, plantableAreas=None):
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
    treeCover = calculate_treeCanopy(treecanopy_1, treecanopy_2, stratRast, tree_canopy, aoi, cellsize)

    # Plantable areas - sum up pixel values for plantable areas in year 2
    arcpy.AddMessage("STEP 2.5: Summing plantable areas by stratification class")
    calculate_plantable(plantableAreas,stratRast,treeCover,aoi,cellsize)

    ###### --------------- Disturbance - tabulate the area -----------------
    arcpy.AddMessage("STEP 3: Cross tabulating disturbance area by stratification class")
    arcpy.AddMessage("Number of disturbance rasters: {}".format(len(disturbanceRasters)))
    disturbance_wide, disturbRast = calculate_disturbances(disturbanceRasters, stratRast)


    # Carbon - zonal sum
    arcpy.AddMessage("STEP 4: Zonal statistics sum for carbon rasters by stratification class")
    carbon = zonal_sum_carbon(stratRast, carbon_ag_bg_us, carbon_sd_dd_lt, carbon_so)


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

    # map the NLCD class to the parent class in forest age dataframe
    forestAge["NLCD_1_ParentClass"] = forestAge["NLCD1_class"].map(nlcdParentRollupCategories)
    forestAge["NLCD_2_ParentClass"] = forestAge["NLCD2_class"].map(nlcdParentRollupCategories)

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
    # fill empty cells with zero for calaculations
    fillNA(forestAge)

    # todo add to input config
    # enter the path to the forest lookup table
    forest_lookup_csv = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters\ForestType\forest_raster_09172020.csv"

    # create list of columns to read
    col_list = ['ForestAgeTypeRegion', 'Nonforest to Forest Removal Factor',
                'Forests Remaining Forest Removal Factor', 'Fire Emissions Factor',
                'Insect Emissions Factor', 'Harvest Emissions Factor']
    # read the forest able using column list
    forest_table = pd.read_csv(forest_lookup_csv, usecols=col_list)
    # merge forestAge and forest lookup table
    forestAge = pd.merge(forestAge, forest_table)

    # calcaulte emissiosn from disturbances, removals from forests remaining forests,
    # and removals from nonforest to forest using function
    calculate_FRF(forestAge, year1, year2)

    # todo add function to calculate emissions using dictionary

    return groupByLanduseChangeDF.sort_values(by=["Hectares"], ascending=False), forestAge.sort_values(
        by=["Hectares"],
        ascending=False)


# # execute the main function
landuse_result, forestType_result = main(**inputConfig)

# save the results
save_results(landuse_result,forestType_result,outputPath,datetime,startTime)
