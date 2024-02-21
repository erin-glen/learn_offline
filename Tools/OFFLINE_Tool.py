# helper functions used in multiple parts of the gp tool
import arcpy
import os
import pandas as pd
from datetime import datetime
from lookups import nlcdParentRollupCategories
from funcs import tabulateAreaByStratification, calculate_category, save_results, \
    calculate_FRF, fillNA, calculate_plantable, calculate_treeCanopy, disturbanceMax, zonal_sum_carbon, \
    landuseStratificationRaster, calculateDisturbances, mergeAgeFactors, calculateFNF, summarize_ghg, \
    summarize_treecanopy, create_matrix, write_dataframes_to_csv

arcpy.env.overwriteOutput = True
# todo add option to select e / r and calculate carbon from TOF

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

    years_list = ["2001", "2004", "2006", "2008", "2011", "2013", "2016", "2019", "2021"]
    year1 = str(input("Year 1: "))
    assert year1 in years_list, f"{year1} is not in the list of valid inputs."
    year2 = str(input("Year 2: "))
    assert year2 in years_list, f"{year2} is not in the list of valid inputs."
    aoi_name = input("Name of AOI? ")
    tree_canopy_list = ["NLCD","CBW", "Local"]
    tree_canopy = input("Tree Canopy source? ")
    assert tree_canopy in tree_canopy_list, f"{tree_canopy} is not in the list of valid inputs."

    cellsize = "30"



    # hardcoded AOI for development - Montgomery County, Maryland
    # aoi = os.path.join(wd, "data", "AOI", "MontgomeryMD.shp")
    aoi_path = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters\AOI"
    aoi = os.path.join(aoi_path, (aoi_name + ".shp"))

    if "NLCD" in tree_canopy:
        treecanopy_path = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters\TreeCanopy\NLCD"
        treecanopy_1=os.path.join(treecanopy_path, "nlcd_tcc_conus_" + year1 + "_v2021-4.tif")
        treecanopy_2=os.path.join(treecanopy_path, "nlcd_tcc_conus_" + year2 + "_v2021-4.tif")

    elif "CBW" in tree_canopy:
        treecanopy_path = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters\TreeCanopy\CBW"
        treecanopy_1 = os.path.join(treecanopy_path, "cbw_2013_treecanopy_Agg30m_int.tif")
        treecanopy_2 = os.path.join(treecanopy_path, "cbw_2018_treecanopy_Agg30m_int.tif") \

    # set paths for custom tree canopy and plantable areas
    else:
        # treecanopy_path = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters\TreeCanopy\Local"
        # filename1 = "Aspen_2016.tif"
        # filename2 = "Aspen_2020.tif"
        # treecanopy_1 = os.path.join(treecanopy_path,aoi_name,filename1)
        # treecanopy_2 = os.path.join(treecanopy_path,aoi_name,filename2)

        treecanopy_1 = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters\TreeCanopy\Local\Aspen\Reclassed\Aspen_2016.tif"
        treecanopy_2 = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters\TreeCanopy\Local\Aspen\Reclassed\Aspen_2020.tif"

    plantableAreas = "None"

    #####ALTER INPUT HERE
    inputConfig = dict(
        nlcd_1=os.path.join(dataFolder, "LandCover", ("NLCD_" + year1 + "_Land_Cover_l48_20210604.tif")),
        nlcd_2=os.path.join(dataFolder, "LandCover", ("NLCD_" + year2 + "_Land_Cover_l48_20210604.tif")),
        forestAgeRaster=os.path.join(dataFolder, "ForestType", "forest_raster_07232020.tif"),
        treecanopy_1= treecanopy_1,
        treecanopy_2= treecanopy_2,
        plantableAreas= plantableAreas,
        carbon_ag_bg_us=os.path.join(dataFolder, "Carbon", "carbon_ag_bg_us.tif"),
        carbon_sd_dd_lt=os.path.join(dataFolder, "Carbon", "carbon_sd_dd_lt.tif"),
        carbon_so=os.path.join(dataFolder, "Carbon", "carbon_so.tif"),
        forest_lookup_csv=os.path.join(dataFolder, "ForestType", "forest_raster_09172020.csv"),
        disturbanceRasters=[
            # os.path.join(dataFolder, "Disturbances", "disturbance_0104.tif"),
            # os.path.join(dataFolder, "Disturbances", "disturbance_0406.tif"),
            # os.path.join(dataFolder, "Disturbances", "disturbance_0608.tif"),
            # os.path.join(dataFolder, "Disturbances", "disturbance_0811.tif"),
            # os.path.join(dataFolder, "Disturbances", "disturbance_1113.tif"),
            # os.path.join(dataFolder, "Disturbances", "disturbance_1316.tif"),
            os.path.join(dataFolder, "Disturbances", "disturbance_1619.tif"),
            #os.path.join(dataFolder, "Disturbances", "disturbance_1921.tif")
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


def main(aoi, nlcd_1, nlcd_2, forestAgeRaster, treecanopy_1, treecanopy_2, plantableAreas, carbon_ag_bg_us,
         carbon_sd_dd_lt, carbon_so, forest_lookup_csv, disturbanceRasters):
    """
    Landuse change stratification summaries for forest age, treecanopy, carbon & disturbance
    :param plantableAreas:
    :param forest_lookup_csv:
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
    # add these results to the output
    treeCover = calculate_plantable(plantableAreas, stratRast, treeCover, aoi, cellsize)

    ###### --------------- Disturbance - tabulate the area -----------------
    arcpy.AddMessage("STEP 3: Cross tabulating disturbance area by stratification class")
    arcpy.AddMessage("Number of disturbance rasters: {}".format(len(disturbanceRasters)))
    disturbance_wide, disturbRast = disturbanceMax(disturbanceRasters, stratRast)

    # Carbon - zonal sum
    arcpy.AddMessage("STEP 4: Zonal statistics sum for carbon rasters by stratification class")
    carbon = zonal_sum_carbon(stratRast, carbon_ag_bg_us, carbon_sd_dd_lt, carbon_so)

    # merge disturbance area, tree cover totals, carbon totals by stratification class
    groupByLanduseChangeDF = carbon.merge(treeCover, how='outer',
                                          on=["StratificationValue", "NLCD1_class", "NLCD2_class"]).merge(
        disturbance_wide, how='outer', on=["StratificationValue", "NLCD1_class", "NLCD2_class"]
    )

    groupByLanduseChangeDF["NLCD_1_ParentClass"] = groupByLanduseChangeDF["NLCD1_class"].map(nlcdParentRollupCategories)
    groupByLanduseChangeDF["NLCD_2_ParentClass"] = groupByLanduseChangeDF["NLCD2_class"].map(nlcdParentRollupCategories)
    groupByLanduseChangeDF['Category'] = groupByLanduseChangeDF.apply(calculate_category, axis=1)
    groupByLanduseChangeDF['Total Emissions Forest to Non Forest CO2'] = groupByLanduseChangeDF.apply(calculateFNF,
                                                                                                      axis=1)

    # Forest Age Type - tabulate the area
    arcpy.AddMessage("STEP 5: Tabulating total area for the forest age types by stratification class")
    forestAge = tabulateAreaByStratification(stratRast, forestAgeRaster, outputName="ForestAgeTypeRegion")

    # merge forestAge area total + disturbance areas by forestAge
    arcpy.AddMessage("STEP 6: Tabulating disturbance area for the forest age types by stratification class for fires")
    forestAge = calculateDisturbances(disturbRast, stratRast, forestAgeRaster, forestAge)

    # fill empty cells with zero for calculations
    forestAge = fillNA(forestAge)

    forestAge = mergeAgeFactors(forestAge, forest_lookup_csv)

    arcpy.AddMessage(
        "STEP 7: Calculating emissions from disturbances and removals from undisturbed forests / non forest to forest")
    # calculate emissions from disturbances, removals from forests remaining forests,
    # and removals from non forest to forest using function
    calculate_FRF(forestAge, year1, year2)

    return groupByLanduseChangeDF.sort_values(by=["Hectares"], ascending=False), forestAge.sort_values(
        by=["Hectares"],
        ascending=False)


# # execute the main function
landuse_result, forestType_result = main(**inputConfig)

# save the results
save_results(landuse_result, forestType_result, outputPath, datetime, startTime)

# summarize the results
years = int(year2) - int(year1)
tc_summary = summarize_treecanopy(landuse_result)
transition_matrix = create_matrix(landuse_result)
ghg_result = summarize_ghg(landuse_result, forestType_result, years)

# save the summary
df_list = [transition_matrix, tc_summary, ghg_result]
csv_file_path = os.path.join(outputPath, "summary.csv")
write_dataframes_to_csv(df_list, csv_file_path, space=5)
