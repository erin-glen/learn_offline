import arcpy
import pprint
import os
arcpy.env.overwriteOutput = True

arcpy.env.workspace = r"U:\eglen\Data\PADUS3_0Geodatabase\PAD_US3_0.gdb"
pprint.pprint(arcpy.ListFeatureClasses())
output = r"U:\eglen\Projects\LEARN Tools\Data\SourceData\Data\Rasters\AOI\NF_Selections"

featureclass = r"U:\eglen\Data\PADUS3_0Geodatabase\PAD_US3_0.gdb\PADUS3_0Combined_Proclamation_Marine_Fee_Designation_Easement"
fields = arcpy.ListFields(featureclass)
for field in fields:
    pprint.pprint(field.name)
print(arcpy.GetCount_management(featureclass))

#selection_list = ["Deschutes", "Gifford Pinchot", "Eldorado", "Cherokee", "Chequamegon-Nicolet"]
#for selection in selection_list:
    #selection_name = selection

NF_selection = os.path.join(output,'NF_selections.shp')
outputFolder = os.path.join(output,"NF_AOIs")
arcpy.analysis.Split(NF_selection,NF_selection,'Unit_Nm',outputFolder)

arcpy.env.workspace = outputFolder
datalist = []
datalist = arcpy.ListFeatureClasses()
print(datalist)


