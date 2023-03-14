import arcpy
import pprint
import os
arcpy.env.overwriteOutput = True

arcpy.env.workspace = r"U:\eglen\Data\Scope3_AOIs"
datalist = []
datalist = arcpy.ListFeatureClasses()
print(datalist)


