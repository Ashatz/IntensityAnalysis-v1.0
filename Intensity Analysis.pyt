#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Andrew Shatz
#
# Created:     16/03/2014
# Copyright:   (c) Andrew Shatz 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os, sys, shutil
import arcpy
from arcpy import env
from arcpy.sa import *
from arcpy.da import *
import numpy


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [Calculate_Intensity,
                      Insert_Legend,
                      Calculate_Transitions,
                      Create_Raster_Workspace]

class Create_Raster_Workspace(object):
    def __init__(self):
        """Definte the tool (tool name is the name of the class)."""
        self.label = "Create Workspace"
        self.description = "This tool allows users to create the File Geodatabase "\
                         + "from which the intensity analysis will be performed."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        in_dir = arcpy.Parameter(
            displayName = "Input workspace location",
            name = "in_dir",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")

        in_gdb = arcpy.Parameter(
            displayName = "Input Workspace Name",
            name = "in_gdb",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")

        in_fc = arcpy.Parameter(
            displayName = "Input Boundary Feature",
            name = "in_fc",
            datatype = "DEFeatureClass",
            parameterType = "Required",
            direction = "Input")

        in_field = arcpy.Parameter(
            displayName = "Input Boundary Field",
            name = "in_field",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")

        in_field.filter.type = "ValueList"
        in_field.filter.list = [" "]

        in_lulc = arcpy.Parameter(
            displayName = "Input Raster Layers",
            name = "in_lulc",
            datatype = "GPValueTable",
            parameterType = "Required",
            direction = "Input")
        in_lulc.columns = [["Raster Layer", "Input Raster"],['Long', 'Year']]

        params = [in_dir, in_gdb, in_fc, in_field, in_lulc]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        #List fields for input boundary field parameter
        if parameters[2].value != None:
            featureClass = parameters[2].valueAsText
            parameters[3].filter.list = [field.name for field in
                                         arcpy.ListFields(featureClass)]

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        if parameters[4].values != None:
            inputRDs = [parameters[4].values[i][0]
                        for i in range(len(parameters[4].values))]
            sRefs = []
            for RD in inputRDs:
                desc = arcpy.Describe(RD)
                sRefs.append(desc.SpatialReference.name)
            sRefTest = list(set(sRefs))
            if len(sRefTest) > 1:
                parameters[4].setErrorMessage("Spatial references for input raster datasets do not match...")
            elif sRefTest[0] == "Unknown":
                parameters[4].setWarningMessage("Spatial references for input raster datasets is unknown...")

        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        #Step 1: Create file geodatabase workspace

        #----------------------------------------------------------------------#
        messages.addMessage("Creating intensity analysis workspace " +
        "environment...")

        #Declare initial variables
        inDir = parameters[0].valueAsText
        inGdb = parameters[1].valueAsText + ".gdb"
        workspace = os.path.join(inDir, inGdb)
        inFeature = parameters[2].valueAsText
        featureFile = os.path.splitext(
            os.path.basename(inFeature))[0]
        inField = parameters[3].valueAsText
        #Sort input image objects ascending by date
        image_dict = {}
        for i in range(len(parameters[4].values)):
            #messages.addMessage(parameters[1].values[i])
            year = parameters[4].values[i][1]
            #messages.addMessage(year)
            image = parameters[4].values[i][0]
            #messages.addMessage(image)
            image_dict[year] = image
        years = image_dict.keys()
        years.sort()
        images = [image_dict[year] for year in years]

        #Create workspace geodatabase
        try:
            arcpy.CreateFileGDB_management(inDir, inGdb)
        except:
            pass

        #Create boundary feature dataset
        arcpy.env.workspace = workspace
        desc = arcpy.Describe(inFeature)
        sr = desc.SpatialReference
        if "Boundaries" not in arcpy.ListDatasets("", "Feature"):
            arcpy.CreateFeatureDataset_management(workspace, "Boundaries", sr)




        #----------------------------------------------------------------------#

        #Step 2: Create reference tables

        #----------------------------------------------------------------------#
        #Remove reference tables if they exist

        for table in ["tInputs", "tBoundaries"]:
            if table in arcpy.ListTables():
                messages.addMessage("Truncating %s table..."%table)
                arcpy.Delete_management(os.path.join(workspace, table))

        messages.addMessage("Creating file reference table...")

        #Create empty input reference table
        inputTable = os.path.join(workspace, "tInputs")
        arcpy.CreateTable_management(workspace, "tInputs")

        #Create fields
        arcpy.AddField_management(inputTable, "File_Name", "TEXT")
        arcpy.AddField_management(inputTable, "File_Year", "SHORT")

        #Add values to the input table in chronological order
        with arcpy.da.InsertCursor(inputTable, ["File_Name", "File_Year"]) as ic:
            for i in range(len(parameters[4].values)):
                baseImage = os.path.splitext(os.path.basename(str(images[i])))[0]
                ic.insertRow((baseImage, years[i]))

        #Create boundaries table
        boundaryTable = os.path.join(workspace, "tBoundaries")
        arcpy.CreateTable_management(workspace, "tBoundaries")

        #Create fields
        arcpy.AddField_management(boundaryTable, "File_Name", "TEXT")
        arcpy.AddField_management(boundaryTable, "Reference_Field", "TEXT")

        #Insert values into boundaries table
        with arcpy.da.InsertCursor(boundaryTable, ["File_Name", "Reference_Field"]) as ic:
            ic.insertRow((featureFile, inField))

        #----------------------------------------------------------------------#

        #Step 3: Load input files to primary geodatabase

        #----------------------------------------------------------------------#
        messages.addMessage("Loading input file geodatabase...")

        #Format zonal feature class input from shapefile

        """
        arcpy.FeatureClassToFeatureClass_conversion(
            parameters[2].value,
            os.path.join(
                parameters[0].value,
                parameters[1].value + ".gdb",
                "Boundaries"),
            featureFile)
        """
        arcpy.Sort_management(inFeature,
            os.path.join(workspace, "Boundaries", featureFile),
            [[inField, "ASCENDING"]])


        arcpy.RasterToGeodatabase_conversion(
            images,
            os.path.join(
                parameters[0].value,
                parameters[1].value + ".gdb"))

        messages.addMessage("Raster datasets loaded to geodatabase...")

        #Overwrite and reformat raster dataset attribute tables
        arcpy.env.workspace = workspace
        rdSets = arcpy.ListDatasets("", "Raster")
        for rdSet in rdSets:
            arcpy.management.BuildRasterAttributeTable(rdSet, True)

        messages.addMessage("Intensity analysis workspace created...")

        return
        #----------------------------------------------------------------------#

class Insert_Legend(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)"""
        self.label = "Insert Legend Categories"
        self.description = ""
        self.canRunInBackground= False

    def getParameterInfo(self):
        """Define parameter definitions"""

        in_gdb = arcpy.Parameter(
            displayName = "Input Workspace Name",
            name = "in_gdb",
            datatype = "DEWorkSpace",
            parameterType = "Required",
            direction = "Input")

        in_cats = arcpy.Parameter(
            displayName = "Input Land Use/Cover Categories",
            name = "in_cats",
            datatype = "GPValueTable",
            parameterType = "Required",
            direction = "Input")

        in_cats.columns = [["Long", "Value"],["String", "Category"]]
        params = [in_gdb, in_cats]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        #Update second parameter value for value ID drop-down list
        if parameters[0].value != None and parameters[1].values == None:
            workspace = parameters[0].valueAsText
            arcpy.env.workspace = workspace
            #Get the first file name (assuming that all files have the same legend)
            with arcpy.da.SearchCursor(
                os.path.join(workspace, "tInputs"),
                "File_Name") as sc:
                    input = [row[0] for row in sc]
                    input = input[0]
            #Get the values from the file
            with arcpy.da.SearchCursor(
                os.path.join(workspace, input),
                "Value") as sc:
                    values = [row[0] for row in sc]
            del sc #Remove any locks from the cursor
            #Generate default values list
            valueList = [[value, None] for value in values]
            for i in range(len(valueList)):
                if valueList[i][0] == 0:
                    valueList[i][1] = "NoData"
            parameters[1].values = valueList

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        #Define initial variables
        workspace = parameters[0].valueAsText
        arcpy.env.workspace = workspace

        legendCats = parameters[1].values

        #Create tLegend table and insert values
        messages.addMessage("Creating legend reference table...")
        legendTable = os.path.join(workspace, "tLegend")
        if legendTable in arcpy.ListTables():
            arcpy.Delete_management(legendTable)
        arcpy.CreateTable_management(workspace, "tLegend")
        arcpy.AddField_management(legendTable, "Value", "LONG")
        arcpy.AddField_management(legendTable, "Category", "TEXT")
        with arcpy.da.InsertCursor(legendTable, ["Value", "Category"]) as ic:
            for i in range(len(legendCats)):
                ic.insertRow((legendCats[i][0], legendCats[i][1]))


class Calculate_Transitions(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Calculate Transitions"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        #First parameter: identify intensity analysis workspace
        in_gdb = arcpy.Parameter(
            displayName = "Input Workspace Name",
            name = "in_gdb",
            datatype = "DEWorkSpace",
            parameterType = "Required",
            direction = "Input")

        #Second Parameter: Input boundary feature class
        in_bfc = arcpy.Parameter(
            displayName = "Input Boundary Feature Class",
            name = "id_bfc",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")
        in_bfc.filter.type = "ValueList"
        in_bfc.filter.list = [" "]

        #Third parameter: select losing categories
        in_loss = arcpy.Parameter(
            displayName = "Input Losing Categories",
            name = "in_loss",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input",
            multiValue = True)

        in_loss.filter.type = "ValueList"
        in_loss.filter.list = [" "]

        #Fourth parameter: select gaining categories
        in_gain = arcpy.Parameter(
            displayName = "Input Gaining Categories",
            name = "in_gain",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input",
            multiValue = True)

        in_gain.filter.list = [" "]
        in_gain.filter.type = "ValueList"

        opt_scratch = arcpy.Parameter(
            displayName = "Keep Scratch Workspace",
            name = "opt_scratch",
            datatype = "GPBoolean",
            parameterType = "Optional",
            direction = "Input")

        opt_scratch.value = False

        params = [in_gdb, in_bfc, in_loss, in_gain, opt_scratch]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        #Update second parameter value for boundary ID drop-down list
        workspace = parameters[0].valueAsText
        arcpy.env.workspace = workspace
        if parameters[0].value != None:
            #Create drop-down list for second parameter
            with arcpy.da.SearchCursor(
                os.path.join(workspace, "tBoundaries"),
                "File_Name") as sc:
                    parameters[1].filter.list = [row[0] for row in sc]
            del sc
            #Create drop-down list for third and fourth parameters
            with arcpy.da.SearchCursor(
                os.path.join(workspace, "tLegend"),
                "Category") as sc:
                    catList = [row[0] for row in sc]
                    parameters[2].filter.list = catList
                    parameters[3].filter.list = catList

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        #-----Variable preparation------#

        #Prepare strings from input parameters
        workspace = parameters[0].valueAsText
        arcpy.env.workspace = workspace
        boundaries = parameters[1].valueAsText
        boundariesPath = os.path.join(workspace, "Boundaries", boundaries)

        #Define scratch workspace
        a = os.path.splitext(workspace)
        scratch = a[0] + "_scratch" + a[1]

        #Declare boundary polygon/field
        with arcpy.da.SearchCursor(os.path.join(workspace, "tBoundaries"),
            ["File_Name", "Reference_Field"]) as sc:
                boundaryDict = {}
                for row in sc:
                    boundaryDict[row[0]] = row[1]
        del sc #Remove any locks
        boundariesID = boundaryDict[boundaries]

        #Create master input list
        with arcpy.da.SearchCursor(os.path.join(workspace, "tInputs"),
            ["File_Name", "File_Year"]) as sc:
                inputs = [row for row in sc]
                #Create loss/gain list
        del sc #Remove any locks

        loss_list = [os.path.join(workspace, inputs[i][0]) for i in range(len(inputs)-1)]
        gain_list = [os.path.join(workspace, inputs[i][0]) for i in range(1, len(inputs))]
        years = [input[1] for input in inputs]
        intervals = [years[i+1] - years[i] for i in range(len(years)-1)]

        #Create master category lists
        with arcpy.da.SearchCursor(os.path.join(workspace, "tLegend"),
            ["Category", "Value"]) as sc:
                catDict = {}
                for row in sc:
                    catDict[row[0]] = row[1]

        lossCats = parameters[2].values
        gainCats = parameters[3].values

        #-------------------------------#

        #--------Begin analysis---------#

        #Create scratch workspace
        messages.addMessage("Creating scratch workspace...")
        try:
            arcpy.CreateFileGDB_management(
                      os.path.dirname(scratch),
                      os.path.basename(scratch))
        except:
            messages.addMessage("Truncating scratch workspace...")
            arcpy.Delete_management(scratch)
            arcpy.CreateFileGDB_management(
                      os.path.dirname(scratch),
                      os.path.basename(scratch))

        #Define remap values for loss/gain classification
        def create_catRemap(value_cats, transition):
            #Get list of all possible category values from category field
            with arcpy.da.SearchCursor(os.path.join(workspace, "tLegend"),
                ["Value", "Category"]) as sc:
                    cats = [list(row) for row in sc]
                    remapTable = []
                    for cat in cats:
                        if cat[1] in value_cats:
                            if transition == "loss":
                                remapTable.append([str(cat[0]), "1"])
                            elif transition== "gain":
                                remapTable.append([str(cat[0]), "2"])
                        else:
                            remapTable.append([str(cat[0]), "0"])
            del sc
            return RemapValue(remapTable)

        #Define loss/gain image reclassification
        def lgReclass(map_list, map_cats, output, transition):
            lossReclass = Reclassify(map_list, "Value",
                create_catRemap(map_cats, transition),
                "NODATA")
            lossReclass.save(os.path.join(scratch, output))

        messages.addMessage("Generating loss and gain images...")
        #Generate loss and gain images
        for i in range(len(intervals)):
            loss_output = "i" + str(i + 1) + "_loss"
            lgReclass(loss_list[i], lossCats, loss_output, "loss")
            gain_output = "i" + str(i + 1) + "_gain"
            lgReclass(gain_list[i], gainCats, gain_output, "gain")

        messages.addMessage("Generating loss-to-gain images...")
        #Generate loss-to-gain images
        for i in range(len(intervals)):
            in_loss = os.path.join(scratch,"i" + str(i + 1) + "_loss")
            in_gain = os.path.join(scratch,"i" + str(i + 1) + "_gain")
            output = os.path.join(scratch, "i" + str(i + 1) + "_ltg")
            arcpy.Plus_3d(in_loss, in_gain, output)

        #Create and load master area table
        messages.addMessage("Generating master border area table...")

        area_table = os.path.join(workspace, "tAreas")

        try:
            arcpy.CreateTable_management(workspace, "tAreas")
        except:
            messages.addMessage("Border area table already exists...truncating...")
            #This is temporary
            arcpy.Delete_management(area_table)
            messages.addMessage("Recreating area table...")
            arcpy.CreateTable_management(workspace, "tAreas")

        arcpy.AddField_management(area_table, boundariesID, "LONG")

        #Calculate polygon areas
        try:
            arcpy.AddField_management(area_table, "Area_Tot", "FLOAT")
        except:
            pass

        outTable = os.path.join(scratch, "tar_master")

        arcpy.env.workspace= scratch
        if "tar_master" not in arcpy.ListTables():
            messages.addMessage("Calculating polygon areas...")
            arcpy.sa.TabulateArea(boundariesPath, boundariesID,
                loss_list[0], "Value", outTable)
        else:
            messages.addMessage("Overwriting polygon area table...")
            arcpy.management.Delete(outTable)
            messages.addMessage("Recalculating polygon areas...")
            arcpy.sa.TabulateArea(boundaries, boundariesID,
                loss_list[0], "Value", outTable)

        id_field = [field.name for field in arcpy.ListFields(outTable)
                    if boundariesID in field.name]
        tarFields = [field.name for field in arcpy.ListFields(outTable)
                     if "VALUE" in field.name
                     and int(field.name.split("_")[-1]) > 0]

        #Calculate town areas and insert them into master table
        with arcpy.da.SearchCursor(outTable, id_field) as sc:
            townIds = [int(row[0]) for row in sc]
        with arcpy.da.SearchCursor(outTable, tarFields) as sc:
            townArea = [sum(row) for row in sc]

        del sc

        #Insert Town IDs and total areas into master area table
        with arcpy.da.InsertCursor(area_table, [boundariesID, "Area_Tot"]) as ic:
            for i in range(len(townIds)):
                ic.insertRow((townIds[i],townArea[i]))
        del ic

        #Define function to:
        #   1) Insert new fields into the master table for each interval
        #   2) Produce tabular area of loss, gain, and loss-to-gain per interval
        #   3) Insert results into fields produced.

        def insertAreas(interval):
            #Step 0: Define local variables
            loss = "i" + str(interval) + "_loss"
            gain = "i" + str(interval) + "_gain"
            ltg = "i" + str(interval) + "_ltg"
            #Tabular area table - tar
            tar = "i" + str(interval) + "_tar"
            tempList = [loss, gain, ltg]

            #Step 1:
            [arcpy.AddField_management(
                os.path.join(workspace, area_table), temp, "FLOAT")
                for temp in tempList]

            #Step 2:
            outTar = os.path.join(scratch, tar)
            try:
                arcpy.sa.TabulateArea(boundariesPath,
                    boundariesID,
                    os.path.join(scratch, ltg),
                    "Value",
                    outTar)
            except:
                arcpy.Delete_management(outTar)
                arcpy.sa.TabulateArea(boundariesPath,
                    boundariesID,
                    os.path.join(scratch, ltg),
                    "Value",
                    outTar)

            fieldList = [field.name for field in arcpy.ListFields(outTar, "VALUE_*")
                         if field.name[-1] in ["1", "2", "3"]]


            with arcpy.da.SearchCursor(outTar, fieldList) as sc:
                valueList = [list(row) for row in sc]
            del sc

            with arcpy.da.UpdateCursor(area_table, tempList) as uc:
                i = 0
                for row in uc:
                    #Loss and gain are part of index 2 and must be considered
                    #in the area calculation to avoid high intensity bias
                    row[0] = valueList[i][0] + valueList[i][2]
                    row[1] = valueList[i][1] + valueList[i][2]
                    row[2] = valueList[i][2]
                    uc.updateRow(row)
                    i += 1
            del uc

        messages.addMessage("Calculating interval losses and gains...")
        [insertAreas(i) for i in range(1, len(years))]

        #Remove intermediary files
        #shutil.rmtree(workspace)
        if parameters[4].value == False:
            arcpy.Delete_management(scratch)

        return


class Calculate_Intensity(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Calculate Intensity"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        #Parameter 1: input intensity workspace
        in_gdb = arcpy.Parameter(
            displayName = "Input Workspace",
            name = "in_gdb",
            datatype = "DEWorkspace",
            parameterType = "Required",
            direction = "Input")


        in_bfc = arcpy.Parameter(
            displayName = "Input Boundary Feature Class",
            name = "in_bfc",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")
        in_bfc.filter.type = "ValueList"


        in_tar = arcpy.Parameter(
            displayName = "Input Area Table",
            name = "in_tar",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")

        in_tar.filter.type = "ValueList"

        """
        int_type = arcpy.Parameter(
            displayName = "Intensity Calculation",
            name = "int_type",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")

        int_type.filter.type = "ValueList"
        int_type.filter.list = ["SPATIAL", "TEMPORAL", "SPATIOTEMPORAL", "ALL"]
        int_type.value = "ALL"
        """

        params = [in_gdb, in_bfc, in_tar]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        #Update area table parameter values for drop-down list
        workspace = parameters[0].valueAsText
        if workspace != None:
            arcpy.env.workspace = workspace

            #Create boundary file drop down list
            with arcpy.da.SearchCursor(
                os.path.join(workspace,"tBoundaries"), "File_Name") as sc:
                parameters[1].filter.list = [row[0] for row in sc]
            del sc

            #Create area table drop down list
            inTar = arcpy.ListTables("tAreas*")
            parameters[2].filter.list = [table for table in inTar]
        elif workspace == None:
            parameters[1].filter.list = [" "]
            parameters[2].filter.list = [" "]
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        #Step 0: Declare initial variables
        messages.addMessage("Declaring input variables...")
        workspace = parameters[0].valueAsText

        #Get boundary field
        boundaryFc = parameters[1].valueAsText
        with arcpy.da.SearchCursor(
            os.path.join(workspace, "tBoundaries"),
            ["File_Name", "Reference_Field"]) as sc:
                boundaryField = [row[1] for row in sc if row[0] == boundaryFc][0]
        del sc
        with arcpy.da.SearchCursor(
            os.path.join(workspace, boundaryFc),
            boundaryField) as sc:
                boundaryIDs = [row[0] for row in sc]

        #Declare intervals
        inputTable = os.path.join(workspace, "tInputs")
        years = [row.File_Year for row in arcpy.SearchCursor(inputTable)]
        intervals = [years[i+1] - years[i] for i in range(len(years)-1)]
        intervalSum = sum(intervals)

        #Declare loss/gain/loss-to-gain matrices
        areaTable = os.path.join(workspace, parameters[2].valueAsText)
        areaFields = [field.name for field in
                      arcpy.ListFields(areaTable)]
        #Loss matrix
        lossFields = [field for field in areaFields if "_loss" in field]
        with arcpy.da.SearchCursor(areaTable, lossFields) as sc:
            lossMatrix = [list(row) for row in sc]
        del sc

        #Gain matrix
        gainFields = [field for field in areaFields if "gain" in field]
        with arcpy.da.SearchCursor(areaTable, gainFields) as sc:
            gainMatrix = [list(row) for row in sc]
        del sc

        #Loss-to-gain matrix
        ltgFields = [field for field in areaFields if "ltg" in field]
        with arcpy.da.SearchCursor(areaTable, ltgFields) as sc:
            ltgMatrix = [list(row) for row in sc]
        del sc


        #Step 1: Design intensity analysis schema.
        arcpy.env.workspace = workspace
        intensityTable = os.path.join(workspace, "tIntensity")

        #Create table
        if "tIntensity" in arcpy.ListTables():
            messages.addMessage("Truncating current intensity table schema...")
            arcpy.Delete_management(intensityTable)
        messages.addMessage("Creating intensity table schema...")
        arcpy.CreateTable_management(workspace, "tIntensity")
        #Add boundary ID field and insert values
        arcpy.AddField_management(intensityTable, boundaryField, "LONG")
        #Insert values into ID field
        with arcpy.da.InsertCursor(intensityTable, boundaryField) as ic:
            for ID in boundaryIDs:
                ic.insertRow((ID,))

        #Add intensity fields to global intensity table
        messages.addMessage("Adding appropriate fields to table...")
        fieldList = []
        directions = ["loss", "gain"]
        for direction in directions:
            #Create intensity fields
            [fieldList.append("i" + str(i) + "t_" + direction)
             for i in range(1, len(years))]
            fieldList.append("tu_" + direction)
            [fieldList.append("i" + str(i) + "tr_" + direction)
             for i in range(1, len(years))]
            fieldList.append("tsr_" + direction)
        [arcpy.AddField_management(intensityTable, field, "FLOAT")
         for field in fieldList]

        #Calculate intensity values
        messages.addMessage("Calculating intensities...")
        inputMatrix = []
        for direction in directions:
            inputFieldList = [field.name for field in arcpy.ListFields(intensityTable)
                              if "t_" + direction in field.name]

            def calculateIntesityRows(direction):
                output = []
                for i in range(len(boundaryIDs)):
                    row = []
                    for interval in intervals:
                        intervalIndex = intervals.index(interval)
                        if direction == "loss":
                            row.append(ltgMatrix[i][intervalIndex]/
                                      (lossMatrix[i][intervalIndex]*interval))
                        elif direction == "gain":
                            row.append(ltgMatrix[i][intervalIndex]/
                                      (gainMatrix[i][intervalIndex]*interval))
                    output.append(row)
                return output

            intensityMatrix = calculateIntesityRows(direction)
            with arcpy.da.UpdateCursor(intensityTable, inputFieldList) as uc:
                i = 0
                for row in uc:
                    for j in range(len(inputFieldList)):
                        row[j] = intensityMatrix[i][j]*100.0
                    uc.updateRow(row)
                    i += 1
            del uc

            #Calculate uniform temporal intensity
            if direction == "loss":
                uniformTIntensity = [sum(ltgMatrix[i])/(intervalSum * lossMatrix[i][0])
                                     for i in range(len(ltgMatrix))]
            elif direction == "gain":
                uniformTIntensity = [sum(ltgMatrix[i])/(intervalSum * gainMatrix[i][-1])
                                     for i in range(len(ltgMatrix))]
            with arcpy.da.UpdateCursor(intensityTable, "tu_" + direction) as uc:
                i = 0
                for row in uc:
                    #Scale into percentages
                    row[0] = uniformTIntensity[i]*100
                    uc.updateRow(row)
                    i += 1
            del uc

            #Calculate temporal ratio
            inputRatioFields = [field.name for field in arcpy.ListFields(intensityTable)
                                if "tr_" + direction in field.name]
            uniformMatrix = []
            with arcpy.da.SearchCursor(intensityTable, "tu_"+ direction) as sc:
                uniformList = [row[0] for row in sc]
                [uniformMatrix.append(uniformList) for i in range(len(intervals))]
                uniformMatrix = numpy.matrix(uniformMatrix)
                uniformTranspose = numpy.matrix.transpose(uniformMatrix)
            del sc
            ratioMatrix = intensityMatrix/uniformTranspose
            with arcpy.da.UpdateCursor(intensityTable, inputRatioFields) as uc:
                i = 0
                for row in uc:
                    for j in range(len(inputRatioFields)):
                        row[j] = ratioMatrix.tolist()[i][j]*100
                        uc.updateRow(row)
                    i += 1
            del uc

            #Calculate global-local spatio-temporal uniformity
            #Gather area variables
            with arcpy.da.SearchCursor(os.path.join(workspace,"tAreas"),
                 "Area_Tot") as sc:
                boundaryArea = [row[0] for row in sc]
                boundarySum = sum(boundaryArea)
            del sc
            #Gather global spatial intensity
            globalUniform = sum([(uniformList[i] * (boundaryArea[i]/boundarySum))
                             for i in range(len(uniformList))])
            with arcpy.da.UpdateCursor(intensityTable, "tsr_" + direction) as uc:
                i = 0
                for row in uc:
                    row[0] = uniformList[i]/globalUniform
                    i += 1
                    uc.updateRow(row)
            del uc





# # # # # Notes # # # # #
#
# Import boundary file into geodatabase using Sort tool (Management) and add
# parameter for input boundary field.