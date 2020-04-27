import arcpy
import os
import time
import json
from . import hydro
from . import tasker

arcpy.CheckOutExtension("Spatial")

dirname = os.path.dirname(__file__)
directory = os.path.join(dirname, '../static/data')
sr_map = arcpy.SpatialReference(4326)
sr_md = arcpy.SpatialReference(26985)
modifieddt_ = "April 9, 2020"
thomasversion = "2020"

def coord_transf(x1,y1,sr1,sr2):
    pt1 = arcpy.Point(x1,y1)
    ptgeo1 = arcpy.PointGeometry(pt1, sr1)
    ptgeo2 = ptgeo1.projectAs(sr2)
    pt2 = ptgeo2.lastPoint
    x2 = float(pt2.X)
    y2 = float(pt2.Y)
    return x2, y2

def data_selection(request):

    coord_0 = float(request.GET['cw'])
    coord_1 = float(request.GET['cs'])
    coord_2 = float(request.GET['ce'])
    coord_3 = float(request.GET['cn'])
    proj_name = str(request.GET['proj']).replace(" ", "_")
    dem_layer = str(request.GET['dem'])
    soil_layer = str(request.GET['soil'])
    land_layer = str(request.GET['lu'])
    hyd_cond = str(request.GET['hyd'])
    acc_thr = int(request.GET['acc'])
    burnopt = request.GET.get('burn') == 'true'

    temp_folder = os.path.join(directory, "example")

    folder_name = time.strftime("%Y%m%d_%H%M%S") + "_" + proj_name
    optfolder = os.path.join(temp_folder, folder_name)
    os.makedirs(optfolder)

    file = open(os.path.join(optfolder, "debug.txt"),"w")
    file.write(str(burnopt))
    file.close()

    arcpy.env.overwriteOutput = True
    arcpy.env.scratchWorkspace = optfolder
    arcpy.env.workspace = optfolder

    coord_t0, coord_t1 = coord_transf(coord_0, coord_1, sr_map, sr_md)
    coord_t2, coord_t3 = coord_transf(coord_2, coord_3, sr_map, sr_md)
    coord_t0 = round(coord_t0,0)
    coord_t1 = round(coord_t1,0)
    coord_t2 = round(coord_t2,0)
    coord_t3 = round(coord_t3,0)

    extent_clip = arcpy.Extent(coord_t0, coord_t1, coord_t2, coord_t3)
    arcpy.env.extent = extent_clip

    array = arcpy.Array()
    array.add(arcpy.Point(extent_clip.XMin, extent_clip.YMin))
    array.add(arcpy.Point(extent_clip.XMin, extent_clip.YMax))
    array.add(arcpy.Point(extent_clip.XMax, extent_clip.YMax))
    array.add(arcpy.Point(extent_clip.XMax, extent_clip.YMin))
    array.add(arcpy.Point(extent_clip.XMin, extent_clip.YMin))
    mask_layer = arcpy.Polygon(array,sr_md)

    dem = os.path.join(directory,"rasters/dems/" + dem_layer)
    dem_clip = os.path.join(optfolder, "dem_clip.tif")
    landuse = os.path.join(directory, "rasters/landuse/" + land_layer)
    soils = os.path.join(directory, "rasters/soils/" + soil_layer)
    nhd = os.path.join(directory, "shp/md_streams/nhd_streamsm.shp")
    streets = os.path.join(directory, "shp/md_roads/mjr-rdsstpm.shp")
    arcpy.env.snapRaster = dem

    mask = os.path.join(optfolder, "mask.shp")
    arcpy.CopyFeatures_management(mask_layer, mask)
    arcpy.Clip_management(dem, "", dem_clip, mask)
    arcpy.Clip_management(dem, "", os.path.join(optfolder, "dem_clip2"), mask)

    arcpy.Clip_analysis(nhd, mask, os.path.join(optfolder, "nhd_strs.shp"))
    arcpy.Clip_analysis(streets, mask, os.path.join(optfolder, "roads.shp"))
    arcpy.PolylineToRaster_conversion(os.path.join(optfolder, "nhd_strs.shp"), "FID", os.path.join(optfolder, "nhd_rast.tif"), "", "", dem)
    arcpy.Clip_management(landuse, "", os.path.join(optfolder, "landuse.tif"), mask,"","ClippingGeometry","MAINTAIN_EXTENT")
    arcpy.Clip_management(soils, "", os.path.join(optfolder, "soils.tif"), mask,"","ClippingGeometry","MAINTAIN_EXTENT")
    arcpy.SymDiff_analysis(os.path.join(directory,"shp/md_extent/big_mask.shp"), mask, os.path.join(optfolder, "big_mask.shp"))

    data_selection2(optfolder, land_layer, soil_layer, hyd_cond)

    burn = 0
    if burnopt:
        burn = 300

    calc1 = arcpy.sa.IsNull(os.path.join(optfolder, "nhd_rast.tif"))
    calc2 = arcpy.sa.Minus(calc1, 1)
    calc1 = arcpy.sa.Times(calc2, -1 * burn)
    calc2 = arcpy.sa.Minus(os.path.join(optfolder, "dem_clip.tif"), calc1)
    burned = arcpy.sa.Plus(calc2, burn)
    burned.save(os.path.join(optfolder, "burned.tif"))

    fill = arcpy.sa.Fill(os.path.join(optfolder, "burned.tif"))
    fill.save(os.path.join(optfolder, "fill.tif"))

    flowdir = arcpy.sa.FlowDirection(fill)
    flowdir.save(os.path.join(optfolder, "flowdir.tif"))

    flowacc = arcpy.sa.FlowAccumulation(flowdir)
    flowacc.save(os.path.join(optfolder, "flowacc.tif"))


    infstr = arcpy.sa.Con(flowacc >= acc_thr, flowacc)
    infstr.save(os.path.join(optfolder, "infstreams.tif"))

    arcpy.Project_management(os.path.join(optfolder, "big_mask.shp"), os.path.join(optfolder, "mask_proj.shp"), sr_map)
    arcpy.Project_management(os.path.join(optfolder, "nhd_strs.shp"), os.path.join(optfolder, "nhd_proj.shp"), sr_map)
    arcpy.Project_management(os.path.join(optfolder, "roads.shp"), os.path.join(optfolder, "roads_proj.shp"), sr_map)

    arcpy.Delete_management(os.path.join(optfolder, "burned.tif"))
    arcpy.Delete_management(os.path.join(optfolder, "fill.tif"))
    arcpy.Delete_management(os.path.join(optfolder, "big_mask.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "nhd_strs.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "roads.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "nhd_rast.tif"))
    arcpy.Delete_management(os.path.join(optfolder, "dir_shift_8"))

    response = {
        'folder': r"static/data/example/" + folder_name,
        'path': optfolder
    }

    return(json.dumps(response))

def showstreams(request):

    optfolder = str(request.GET['path'])
    arcpy.env.overwriteOutput = True
    arcpy.env.scratchWorkspace = optfolder
    arcpy.env.workspace = optfolder

    dem = arcpy.Raster(os.path.join(optfolder, "infstreams.tif"))
    regpoints = dem.extent
    regx, regy = coord_transf(regpoints.XMin, regpoints.YMin, sr_md, sr_map)
    regpoints = str(regx) + " " + str(regy)
    arcpy.ProjectRaster_management(os.path.join(optfolder, "infstreams.tif"), os.path.join(optfolder, "str_proj.tif"), sr_map, "", "", "", regpoints, sr_md)

def delineation(request):

    optfolder = str(request.GET['path'])
    arcpy.env.overwriteOutput = True
    arcpy.env.scratchWorkspace = optfolder
    arcpy.env.workspace = optfolder
    arcpy.env.snapRaster = os.path.join(optfolder, "dem_clip.tif")

    wshed = arcpy.sa.Watershed(os.path.join(optfolder, "flowdir.tif"), os.path.join(optfolder, "pour_point.tif"))
    watershed = arcpy.sa.Con(wshed >= 0, 1, arcpy.sa.IsNull(wshed))
    watershed.save(os.path.join(optfolder, "basingrid.tif"))
    arcpy.RasterToPolygon_conversion(watershed, os.path.join(optfolder, "wshed.shp"),"NO_SIMPLIFY","VALUE")
    arcpy.Project_management (os.path.join(optfolder, "wshed.shp"), os.path.join(optfolder, "wshed_proj.shp"), sr_map)

    # *******************************************************************************************************
    # gage checking and selection
    # *******************************************************************************************************

    usgsgages = os.path.join(directory, "shp/md_gages/usgsgagesm.shp")
    mdgages = os.path.join(directory, "shp/md_gages/mdgagedstreams2016.shp")
    outletpoint = os.path.join(optfolder, "pour_point.shp")
    gagefound = hydro.CheckGages(usgsgages, mdgages, outletpoint, os.path.join(optfolder, "mask.shp"), optfolder)

    # Creates layers for delineating subwatersheds in the future
    arcpy.CreateFeatureclass_management(optfolder, "addasstreams.shp", "POLYLINE", "", "ENABLED", "DISABLED", sr_md)
    arcpy.CreateFeatureclass_management(optfolder, "addasoutlets.shp", "POINT", "", "ENABLED", "DISABLED", sr_md)
    arcpy.CreateFeatureclass_management(optfolder, "addasreservoir.shp", "POINT", "", "ENABLED", "DISABLED", sr_md)

    arcpy.Delete_management(os.path.join(optfolder, "mask_proj.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "nhd_proj.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "roads_proj.shp"))

    response = json.dumps(gagefound)
    return(response)

def delineation_check(request):

    optfolder = str(request.GET['path'])
    mouse_lat_proj = float(request.GET['mouse_lat'])
    mouse_lon_proj = float(request.GET['mouse_lon'])
    x, y = coord_transf(mouse_lon_proj, mouse_lat_proj, sr_map, sr_md)

    arcpy.env.overwriteOutput = True
    arcpy.env.scratchWorkspace = optfolder
    arcpy.env.workspace = optfolder
    arcpy.env.snapRaster = os.path.join(optfolder, "dem_clip.tif")

    rast = arcpy.Raster(os.path.join(optfolder, "infstreams.tif"))
    cellsize = rast.meanCellWidth
    point = arcpy.PointGeometry(arcpy.Point(x, y), sr_md)
    outSnapPour = arcpy.sa.SnapPourPoint(point, os.path.join(optfolder, "infstreams.tif"), cellsize)
    pour_point = arcpy.RasterToPoint_conversion(outSnapPour, os.path.join(optfolder, "pour_point.shp"))
    for row in arcpy.da.SearchCursor(pour_point, ["SHAPE@XY"]):
        x, y = row[0]
    outletxy = arcpy.management.GetCellValue(os.path.join(optfolder, "infstreams.tif"), "{} {}".format(x, y))
    if outletxy.getOutput(0).isnumeric():
        outSnapPour.save(os.path.join(optfolder, "pour_point.tif"))
        return(True)
    else:
        return(False)

def data_selection2(optfolder, landuse, soil, hyd):

    # *******************************************************************************************************
    # Select, clip, and generate CN using soil and landuse rasters
    # *******************************************************************************************************

    arcpy.env.overwriteOutput = True
    arcpy.env.scratchWorkspace = optfolder
    arcpy.env.workspace = optfolder
    arcpy.env.snapRaster = os.path.join(optfolder, "dem_clip.tif")

    landuse_path = os.path.join(optfolder, "landuse.tif")
    soils_path = os.path.join(optfolder, "soils.tif")
    cn_path = os.path.join(optfolder, "curvenumber.tif")

    # Ragan data
    if soil == "ragan":

        # Landuse and hydrologic condition for CN calculation
        if (landuse == "nlcd2011" and hyd == "Fair"):
            outGrid = hydro.nlcdlookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "nlcd2011" and hyd == "Good"):
            outGrid = hydro.nlcdlookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "nlcd2011" and hyd == "Poor"):
            outGrid = hydro.nlcdlookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "nlcd2006" and hyd == "Fair"):
            outGrid = hydro.nlcdlookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "nlcd2006" and hyd == "Good"):
            outGrid = hydro.nlcdlookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "nlcd2006" and hyd == "Poor"):
            outGrid = hydro.nlcdlookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "nlcd2001" and hyd == "Fair"):
            outGrid = hydro.nlcdlookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "nlcd2001" and hyd == "Good"):
            outGrid = hydro.nlcdlookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "nlcd2001" and hyd == "Poor"):
            outGrid = hydro.nlcdlookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "lu2010" and hyd == "Fair"):
            outGrid = hydro.andlookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "lu2010" and hyd == "Good"):
            outGrid = hydro.andlookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "lu2010" and hyd == "Poor"):
            outGrid = hydro.andlookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "mdplu2002" and hyd == "Fair"):
            outGrid = hydro.andlookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "mdplu2002" and hyd == "Good"):
            outGrid = hydro.andlookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "mdplu2002" and hyd == "Poor"):
            outGrid = hydro.andlookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "lu97m" and hyd == "Fair"):
            outGrid = hydro.andlookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "lu97m" and hyd == "Good"):
            outGrid = hydro.andlookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "lu97m" and hyd == "Poor"):
            outGrid = hydro.andlookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "mdde2002" and hyd == "Fair"):
            outGrid = hydro.mddelookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "mdde2002" and hyd == "Good"):
            outGrid = hydro.mddelookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "mdde2002" and hyd == "Poor"):
            outGrid = hydro.mddelookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "luult" and hyd == "Fair"):
            outGrid = hydro.zoninglookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "luult" and hyd == "Good"):
            outGrid = hydro.zoninglookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "luult" and hyd == "Poor"):
            outGrid = hydro.zoninglookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "mrlc" and hyd == "Fair"):
            outGrid = hydro.mrlclookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "mrlc" and hyd == "Good"):
            outGrid = hydro.mrlclookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "mrlc" and hyd == "Poor"):
            outGrid = hydro.mrlclookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "lu70" and hyd == "Fair"):
            outGrid = hydro.usgslookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "lu70" and hyd == "Good"):
            outGrid = hydro.usgslookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "lu70" and hyd == "Poor"):
            outGrid = hydro.usgslookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

    # SSURGO data
    elif soil == "ssurgo_old" or soil == "ssurgo_2018":

        # Landuse and hydrologic condition for CN calculation
        if (landuse == "nlcd2011" and hyd == "Fair"):
            outGrid = hydro.nlcdlookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "nlcd2011" and hyd == "Good"):
            outGrid = hydro.nlcdlookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "nlcd2011" and hyd == "Poor"):
            outGrid = hydro.nlcdlookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "nlcd2006" and hyd == "Fair"):
            outGrid = hydro.nlcdlookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "nlcd2006" and hyd == "Good"):
            outGrid = hydro.nlcdlookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "nlcd2006" and hyd == "Poor"):
            outGrid = hydro.nlcdlookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "nlcd2001" and hyd == "Fair"):
            outGrid = hydro.nlcdlookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "nlcd2001" and hyd == "Good"):
            outGrid = hydro.nlcdlookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "nlcd2001" and hyd == "Poor"):
            outGrid = hydro.nlcdlookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "lu2010" and hyd == "Fair"):
            outGrid = hydro.andlookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "lu2010" and hyd == "Good"):
            outGrid = hydro.andlookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "lu2010" and hyd == "Poor"):
            outGrid = hydro.andlookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "mdplu2002" and hyd == "Fair"):
            outGrid = hydro.andlookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "mdplu2002" and hyd == "Good"):
            outGrid = hydro.andlookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "mdplu2002" and hyd == "Poor"):
            outGrid = hydro.andlookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "lu97m" and hyd == "Fair"):
            outGrid = hydro.andlookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "lu97m" and hyd == "Good"):
            outGrid = hydro.andlookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "lu97m" and hyd == "Poor"):
            outGrid = hydro.andlookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "mdde2002" and hyd == "Fair"):
            outGrid = hydro.mddelookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "mdde2002" and hyd == "Good"):
            outGrid = hydro.mddelookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "mdde2002" and hyd == "Poor"):
            outGrid = hydro.mddelookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "luult" and hyd == "Fair"):
            outGrid = hydro.zoninglookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "luult" and hyd == "Good"):
            outGrid = hydro.zoninglookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "luult" and hyd == "Poor"):
            outGrid = hydro.zoninglookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "mrlc" and hyd == "Fair"):
            outGrid = hydro.mrlclookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "mrlc" and hyd == "Good"):
            outGrid = hydro.mrlclookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "mrlc" and hyd == "Poor"):
            outGrid = hydro.mrlclookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

        if (landuse == "lu70" and hyd == "Fair"):
            outGrid = hydro.usgslookupfair(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "lu70" and hyd == "Good"):
            outGrid = hydro.usgslookupgood(landuse_path, soils_path)
            outGrid.save(cn_path)
        elif (landuse == "lu70" and hyd == "Poor"):
            outGrid = hydro.usgslookuppoor(landuse_path, soils_path)
            outGrid.save(cn_path)

def basincomp(request):

    optfolder = str(request.GET['path'])
    landuse = str(request.GET['lu'])
    hyd = str(request.GET['hyd'])

    arcpy.env.overwriteOutput = True
    arcpy.env.scratchWorkspace = optfolder
    arcpy.env.workspace = optfolder
    arcpy.env.snapRaster = os.path.join(optfolder, "dem_clip.tif")

    landuse_path = os.path.join(optfolder, "landuse.tif")
    soils_path = os.path.join(optfolder, "soils.tif")
    wshed_path = os.path.join(optfolder, "basingrid.tif")

    ## extract by mask to extent of watershed
    lu_out = arcpy.sa.ExtractByMask(landuse_path, wshed_path)
    soil_out = arcpy.sa.ExtractByMask(soils_path, wshed_path)

    ## convert all clipped rasters to polygons
    lu_poly = arcpy.RasterToPolygon_conversion(lu_out, os.path.join(optfolder, "lu_poly.shp"), "NO_SIMPLIFY", "VALUE")
    soil_poly = arcpy.RasterToPolygon_conversion(soil_out, os.path.join(optfolder, "soil_poly.shp"), "NO_SIMPLIFY", "VALUE")

    ## intersect land use and soil to prepare two polygons: "lu_soil" and "lu_cn"
    arcpy.Intersect_analysis([lu_poly, soil_poly], os.path.join(optfolder, "lu_soil.shp"), "ALL", "#", "INPUT")

    ## dissolve above intersected polygon
    arcpy.Dissolve_management(os.path.join(optfolder, "lu_soil.shp"), os.path.join(optfolder, "lu_soil_diss.shp"), "GRIDCODE;GRIDCODE_1", "#", "MULTI_PART", "DISSOLVE_LINES")

    ## add filed to both of above dissolved polygons and compute area in acres
    if not len(arcpy.ListFields(os.path.join(optfolder, "lu_soil_diss.shp"), "area")) > 0:
        arcpy.AddField_management(os.path.join(optfolder, "lu_soil_diss.shp"), "area", "FLOAT", 15, 4)
    arcpy.CalculateField_management(os.path.join(optfolder, "lu_soil_diss.shp"), "area", "!shape.area@acres!", "PYTHON")

    # prepre a list of lu codes to feed into lu_description function in order to obtain matching descriptions list
    lu_match = []
    sc = arcpy.SearchCursor(lu_out, "", "", "VALUE", "")
    for i in sc:
        v = i.getValue("VALUE")
        lu_match.append(v)

    # create list of lists with zeroes
    soil_acre_lists = [[0, 0, 0, 0] for i in range(len(lu_match))]

    # preapre a list of soil acreage using lu_match list
    lc_soil_diss = []
    soil_lc_diss = []
    lc_soil_aa = []
    lu_soil_sc = arcpy.SearchCursor(os.path.join(optfolder, "lu_soil_diss.shp"), "", "", "GRIDCODE;GRIDCODE_1;area", "")
    for s in lu_soil_sc:
        lc = s.getValue("GRIDCODE")
        lc_soil_diss.append(lc)
        sc = s.getValue("GRIDCODE_1")
        soil_lc_diss.append(sc)
        aa = s.getValue("area")
        lc_soil_aa.append(round(aa, 2))

    for idx, lu in enumerate(lu_match):
        for l, s, a in zip(lc_soil_diss, soil_lc_diss, lc_soil_aa):
            if l == lu:
                soil_acre_lists[idx][int(s) - 1] = a

    # prepare matching list of lu description using lu codes from lu raster of watershed
    if landuse == "nlcd2011":
        if hyd == "Fair":
            lut_file = os.path.join(directory, "lookup/nlcdlookupfair.txt")
        elif hyd == "Good":
            lut_file = os.path.join(directory, "lookup/nlcdlookupgood.txt")
        elif hyd == "Poor":
            lut_file = os.path.join(directory, "lookup/nlcdlookuppoor.txt")

    if landuse == "nlcd2006":
        if hyd == "Fair":
            lut_file = os.path.join(directory, "lookup/nlcdlookupfair.txt")
        elif hyd == "Good":
            lut_file = os.path.join(directory, "lookup/nlcdlookupgood.txt")
        elif hyd == "Poor":
            lut_file = os.path.join(directory, "lookup/nlcdlookuppoor.txt")

    if landuse == "nlcd2001":
        if hyd == "Fair":
            lut_file = os.path.join(directory, "lookup/nlcdlookupfair.txt")
        elif hyd == "Good":
            lut_file = os.path.join(directory, "lookup/nlcdlookupgood.txt")
        elif hyd == "Poor":
            lut_file = os.path.join(directory, "lookup/nlcdlookuppoor.txt")

    if landuse == "lu97m":
        if hyd == "Fair":
            lut_file = os.path.join(directory, "lookup/andlookupfair.txt")
        elif hyd == "Good":
            lut_file = os.path.join(directory, "lookup/andlookupgood.txt")
        elif hyd == "Poor":
            lut_file = os.path.join(directory, "lookup/andlookuppoor.txt")

    if landuse == "mdplu2002":
        if hyd == "Fair":
            lut_file = os.path.join(directory, "lookup/andlookupfair.txt")
        elif hyd == "Good":
            lut_file = os.path.join(directory, "lookup/andlookupgood.txt")
        elif hyd == "Poor":
            lut_file = os.path.join(directory, "lookup/andlookuppoor.txt")

    if landuse == "lu2010":
        if hyd == "Fair":
            lut_file = os.path.join(directory, "lookup/andlookupfair.txt")
        elif hyd == "Good":
            lut_file = os.path.join(directory, "lookup/andlookupgood.txt")
        elif hyd == "Poor":
            lut_file = os.path.join(directory, "lookup/andlookuppoor.txt")

    if landuse == "mdde2002":
        if hyd == "Fair":
            lut_file = os.path.join(directory, "lookup/mddelookupfair.txt")
        elif hyd == "Good":
            lut_file = os.path.join(directory, "lookup/mddelookupgood.txt")
        elif hyd == "Poor":
            lut_file = os.path.join(directory, "lookup/mddelookuppoor.txt")

    if landuse == "luult":
        if hyd == "Fair":
            lut_file = os.path.join(directory, "lookup/zoninglookupfair.txt")
        elif hyd == "Good":
            lut_file = os.path.join(directory, "lookup/zoninglookupgood.txt")
        elif hyd == "Poor":
            lut_file = os.path.join(directory, "lookup/zoninglookuppoor.txt")

    if landuse == "mrlc":
        if hyd == "Fair":
            lut_file = os.path.join(directory, "lookup/mrlclookupfair.txt")
        elif hyd == "Good":
            lut_file = os.path.join(directory, "lookup/mrlclookupgood.txt")
        elif hyd == "Poor":
            lut_file = os.path.join(directory, "lookup/mrlclookuppoor.txt")

    if landuse == "lu70":
        if hyd == "Fair":
            lut_file = os.path.join(directory, "lookup/usgslookupfair.txt")
        elif hyd == "Good":
            lut_file = os.path.join(directory, "lookup/usgslookupgood.txt")
        elif hyd == "Poor":
            lut_file = os.path.join(directory, "lookup/usgslookuppoor.txt")

    # run hydro function to obtain land use description of categories present in watershed
    lu_desc = hydro.lu_description(lut_file, lu_match)

    # sum list of lists separately and cat at the end of lu description
    total_area = [round(sum(i),2) for i in zip(*soil_acre_lists)]

    # loop over land use, related total acreage, percent of land covered by this lu category, and A-B-C-D curve numbers
    curve_num = []
    for l in lu_match:
        with open(lut_file, "r") as f:
            next(f)
            for line in f:
                luc = line.split("\t")[0]
                if int(l) == int(luc):
                    temp = []
                    A = line.split("\t")[2]  # CN A
                    temp.append(A)
                    B = line.split("\t")[3]  # CN B
                    temp.append(B)
                    C = line.split("\t")[4]  # CN C
                    temp.append(C)
                    D = line.split("\t")[5]  # CN D
                    temp.append(D)
                    curve_num.append(temp)

    # sum areas for each sub-list individually
    acres = [round(sum(i),2) for i in soil_acre_lists]
    total_all = sum(total_area)
    percent = [round(float(ac_x / total_all) * 100, 2) for ac_x in acres]

    arcpy.Delete_management(os.path.join(optfolder, "lu_poly.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "lu_soil.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "lu_soil_diss.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "soil_poly.shp"))

    response = {
        'lu_type': lu_desc,
        'soil_list': soil_acre_lists,
        'soil_ac': total_area,
        'cn_ac': acres,
        'cn_per': percent,
        'cn_list': curve_num
    }

    return(json.dumps(response))


def basinstats(request):

    optfolder = str(request.GET['path'])
    landuse = str(request.GET['lu'])
    hyd = str(request.GET['hyd'])
    dem_layer = str(request.GET['dem'])
    soil_layer = str(request.GET['soil'])

    arcpy.env.overwriteOutput = True
    arcpy.env.scratchWorkspace = optfolder
    arcpy.env.workspace = optfolder
    arcpy.env.snapRaster = os.path.join(optfolder, "dem_clip.tif")

    # *******************************************************************************************************
    # Warning messages
    # *******************************************************************************************************
    Impwarntext = """
            IMPERVIOUS AREA IN WATERSHED EXCEEDS 10%.
            Calculated discharges from USGS Regression
            Equations may not be appropriate.
                     """
    provwarntext = """
            Watershed is within 5km of physiographic
            province boundary.  You should consider
            sensitivity of discharges to region location.
                     """
    limewarntext = """
            Watershed is within 1km of underlying limestone
            geology.  You should consider sensitivity
            of discharges to percent limestone calculated.
                     """

    # *******************************************************************************************************
    # Get outlet coordinates and prepare masked grids for calculations
    # *******************************************************************************************************
    out_rast = arcpy.Raster(os.path.join(optfolder, "pour_point.tif"))
    cellsize = out_rast.meanCellWidth
    cellsq = cellsize * cellsize

    basingrid = os.path.join(optfolder, "basingrid.tif")
    dirgrid = arcpy.sa.Times(os.path.join(optfolder, "flowdir.tif"), basingrid)
    elevgrid = arcpy.sa.Times(os.path.join(optfolder, "dem_clip.tif"), basingrid)
    lantype = arcpy.sa.Times(os.path.join(optfolder, "landuse.tif"), basingrid)

    # Get basingrid count [number of pixels]
    shedtab = arcpy.SearchCursor(basingrid, "", "", "Count", "")
    for row in shedtab:
        basinarea = row.getValue("Count")

    # *******************************************************************************************************
    # Compute channel and land slope
    # *******************************************************************************************************
    theslope = hydro.channelslope(dirgrid, elevgrid, optfolder)  # already converted into feet/mile
    theslope_feet = float(theslope / 5280.0)

    maxlength = float(hydro.maxlength)  # already converted into miles for use in channel slope function

    if not os.path.exists(os.path.join(optfolder, "slope_calc")):
        os.mkdir(os.path.join(optfolder, "slope_calc"), 0o755)

    # SLOPE GRID:
    dlgrid_temp1 = arcpy.sa.Log2(dirgrid)
    dlgrid_temp2 = dlgrid_temp1 % 2
    dlgrid_temp3 = arcpy.sa.Con(dlgrid_temp2 > 0, pow(2, 0.5), 1)
    dlgrid_temp4 = arcpy.sa.Times(dlgrid_temp3, cellsize)
    dlgrid = arcpy.sa.Times(dlgrid_temp4, basingrid)
    dir_shift_1 = arcpy.Shift_management(elevgrid, os.path.join(optfolder, "dir_shift_1.tif"), -cellsize, 0)
    dir_shift_2 = arcpy.Shift_management(elevgrid, os.path.join(optfolder, "dir_shift_2.tif"), -cellsize, cellsize)
    dir_shift_4 = arcpy.Shift_management(elevgrid, os.path.join(optfolder, "dir_shift_4.tif"), 0, cellsize)
    dir_shift_8 = arcpy.Shift_management(elevgrid, os.path.join(optfolder, "dir_shift_8.tif"), cellsize, cellsize)
    dir_shift_16 = arcpy.Shift_management(elevgrid, os.path.join(optfolder, "dir_shift_16.tif"), cellsize, 0)
    dir_shift_32 = arcpy.Shift_management(elevgrid, os.path.join(optfolder, "dir_shift_32.tif"), cellsize, -cellsize)
    dir_shift_64 = arcpy.Shift_management(elevgrid, os.path.join(optfolder, "dir_shift_64.tif"), 0, -cellsize)
    dir_shift_128 = arcpy.Shift_management(elevgrid, os.path.join(optfolder, "dir_shift_128.tif"), -cellsize, -cellsize)

    shift_temp1 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_1), dlgrid)
    shift_temp2 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_2), dlgrid)
    shift_temp3 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_4), dlgrid)
    shift_temp4 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_8), dlgrid)
    shift_temp5 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_16), dlgrid)
    shift_temp6 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_32), dlgrid)
    shift_temp7 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_64), dlgrid)
    shift_temp8 = arcpy.sa.Divide(arcpy.sa.Minus(elevgrid, dir_shift_128), dlgrid)

    slope = arcpy.sa.Con(dirgrid == 1, shift_temp1,
                         arcpy.sa.Con(dirgrid == 2, shift_temp2,
                                      arcpy.sa.Con(dirgrid == 4, shift_temp3,
                                                   arcpy.sa.Con(dirgrid == 8, shift_temp4,
                                                                arcpy.sa.Con(dirgrid == 16, shift_temp5,
                                                                             arcpy.sa.Con(dirgrid == 32,
                                                                                          shift_temp6,
                                                                                          arcpy.sa.Con(
                                                                                              dirgrid == 64,
                                                                                              shift_temp7,
                                                                                              arcpy.sa.Con(
                                                                                                  dirgrid == 128,
                                                                                                  shift_temp8,
                                                                                                  0))))))))
    slope.save(os.path.join(optfolder, "slope_calc/landslope.tif"))
    slope = arcpy.Raster(os.path.join(optfolder, "slope_calc/landslope.tif"))

    # Slope value was going below 0 which isn't correct (raw dem used instead of filled). Condition is added to at least have 0.01
    slopegrid = arcpy.sa.Con(slope > 0, slope, 0.001)
    landsloperesult = arcpy.GetRasterProperties_management(slopegrid, "MEAN")
    landslopevalue = float(landsloperesult.getOutput(0))
    landslope = float(landslopevalue) / 3.28084  # modified: 09/06/2018 (divided by 3.28084)

    # *******************************************************************************************************
    # Determine province of outlet location
    # *******************************************************************************************************
    prov = os.path.join(directory, "shp/md_prov/Mdprov.shp")
    arcpy.Intersect_analysis([prov, os.path.join(optfolder, "pour_point.shp")], os.path.join(optfolder, "prov_inter.shp"))

    ProvTab = arcpy.SearchCursor(os.path.join(optfolder, "prov_inter.shp"), "", "", "PROVINCE", "")
    for row in ProvTab:
        if row.getValue("PROVINCE").upper() == "A":
            provstring = "Appalachian Plateaus and Allegheny Ridges"
        elif row.getValue("PROVINCE").upper() == "B":
            provstring = "Blue Ridge and Great Valley"
        elif row.getValue("PROVINCE").upper() == "P":
            provstring = "Piedmont"
        elif row.getValue("PROVINCE").upper() == "W":
            provstring = "Western Coastal Plain"
        elif row.getValue("PROVINCE").upper() == "E":
            provstring = "Eastern Coastal Plain"
        else:
            row.getValue("PROVINCE") == "Unknown"
            provstring == "Unknown"

    # *******************************************************************************************************
    # Get percent soil types from soil dataset
    # *******************************************************************************************************
    ssurgo = arcpy.sa.Times(os.path.join(optfolder, "soils.tif"), basingrid)
    pct = hydro.SSURGOPct(basinarea, ssurgo)
    pctSoil = list(map(float, pct))
    pctAsoil = float(pctSoil[0])
    pctBsoil = float(pctSoil[1])
    pctCsoil = float(pctSoil[2])
    pctDsoil = float(pctSoil[3])
    pctWsoil = float(pctSoil[4])

    # *******************************************************************************************************
    # Get LU count for Urban, Nil, Forest, and Storage -- it will be different for
    # MOP, MD/DE, MRLC, and USGS
    # *******************************************************************************************************
    # *******************************************************************************************************
    # Get Impervious count -- count will vary depending upon choice of input Landuse and Hyd condition
    # *******************************************************************************************************

    if landuse == "nlcd2011":
        count = hydro.GetLUCountNLCD(basinarea, lantype)
        LUcount = list(map(float, count))
        UrbPct = float(LUcount[0])
        FC = float(LUcount[2])
        ST = float(LUcount[3])
        if hyd == "Fair":
            Impcount = hydro.GetImpCountNLCDFair(lantype)
        elif hyd == "Good":
            Impcount = hydro.GetImpCountNLCDGood(lantype)
        elif hyd == "Poor":
            Impcount = hydro.GetImpCountNLCDPoor(lantype)
    elif landuse == "nlcd2006":
        count = hydro.GetLUCountNLCD(basinarea, lantype)
        LUcount = list(map(float, count))
        UrbPct = float(LUcount[0])
        FC = float(LUcount[2])
        ST = float(LUcount[3])
        if hyd == "Fair":
            Impcount = hydro.GetImpCountNLCDFair(lantype)
        elif hyd == "Good":
            Impcount = hydro.GetImpCountNLCDGood(lantype)
        elif hyd == "Poor":
            Impcount = hydro.GetImpCountNLCDPoor(lantype)

    elif landuse == "nlcd2001":
        count = hydro.GetLUCountNLCD(basinarea, lantype)
        LUcount = list(map(float, count))
        UrbPct = float(LUcount[0])
        FC = float(LUcount[2])
        ST = float(LUcount[3])
        if hyd == "Fair":
            Impcount = hydro.GetImpCountNLCDFair(lantype)
        elif hyd == "Good":
            Impcount = hydro.GetImpCountNLCDGood(lantype)
        elif hyd == "Poor":
            Impcount = hydro.GetImpCountNLCDPoor(lantype)
    elif landuse == "lu2010":
        count = hydro.GetLUCountAnderson(basinarea, lantype)
        LUcount = list(map(float, count))
        UrbPct = float(LUcount[0])
        FC = float(LUcount[2])
        ST = float(LUcount[3])
        if hyd == "Fair":
            Impcount = hydro.GetImpCountAndersonFair(lantype)
        elif hyd == "Good":
            Impcount = hydro.GetImpCountAndersonGood(lantype)
        elif hyd == "Poor":
            Impcount = hydro.GetImpCountAndersonPoor(lantype)
    elif landuse == "mdplu2002":
        count = hydro.GetLUCountAnderson(basinarea, lantype)
        LUcount = list(map(float, count))
        UrbPct = float(LUcount[0])
        FC = float(LUcount[2])
        ST = float(LUcount[3])
        if hyd == "Fair":
            Impcount = hydro.GetImpCountAndersonFair(lantype)
        elif hyd == "Good":
            Impcount = hydro.GetImpCountAndersonGood(lantype)
        elif hyd == "Poor":
            Impcount = hydro.GetImpCountAndersonPoor(lantype)
    elif landuse == "lu97m":
        count = hydro.GetLUCountAnderson(basinarea, lantype)
        LUcount = list(map(float, count))
        UrbPct = float(LUcount[0])
        FC = float(LUcount[2])
        ST = float(LUcount[3])
        if hyd == "Fair":
            Impcount = hydro.GetImpCountAndersonFair(lantype)
        elif hyd == "Good":
            Impcount = hydro.GetImpCountAndersonGood(lantype)
        elif hyd == "Poor":
            Impcount = hydro.GetImpCountAndersonPoor(lantype)
    elif landuse == "mdde2002":
        count = hydro.GetLUCountMDDE(basinarea, lantype)
        LUcount = list(map(float, count))
        UrbPct = float(LUcount[0])
        FC = float(LUcount[2])
        ST = float(LUcount[3])
        if hyd == "Fair":
            Impcount = hydro.GetImpCountMDDEFair(lantype)
        elif hyd == "Good":
            Impcount = hydro.GetImpCountMDDEGood(lantype)
        elif hyd == "Poor":
            Impcount = hydro.GetImpCountMDDEPoor(lantype)
    elif landuse == "luult":
        count = hydro.GetLUCountUltimate(basinarea, lantype)
        LUcount = list(map(float, count))
        UrbPct = float(LUcount[0])
        FC = float(LUcount[2])
        ST = float(LUcount[3])
        if hyd == "Fair":
            Impcount = hydro.GetImpCountUltimateFair(lantype)
        elif hyd == "Good":
            Impcount = hydro.GetImpCountUltimateGood(lantype)
        elif hyd == "Poor":
            Impcount = hydro.GetImpCountUltimatePoor(lantype)
    elif landuse == "mrlc":
        count = hydro.GetLUCountMRLC(basinarea, lantype)
        LUcount = list(map(float, count))
        UrbPct = float(LUcount[0])
        FC = float(LUcount[2])
        ST = float(LUcount[3])
        if hyd == "Fair":
            Impcount = hydro.GetImpCountMRLCFair(lantype)
        elif hyd == "Good":
            Impcount = hydro.GetImpCountMRLCGood(lantype)
        elif hyd == "Poor":
            Impcount = hydro.GetImpCountMRLCPoor(lantype)
    elif landuse == "lu70":
        count = hydro.GetLUCountUSGS(basinarea, lantype)
        LUcount = list(map(float, count))
        UrbPct = float(LUcount[0])
        FC = float(LUcount[2])
        ST = float(LUcount[3])
        if hyd == "Fair":
            Impcount = hydro.GetImpCountUSGSFair(lantype)
        elif hyd == "Good":
            Impcount = hydro.GetImpCountUSGSGood(lantype)
        elif hyd == "Poor":
            Impcount = hydro.GetImpCountUSGSPoor(lantype)

    IA = float((Impcount / basinarea) * 100)

    # *******************************************************************************************************
    # Get Limestone percent count
    # *******************************************************************************************************

    limestonem = os.path.join(directory, "shp/md_limestone/limestonem.shp")
    LIcnt = 0
    try:
        limegrid = arcpy.sa.ExtractByMask(basingrid, limestonem)
        limegrid.save(os.path.join(optfolder, "limegrid.tif"))
        arcpy.BuildRasterAttributeTable_management(os.path.join(optfolder, "limegrid.tif"), "Overwrite")
        with arcpy.da.SearchCursor(os.path.join(optfolder, "limegrid.tif"), "Count") as rows:
            for row in rows:
                LIcnt += row[0] or 0
        LI = float((float(LIcnt) / basinarea) * 100)  # 10-23-2013
    except:
        LI = 0

    areami2 = float((basinarea * cellsq) / 2588881)  # conversion into sq miles

    # *******************************************************************************************************
    # Get basi relief [it is difference of mean elevation and outlet elevation]
    # *******************************************************************************************************
    elev1 = arcpy.GetRasterProperties_management(elevgrid, "MEAN")
    mean_elev = float(elev1.getOutput(0))

    for row in arcpy.da.SearchCursor(os.path.join(optfolder, "pour_point.shp"), ["SHAPE@XY"]):
        x, y = row[0]
    outlet_elev = arcpy.management.GetCellValue(os.path.join(optfolder, "dem_clip.tif"), "{} {}".format(x, y))
    outletelev = float(outlet_elev.getOutput(0))
    basinrelief = float(mean_elev - outletelev)  # Assuming it is already converted into feets

    # *******************************************************************************************************
    # Average CN number using above outgrid depending upon user choice of HydCon
    # *******************************************************************************************************
    cnGrid = arcpy.sa.Times(os.path.join(optfolder, "curvenumber.tif"), basingrid)
    avgCN_val = arcpy.GetRasterProperties_management(cnGrid, "MEAN")  # figure out outgrid source
    avgCN = float(avgCN_val.getOutput(0))

    ssurgan = arcpy.sa.Times(os.path.join(optfolder, "soils.tif"), basingrid)
    pctR = hydro.SoilPct(basinarea, ssurgan)
    pctSoil = list(map(float, pctR))
    pctAR = float(pctSoil[0])
    pctAsoilR = float((pctAR / basinarea) * 100)
    pctBR = float(pctSoil[1])
    pctBsoilR = float((pctBR / basinarea) * 100)
    pctCR = float(pctSoil[2])
    pctCsoilR = float((pctCR / basinarea) * 100)
    pctDR = float(pctSoil[3])
    pctDsoilR = float((pctDR / basinarea) * 100)

    """
    The following code calculates the Time of Concentration. If multiple provinces
    are involved, tc is weighted average of area of watershed in each province.
    More correct would be to perform weighted average based on length of channel
    in each province.  This modification will be performed at a later time.  (GEM - 12/01/99)
    """

    theVTab = os.path.join(optfolder, "theVTab.dbf")

    # *******************************************************************************************************
    # don"t add "theVTab" to TOC -- it will change list by drawing order to list
    # by source which will prohibit addition of new layers
    # *******************************************************************************************************
    arcpy.sa.ZonalStatisticsAsTable(prov, "PROVINCE", basingrid, "theVTab", "DATA", "ALL")
    arcpy.DeleteField_management("theVTab", "ZONE_CODE;MIN;MAX;RANGE;MEAN;STD;SUM;VARIETY;MAJORITY;MINORITY;MEDIAN")
    addFieldNameList = ["Q1.25", "Q1.50", "Q1.75", "Q2", "Q5", "Q10", "Q25", "Q50", "Q100", "Q200", "Q500"]
    for each in addFieldNameList:
        if not len(arcpy.ListFields("theVTab", each)) > 0:
            arcpy.AddField_management("theVTab", each, "FLOAT", 10, 3)

    # *******************************************************************************************************
    # Create regionlist and regionarea and declare them as global for use in
    #  Thomas Discharge script (is it used in Thomas Discharge script?)
    # *******************************************************************************************************
    sumarea = 0
    theVTab = arcpy.SearchCursor("theVTab", "", "", "Count", "")
    for each in theVTab:
        count = each.getValue("Count")
        sumarea = sumarea + count
    sumArea = sumarea
    del each

    regionlist = []
    regionarea = []
    breakstring = [""]
    theVTab = arcpy.SearchCursor("theVTab", "", "", "Province;Count", "")
    for row in theVTab:
        theProv = row.getValue("Province")
        theArea = float(row.getValue("Count"))
        areapercent = float((theArea / sumArea) * 100)
        areapercent = "{0:.2f}".format(areapercent)
        regionlist.append(theProv)
        regionarea.append(areapercent)
        if row.getValue("Province") == "A":
            breakstring.append("Appalachian Plateaus and Allegheny Ridges %s percent of area" % (areapercent))
        elif row.getValue("Province") == "B":
            breakstring.append("Blue Ridge and Great Valley %s percent of area" % (areapercent))
        elif row.getValue("Province") == "P":
            breakstring.append("Piedmont %s percent of area" % (areapercent))
        elif row.getValue("Province") == "W":
            breakstring.append("Western Coastal Plain %s percent of area" % (areapercent))
        elif row.getValue("Province") == "E":
            breakstring.append("Eastern Coastal Plain %s percent of area"  % (areapercent))


    # *******************************************************************************************************
    # Compute Time of Concentration:
    #                               1]  W.O. Thomas, Jr. Equation   [tc]
    #                               2]  SCS Lag equation * 1.67     [lagtime]
    # *******************************************************************************************************
    sumtc = 0
    theVTab = arcpy.SearchCursor("theVTab", "", "", "Province;Count", "")
    for row in theVTab:
        theProv = row.getValue("Province")
        theArea = row.getValue("Count")

        if row.getValue("Province") == "A":
            temptc = 0.133 * ((maxlength) ** (0.475)) * ((theslope) ** (-0.187)) * ((101 - FC) ** (-0.144)) * (
                    (101 - IA) ** (0.861)) * ((ST + 1) ** (0.154)) * ((10) ** (0.194))
        elif row.getValue("Province") == "W":
            temptc = 0.133 * ((maxlength) ** (0.475)) * ((theslope) ** (-0.187)) * ((101 - FC) ** (-0.144)) * (
                    (101 - IA) ** (0.861)) * ((ST + 1) ** (0.154)) * ((10) ** (0.366))
        elif row.getValue("Province") == "E":
            temptc = 0.133 * ((maxlength) ** (0.475)) * ((theslope) ** (-0.187)) * ((101 - FC) ** (-0.144)) * (
                    (101 - IA) ** (0.861)) * ((ST + 1) ** (0.154)) * ((10) ** (0.366))
        else:
            temptc = 0.133 * ((maxlength) ** (0.475)) * ((theslope) ** (-0.187)) * ((101 - FC) ** (-0.144)) * (
                    (101 - IA) ** (0.861)) * ((ST + 1) ** (0.154))
        sumtc = sumtc + (temptc * theArea)
    del row

    tc = (sumtc / basinarea)

    # *******************************************************************************************************
    # Calculate lagtime
    # *******************************************************************************************************

    lagtime = ((float(100 * ((maxlength * 5280) ** (0.8)) * (((1000 / avgCN) - 9) ** (0.7))) / float(
            1900 * ((abs(landslope) * 100) ** (0.5)))) / 60)

    # *******************************************************************************************************
    # Calculate Mean Annual Precipitation
    # *******************************************************************************************************
    mapstpm = os.path.join(directory, "rasters/map/mapstpm")
    prec_grid = os.path.join(directory, "prec/p2-24m")
    arcpy.env.cellSize = str(cellsize)
    basingrid_p = arcpy.Raster(basingrid)

    maprecbasin = arcpy.sa.Times(mapstpm, basingrid_p)  # Make sure basingrid has value 1 otherwise all precip will be 0
    theprec = arcpy.sa.Times(prec_grid, basingrid_p)
    precavg = arcpy.GetRasterProperties_management(maprecbasin, "MEAN")
    precavg = float(precavg.getOutput(0))
    avgprec = arcpy.GetRasterProperties_management(theprec, "MEAN")
    avgprec = float(avgprec.getOutput(0))
    maprec = float(precavg / (1000 * 2.54))
    p2yr = float(avgprec / 1000)

    # *******************************************************************************************************
    # format precision before text file string settings
    # *******************************************************************************************************
    maxlength = "{0:.2f}".format(maxlength)
    theslope = "{0:.8f}".format(theslope)
    theslope_feet = "{0:.8f}".format(theslope_feet)
    landslope = "{0:.8f}".format(landslope)
    pctAsoil = "{0:.1f}".format(pctAsoil)
    pctBsoil = "{0:.1f}".format(pctBsoil)
    pctCsoil = "{0:.1f}".format(pctCsoil)
    pctDsoil = "{0:.1f}".format(pctDsoil)
    pctWsoil = "{0:.1f}".format(pctWsoil)
    UrbPct = "{0:.1f}".format(UrbPct)
    FC = "{0:.1f}".format(FC)
    ST = "{0:.1f}".format(ST)
    IA = "{0:.1f}".format(IA)
    LI = "{0:.1f}".format(LI)
    areami2 = "{0:.2f}".format(areami2)
    basinrelief = "{0:.2f}".format(basinrelief)
    avgCN = "{0:.1f}".format(avgCN)
    pctAsoilR = "{0:.1f}".format(pctAsoilR)
    pctBsoilR = "{0:.1f}".format(pctBsoilR)
    pctCsoilR = "{0:.1f}".format(pctCsoilR)
    pctDsoilR = "{0:.1f}".format(pctDsoilR)
    tc = "{0:.2f}".format(tc)
    lagtime = "{0:.2f}".format(lagtime)
    maprec = "{0:.2f}".format(maprec)
    p2yr = "{0:.2f}".format(p2yr)


    # *******************************************************************************************************
    # Print out Impervious area warning message
    # *** warning message is included in for loop despite the fact that technically it could be printed
    #     twice. Since both "Appalachian Plateau" and "Eastern Coastal Plain" are far apart so it is
    #     impossible to have that big watershed while doing analysis with GISHydroNXT
    # *******************************************************************************************************
    html_warning = ""
    theVTab = arcpy.SearchCursor("theVTab", "", "", "Province", "")
    for row in theVTab:
        if (row.getValue("Province") == "A") or (row.getValue("Province") == "E"):
            if float(IA) >= 10:
                html_warning = html_warning + Impwarntext + "\n"

    # *******************************************************************************************************
    # Close to boundary condition for provinces -- Near tool isn"t available with
    # basic level license therefore a more crude method was emplyed here. It can
    # be improved in future by using "arcpy.Geometry()" tool to get distance
    # *******************************************************************************************************
    arcpy.Buffer_analysis(os.path.join(optfolder, "wshed.shp"), os.path.join(optfolder, "wats_prov.shp"), "5000", "#", "#", "ALL", "FID")
    province = os.path.join(directory, "shp/md_provlines/provlines.shp")
    arcpy.Intersect_analysis([province, os.path.join(optfolder, "wats_prov.shp")], os.path.join(optfolder, "prov_int.shp"))
    prov_cursor = arcpy.SearchCursor(os.path.join(optfolder, "prov_int.shp"), "", "", "FID", "")
    prov = prov_cursor.next()
    if prov != None:
        html_warning = html_warning + provwarntext + "\n"

    # *******************************************************************************************************
    # Close to boundary condition for limestone -- Near tool isn"t available with
    # basic level license therefore a more crude method was emplyed here. It can
    # be improved in future by using "arcpy.Geometry()" tool to get distance
    # *******************************************************************************************************
    arcpy.Buffer_analysis(os.path.join(optfolder, "wshed.shp"), os.path.join(optfolder, "wats_lime.shp"), "1000", "#", "#", "ALL", "FID")
    wats_lime = os.path.join(optfolder, "wats_lime.shp")
    lime_int = os.path.join(optfolder, "lime_int.shp")
    arcpy.Intersect_analysis([limestonem, wats_lime], lime_int, "ALL", "#", "INPUT")
    lime_cursor = arcpy.SearchCursor(lime_int, "", "", "FID", "")
    lime = lime_cursor.next()
    if lime != None:
        html_warning = html_warning + limewarntext + "\n"

    arcpy.Delete_management(os.path.join(optfolder, "dir_shift_1.tif"))
    arcpy.Delete_management(os.path.join(optfolder, "dir_shift_2.tif"))
    arcpy.Delete_management(os.path.join(optfolder, "dir_shift_4.tif"))
    arcpy.Delete_management(os.path.join(optfolder, "dir_shift_8.tif"))
    arcpy.Delete_management(os.path.join(optfolder, "dir_shift_16.tif"))
    arcpy.Delete_management(os.path.join(optfolder, "dir_shift_32.tif"))
    arcpy.Delete_management(os.path.join(optfolder, "dir_shift_64.tif"))
    arcpy.Delete_management(os.path.join(optfolder, "dir_shift_128.tif"))
    arcpy.Delete_management(os.path.join(optfolder, "lime_int.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "prov_int.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "prov_inter.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "wats_lime.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "wats_prov.shp"))

    response = {
        'html_warning': html_warning,
        'dem_layer': dem_layer,
        'landuse': landuse,
        'soil_layer': soil_layer,
        'hyd': hyd,
        'x': int(x),
        'y': int(y),
        'provstring': provstring,
        'areami2': areami2,
        'breakstring': breakstring,
        'theslope': theslope,
        'theslope_feet': theslope_feet,
        'landslope': landslope,
        'UrbPct': UrbPct,
        'IA': IA,
        'tc': tc,
        'lagtime': lagtime,
        'maxlength': maxlength,
        'basinrelief': basinrelief,
        'avgCN': avgCN,
        'FC': FC,
        'ST': ST,
        'LI': LI,
        'pctAsoilR': pctAsoilR,
        'pctBsoilR': pctBsoilR,
        'pctCsoilR': pctCsoilR,
        'pctCsoilR': pctCsoilR,
        'pctDsoilR': pctDsoilR,
        'pctAsoil': pctAsoil,
        'pctBsoil': pctBsoil,
        'pctCsoil': pctCsoil,
        'pctDsoil': pctDsoil,
        'p2yr': p2yr,
        'maprec': maprec,
    }

    return(json.dumps(response))

def discharge(request):

    optfolder = str(request.GET['path'])
    areami2 = float(request.GET['area'])
    landslope = float(request.GET['ls'])
    IA = float(request.GET['ia'])
    LI = float(request.GET['li'])
    gageid = str(request.GET['gageid'])

    arcpy.env.overwriteOutput = True
    arcpy.env.scratchWorkspace = optfolder
    arcpy.env.workspace = optfolder
    arcpy.env.snapRaster = os.path.join(optfolder, "dem_clip.tif")


    # ******************************************************************************************************
    # read global variables (defined in basin stat) into variable names for use
    # in Thomas peak discharge analysis
    # ******************************************************************************************************
    FC = float(hydro.FC)  # LI is already declared as global variable
    FC = "{0:.2f}".format(FC)
    DA = areami2
    HA = float(hydro.pctAsoil)
    HC = float(hydro.pctCsoil)
    HD = float(hydro.pctDsoil)
    SLL = float(landslope)
    HCD = float(HC + HD)
    ImpA = float(IA)

    # ******************************************************************************************************
    # define initial values and lists. add and index those lists as part of dictionary
    # ******************************************************************************************************
    sumarea = 0
    provstring = ""
    sQ1p25 = 0
    sQ1p50 = 0
    sQ2 = 0
    sQ5 = 0
    sQ10 = 0
    sQ25 = 0
    sQ50 = 0
    sQ100 = 0
    sQ200 = 0
    sQ500 = 0
    Q1p25list = [0, 0, 0, 0, 0, 0, 0, 0]
    Q1p50list = [0, 0, 0, 0, 0, 0, 0, 0]
    Q2list = [0, 0, 0, 0, 0, 0, 0, 0]
    Q5list = [0, 0, 0, 0, 0, 0, 0, 0]
    Q10list = [0, 0, 0, 0, 0, 0, 0, 0]
    Q25list = [0, 0, 0, 0, 0, 0, 0, 0]
    Q50list = [0, 0, 0, 0, 0, 0, 0, 0]
    Q100list = [0, 0, 0, 0, 0, 0, 0, 0]
    Q200list = [0, 0, 0, 0, 0, 0, 0, 0]
    Q500list = [0, 0, 0, 0, 0, 0, 0, 0]

    qlist = {"1": [], "2": [], "3": [], "4": [], "5": [], "6": [], "7": [], "8": [], "9": [], "10": []}
    qlist["1"].extend(Q1p25list)
    qlist["2"].extend(Q1p50list)
    qlist["3"].extend(Q2list)
    qlist["4"].extend(Q5list)
    qlist["5"].extend(Q10list)
    qlist["6"].extend(Q25list)
    qlist["7"].extend(Q50list)
    qlist["8"].extend(Q100list)
    qlist["9"].extend(Q200list)
    qlist["10"].extend(Q500list)

    # ******************************************************************************************************
    # read zonal stat table (theVTab), declare fields into variables, and
    # loop through theVTab fields to getValue
    # ******************************************************************************************************
    theVTab = os.path.join(optfolder, "theVTab.dbf")
    theVTab = arcpy.SearchCursor("theVTab", "", "", "Count", "")

    # ******************************************************************************************************
    # loop to get total count of pixels of basingrid
    # ******************************************************************************************************
    for each in theVTab:
        count = each.getValue("Count")
        sumarea = sumarea + count
    sumArea = sumarea

    theVTab = os.path.join(optfolder, "theVTab.dbf")
    theVTab = arcpy.SearchCursor("theVTab", "", "", "Province;Count;Q1.25;Q1.50;Q1.75;Q2;Q5;Q10;Q25;Q50;Q100;Q200;Q500", "")

    # loop to begin tasker handling and area weighted analysis
    discharge_values = []
    flood_intervals = []
    for row in theVTab:
        AreaField = float(row.getValue("Count"))
        areapercent = float((AreaField / sumArea) * 100)
        if row.getValue("Province") == "A":
            provstring = "Appalachian Plateaus and Allegheny Ridges %s percent of area" "\n" % ("{0:.2f}".format(areapercent))
        elif row.getValue("Province") == "B":
            provstring = "Blue Ridge and Great Valley %s percent of area" "\n" % ("{0:.2f}".format(areapercent))
        elif row.getValue("Province") == "P":
            provstring = "Piedmont %s percent of area" "\n" % ("{0:.2f}".format(areapercent))
        elif row.getValue("Province") == "W":
            provstring = "Western Coastal Plain %s percent of area" "\n" % ("{0:.2f}".format(areapercent))
        elif row.getValue("Province") == "E":
            provstring = "Eastern Coastal Plain %s percent of area" "\n" % ("{0:.2f}".format(areapercent))
        else:
            provstring = "No Province Selected"

        intasker = []
        if row.getValue("Province") == "A":
            intasker.append("A")
            intasker.append("%s" % (DA))
            intasker.append("%s" % (SLL))
        elif row.getValue("Province") == "B":
            intasker.append("p" )
            intasker.append("%s" % (DA))
            intasker.append("%s" % (FC))
            intasker.append("%s" % (LI))
            intasker.append("%s" % (ImpA))
        elif row.getValue("Province") == "P":
            intasker.append("P" )
            intasker.append("%s" % (DA))
            intasker.append("%s" % (FC))
            intasker.append("%s" % (LI))
            intasker.append("%s" % (ImpA))

        elif row.getValue("Province") == "W":
            intasker.append("WC")
            intasker.append("%s" % (DA))
            intasker.append("%s" % (ImpA))
            intasker.append("%s" % (HCD))
        elif row.getValue("Province") == "E":
            intasker.append("EC")
            intasker.append("%s" % (DA))
            intasker.append("%s" % (SLL))
            intasker.append("%s" % (HA))
        intasker.append("%s" % (gageid))

        tasker_response = tasker.RRE(directory, optfolder, intasker)

        yhat_Q = tasker_response['yhat_list']

        Q1p25 = yhat_Q[0]
        discharge_values.append(Q1p25)
        Q1p50 = yhat_Q[1]
        discharge_values.append(Q1p50)
        Q2 = yhat_Q[2]
        discharge_values.append(Q2)
        Q5 = yhat_Q[3]
        discharge_values.append(Q5)
        Q10 = yhat_Q[4]
        discharge_values.append(Q10)
        Q25 = yhat_Q[5]
        discharge_values.append(Q25)
        Q50 = yhat_Q[6]
        discharge_values.append(Q50)
        Q100 = yhat_Q[7]
        discharge_values.append(Q100)
        Q200 = yhat_Q[8]
        discharge_values.append(Q200)
        Q500 = yhat_Q[9]
        discharge_values.append(Q500)

        # *****************************************************************************
        # compute discharge and assign confidence intervals to qlist entries
        # *****************************************************************************
        Q1p25 = float(Q1p25)
        sQ1p25 = sQ1p25 + (Q1p25 * AreaField)
        Q1p50 = float(Q1p50)
        sQ1p50 = sQ1p50 + (Q1p50 * AreaField)
        Q2 = float(Q2)
        sQ2 = sQ2 + (Q2 * AreaField)
        Q5 = float(Q5)
        sQ5 = sQ5 + (Q5 * AreaField)
        Q10 = float(Q10)
        sQ10 = sQ10 + (Q10 * AreaField)
        Q25 = float(Q25)
        sQ25 = sQ25 + (Q25 * AreaField)
        Q50 = float(Q50)
        sQ50 = sQ50 + (Q50 * AreaField)
        Q100 = float(Q100)
        sQ100 = sQ100 + (Q100 * AreaField)
        Q200 = float(Q200)
        sQ200 = sQ200 + (Q200 * AreaField)
        Q500 = float(Q500)
        sQ500 = sQ500 + (Q500 * AreaField)

        cl = tasker_response['cl']
        cu = tasker_response['cu']
        for i in range(10):
            c1 = [(float(cl[i][0]) * areapercent) / 100]
            c2 = [(float(cu[i][0]) * areapercent) / 100]
            c3 = [(float(cl[i][1]) * areapercent) / 100]
            c4 = [(float(cu[i][1]) * areapercent) / 100]
            c5 = [(float(cl[i][2]) * areapercent) / 100]
            c6 = [(float(cu[i][2]) * areapercent) / 100]
            c7 = [(float(cl[i][3]) * areapercent) / 100]
            c8 = [(float(cu[i][3]) * areapercent) / 100]
            qlist = [c1, c2, c3, c4, c5, c6, c7, c8]
            flood_intervals.append(qlist)

    lists = {i: [el[0] for el in v] for i, v in enumerate(flood_intervals, start=1)}

    # ******************************************************************************************************
    # index and count number of sub-lists -- prepare for TaskerString function
    # ******************************************************************************************************
    number = hydro.FloodIntervalLists(lists)
    if number > 10:
        q_list1 = [sum(i) for i in zip(lists[1], lists[11])]
        q_list2 = [sum(i) for i in zip(lists[2], lists[12])]
        q_list3 = [sum(i) for i in zip(lists[3], lists[13])]
        q_list4 = [sum(i) for i in zip(lists[4], lists[14])]
        q_list5 = [sum(i) for i in zip(lists[5], lists[15])]
        q_list6 = [sum(i) for i in zip(lists[6], lists[16])]
        q_list7 = [sum(i) for i in zip(lists[7], lists[17])]
        q_list8 = [sum(i) for i in zip(lists[8], lists[18])]
        q_list9 = [sum(i) for i in zip(lists[9], lists[19])]
        q_list10 = [sum(i) for i in zip(lists[10], lists[20])]
    else:
        q_list1 = lists[1]
        q_list2 = lists[2]
        q_list3 = lists[3]
        q_list4 = lists[4]
        q_list5 = lists[5]
        q_list6 = lists[6]
        q_list7 = lists[7]
        q_list8 = lists[8]
        q_list9 = lists[9]
        q_list10 = lists[10]
    q_list_all = [q_list1, q_list2, q_list3, q_list4, q_list5, q_list6, q_list7, q_list8, q_list9, q_list10]
    it_values = ['1.25', '1.50', '2', '5', '10', '25', '50', '100', '200', '500']

    # ******************************************************************************************************
    # discharge computation based on province
    # ******************************************************************************************************
    Q1p25 = int(sQ1p25 / sumArea)
    Q1p50 = int(sQ1p50 / sumArea)
    Q2 = int(sQ2 / sumArea)
    Q5 = int(sQ5 / sumArea)
    Q10 = int(sQ10 / sumArea)
    Q25 = int(sQ25 / sumArea)
    Q50 = int(sQ50 / sumArea)
    Q100 = int(sQ100 / sumArea)
    Q200 = int(sQ200 / sumArea)
    Q500 = int(sQ500 / sumArea)
    Qcfs = [Q1p25,Q1p50,Q2,Q5,Q10,Q25,Q50,Q100,Q200,Q500]

    response = {
        'it_values': it_values,
        'q_list_all': q_list_all,
        'provstring': provstring,
        'Qcfs': Qcfs,
    }
    response.update(tasker_response)

    return(json.dumps(response))

def flowpath(request):

    optfolder = str(request.GET['path'])
    x = float(request.GET['x'])
    y = float(request.GET['y'])
    i = int(request.GET['i'])

    arcpy.env.overwriteOutput = True
    arcpy.env.scratchWorkspace = optfolder
    arcpy.env.workspace = optfolder
    arcpy.env.snapRaster = os.path.join(optfolder, "dem_clip.tif")

    try:
        arcpy.Delete_management(os.path.join(optfolder, "line" + str(i-1) +".shp"))
        arcpy.Delete_management(os.path.join(optfolder, "line" + str(i+1) +".shp"))
        arcpy.Delete_management(os.path.join(optfolder, "line" + str(i) +".shp"))
    except:
        pass

    xmd, ymd = coord_transf(x, y, sr_map, sr_md)

    point = arcpy.Point(xmd, ymd)
    ptGeometry = arcpy.PointGeometry(point)
    theLine = arcpy.sa.CostPath(ptGeometry, os.path.join(optfolder, "dem_clip.tif"), os.path.join(optfolder, "flowdir.tif"), "BEST_SINGLE")
    theLine_masked = arcpy.sa.Times(theLine, os.path.join(optfolder, "basingrid.tif"))
    arcpy.sa.StreamToFeature(theLine_masked, os.path.join(optfolder, "flowdir.tif"), os.path.join(optfolder, "line" + str(i) +".shp"), "NO_SIMPLIFY")

    searchcursor = arcpy.da.SearchCursor(os.path.join(optfolder, "line" + str(i) +".shp"), ("OID@", "SHAPE@"))
    firstrow = searchcursor.next()
    insertcursor = arcpy.da.InsertCursor(os.path.join(optfolder, "addasstreams.shp"), ("OID@", "SHAPE@"))
    insertcursor.insertRow(firstrow)

    return

def clearflowpath(request):
    optfolder = str(request.GET['path'])
    i = int(request.GET['i'])
    arcpy.DeleteRows_management(os.path.join(optfolder, "addasstreams.shp"))
    arcpy.Delete_management(os.path.join(optfolder, "line" + str(i) +".shp"))
    return

def outlet(request):

    optfolder = str(request.GET['path'])
    x = float(request.GET['x'])
    y = float(request.GET['y'])

    arcpy.env.overwriteOutput = True
    arcpy.env.scratchWorkspace = optfolder
    arcpy.env.workspace = optfolder
    arcpy.env.snapRaster = os.path.join(optfolder, "dem_clip.tif")

    xmd, ymd = coord_transf(x, y, sr_map, sr_md)
    xy = (xmd, ymd)

    outletxy = arcpy.sa.ExtractByPoints(os.path.join(optfolder, "infstreams.tif"), [arcpy.Point(xmd, ymd)], "INSIDE")
    outletxy.save(os.path.join(optfolder, "outletxy.tif"))

    aux_max = arcpy.sa.Raster(os.path.join(optfolder, "outletxy.tif")).maximum
    arcpy.Delete_management(os.path.join(optfolder, "outletxy.tif"), "")

    if aux_max is None:
        return(False)

    cursor = arcpy.da.InsertCursor(os.path.join(optfolder, "addasoutlets.shp"), ("SHAPE@XY"))
    cursor.insertRow([xy])
    return(True)

def clearoutlet(request):
    optfolder = str(request.GET['path'])
    arcpy.DeleteRows_management(os.path.join(optfolder, "addasoutlets.shp"))
    return

def subsheds(request):

    optfolder = str(request.GET['path'])
    x = float(request.GET['x'])
    y = float(request.GET['y'])

    arcpy.env.overwriteOutput = True
    arcpy.env.scratchWorkspace = optfolder
    arcpy.env.workspace = optfolder
    arcpy.env.snapRaster = os.path.join(optfolder, "dem_clip.tif")

    arcpy.PointToRaster_conversion(os.path.join(optfolder, "addasoutlets.shp"), "FID", os.path.join(optfolder, "Outlets_temp"), "MOST_FREQUENT", "NONE", 30)
    outlets_adj = arcpy.sa.Plus(os.path.join(optfolder, "Outlets_temp"), 1)
    outlets_adj.save(os.path.join(optfolder, "AddOutlets"))
    outlets_custom = arcpy.sa.Con(arcpy.sa.IsNull(os.path.join(optfolder, "AddOutlets")), 0, os.path.join(optfolder, "AddOutlets"))
    outlets_custom = arcpy.sa.SetNull(outlets_custom, outlets_custom, "VALUE = 0")
    outlets_custom.save(os.path.join(optfolder, "outlets_user"))

    arcpy.env.addOutputsToMap = True
    arcpy.env.extent = "MAXOF"
    arcpy.Merge_management(os.path.join(optfolder, "AddasStreams.shp"), os.path.join(optfolder, "StrmMerge.shp"))  # no apparent benefit of using merge

    arcpy.FeatureToRaster_conversion(os.path.join(optfolder, "StrmMerge.shp"), "Id", os.path.join(optfolder, "ModStr"), "#")
    ModStr = arcpy.sa.Times(os.path.join(optfolder, "ModStr"), os.path.join(optfolder, "basingrid"))
    Streams = arcpy.sa.Con(ModStr == 0, 1)
    traced = arcpy.sa.Times(Streams, os.path.join(optfolder, "basingrid"))
    arcpy.env.addOutputsToMap = False
    traced.save(os.path.join(optfolder, "ModStreams"))  # add it to TOC


    # *******************************************************************************************************
    # add subwatersheds to view -- addOutputsToMap

    basingrid = os.path.join(optfolder, "basingrid")
    flowacc = os.path.join(optfolder, "flowacc")
    flwdir = arcpy.sa.Times(os.path.join(optfolder, "flowdir"), basingrid)
    strlnk = arcpy.sa.StreamLink(os.path.join(optfolder, "ModStreams"), flwdir)
    zonemax = arcpy.sa.ZonalStatistics(strlnk, "Value", flowacc, "MAXIMUM", "NODATA")
    faccmax = arcpy.sa.Con(flowacc == zonemax, zonemax, arcpy.sa.IsNull(zonemax))
    outlets = arcpy.sa.Con(faccmax > 0, faccmax)
    outlets.save(os.path.join(optfolder, "outlets"))

    # if outlets are added by user then add mod stream and user specified outlets
    outlets_user   = os.path.join(optfolder, "outlets_user")
    if os.path.exists(outlets_user):
        aux_out1 = arcpy.RasterToPoint_conversion(outlets)
        aux_out2 = arcpy.RasterToPoint_conversion(outlets_user)
        aux_out3 = arcpy.Merge_management([aux_out1,aux_out2])
        outlets = arcpy.PointToRaster_conversion(aux_out3, "FID", os.path.join(optfolder, "outlets2"), "#", "#", basingrid)

    arcpy.env.extent = "MAXOF"
    subwshed = arcpy.sa.Watershed(flwdir, outlets, "VALUE")
    subwshed = arcpy.sa.Times(subwshed, basingrid)
    arcpy.env.addOutputsToMap = True
    arcpy.RasterToPolygon_conversion(subwshed, os.path.join(optfolder, "tmpsubwshd"),"NO_SIMPLIFY","VALUE")

    arcpy.env.extent = "MAXOF"
    tmpsub = os.path.join(optfolder, "tmpsubwshd.shp")
    subshed = os.path.join(optfolder, "subshed.shp")
    arcpy.Dissolve_management(tmpsub,subshed,"GRIDCODE","#","MULTI_PART","DISSOLVE_LINES")

    subrivers = os.path.join(optfolder, "subrivers.shp")
    # new stream link raster to handle added outlets
    if os.path.exists(outlets_user):
        StrmFDRgrid = arcpy.sa.Times(os.path.join(optfolder, "ModStreams"),flwdir)
        newLnkGrid = arcpy.sa.Watershed(StrmFDRgrid, outlets, "VALUE")
        newLnkGrid = arcpy.sa.Times(os.path.join(optfolder, "ModStreams"),newLnkGrid)
        arcpy.sa.StreamToFeature(newLnkGrid, flwdir, subrivers, "NO_SIMPLIFY")
        del aux_out1,aux_out2,aux_out3
    else:
        arcpy.sa.StreamToFeature(os.path.join(optfolder, "ModStreams"), flwdir, subrivers, "NO_SIMPLIFY")

    # *******************************************************************************************************
    # correct sub-basin indexing -- order of sub-basins has to be the same as
    # appearing in "subriver.shp"
    # *******************************************************************************************************
    target1 = os.path.join(optfolder, "subrivers.shp")
    joint1 = os.path.join(optfolder, "subshed.shp")
    spatial_join1 = os.path.join(optfolder, "subriver_subshed.shp")

    # add table using spatial join (target = subriver, join = subshed)
    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable(target1)
    fieldmappings.addTable(joint1)
    arcpy.SpatialJoin_analysis(target1, joint1, spatial_join1, "#", "#", fieldmappings, "HAVE_THEIR_CENTER_IN", "#", "#")

    # get subriver FID and centroid (or midpoint) XY values
    sub_cur = arcpy.SearchCursor(spatial_join1, "", "", "Shape;ARCID;GRIDCODE", "")
    x = []
    y = []
    for row in sub_cur:
        XMidPoint = row.shape.positionAlongLine(0.50, True).firstPoint.X
        x.append(XMidPoint)
        YMidPoint = row.shape.positionAlongLine(0.50, True).firstPoint.Y
        y.append(YMidPoint)

    del sub_cur

    xy = list(zip(x, y))
    # create a point shapefile, add fields "ARCID" & "GRIDCODE", and add midpoint XY values to field "SHAPE@"
    spatial_reference = arcpy.Describe(os.path.join(optfolder, "subriver_subshed.shp")).spatialReference
    arcpy.CreateFeatureclass_management(optfolder, "subriver_xy.shp", "POINT", "", "ENABLED", "DISABLED", spatial_reference)
    subriver_xy = os.path.join(optfolder, "subriver_xy.shp")

    # added on 09-30-2014 to debug field type error -- changed field type to "Double" with precision "10"
    if not len(arcpy.ListFields(subriver_xy, "ARCID")) > 0:
        arcpy.AddField_management(subriver_xy, "ARCID", "DOUBLE", 10, "")
        arcpy.AddField_management(subriver_xy, "GRIDCODE", "DOUBLE", 10, "")
    sr_xy = arcpy.da.InsertCursor(subriver_xy, ("SHAPE@XY"))
    for i in xy:
        sr_xy.insertRow([i])

    # update "ARCID" and "GRIDCODE" row values using "subriver_subshed" shapefile
    sub_cur = arcpy.SearchCursor(spatial_join1, "", "", "ARCID;GRIDCODE", "")
    ic = arcpy.UpdateCursor(subriver_xy, "", "", "ARCID;GRIDCODE", "")

    row1 = sub_cur.next()
    row2 = ic.next()

    while row1:
        row2.setValue("ARCID", row1.getValue("ARCID"))
        ic.updateRow(row2)
        row2.setValue("GRIDCODE", row1.getValue("GRIDCODE"))
        ic.updateRow(row2)
        row1 = sub_cur.next()
        row2 = ic.next()
    del sub_cur, ic

    # perform spatial join (target = subshed, join = subriver_xy) and lists extraction
    target2 = os.path.join(optfolder, "subshed.shp")
    joint2 = os.path.join(optfolder, "subriver_xy.shp")
    spatial_join2 = os.path.join(optfolder, "subshed_subriver.shp")

    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable(target2)
    fieldmappings.addTable(joint2)
    arcpy.SpatialJoin_analysis(target2, joint2, spatial_join2, "#", "#", fieldmappings, "INTERSECT", "#", "#")

    # create a polygon which will contain subriver FID, ARCID, and GRIDCODE (obtained from spatially joined shapefile)
    spatial_reference = arcpy.Describe(os.path.join(optfolder, "subshed_subriver.shp")).spatialReference
    arcpy.CreateFeatureclass_management(optfolder, "subshed.shp", "POLYGON", "", "DISABLED", "DISABLED",
                                        spatial_reference)

    sc = arcpy.da.SearchCursor(spatial_join2, ("ARCID", "GRIDCODE", "SHAPE@"))

    storage = []
    for row in sc:  # changed sF to sc [could be a mistake to have "sF"]
        storage.append(row)

    sortedlist = sorted(storage, key=lambda x: x[0])  # list sort is based on "ARCID"

    poly = os.path.join(optfolder, "subshed.shp")

    # added on 09-30-2014 to debug field type error -- changed field type to "Double" with precision "10"
    if not len(arcpy.ListFields(poly, "ARCID")) > 0:
        arcpy.AddField_management(poly, "ARCID", "DOUBLE", 10, "")
        arcpy.AddField_management(poly, "GRIDCODE", "DOUBLE", 10, "")

    arcpy.DeleteField_management(poly, "Id")
    fcout = arcpy.da.InsertCursor(poly, ("ARCID", "GRIDCODE", "SHAPE@"))

    for i in sortedlist:
        fcout.insertRow(i)
    del fcout, sc

    # convert newly created and indexed sub-watershed shapefile to raster
    arcpy.PolygonToRaster_conversion(subshed, "ARCID", os.path.join(optfolder, "subsheds"), "", "ARCID",
                                     os.path.join(optfolder, "dem"))  # priority field is "ARCID"


