-- -*- lua -*-  
-- Variable du modulefile  

local version       = "20180125" 
local qtispackdir = "/softs/projets/oc/qtis_pack_64_RH7.20180125/"


local currentdir = lfs.currentdir()
execute {cmd="cd "..qtispackdir,modeA={"load"}}

-- commandes
execute {cmd="source ./init.sh",modeA={"load"}}
--if (mode() == "unload") then
--  LmodError("source ./init.sh : commande non reversible ")
--end


execute {cmd="euclidium_v12_0_init",modeA={"load"}}
execute {cmd="proj4_init",modeA={"load"}}
execute {cmd="stdlib_init",modeA={"load"}}
execute {cmd="gdlib_init",modeA={"load"}}
execute {cmd="xerces_init",modeA={"load"}}
execute {cmd="hdf_init",modeA={"load"}}
execute {cmd="hdf5_init",modeA={"load"}}
execute {cmd="cag_init",modeA={"load"}}

execute {cmd="cd "..currentdir,modeA={"load"}}

LmodMessage("ATTENTION les dépencdances à gdal gdal_plugins et geotiff on été désactivées") 

--depend("gdal/oc")
--depend("gdal_plugins/oc")
--depend("geotiff/oc")
