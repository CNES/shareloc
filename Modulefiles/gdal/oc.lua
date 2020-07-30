-- -*- lua -*-  
-- Aide du module accessible avec la commande module help 
help(   
[[ 
	gdal des outils communs
]]) 
 
-- Information du modulefile 
local os_disponible = "rh7" 
local nom           = "gda" 
local version       = "20180125" 
local installation  = "29/07/2019"
local informations  = system_information() 
local rhos          = informations.os 


-- Variable du modulefile  
local qtispackdir = "/softs/projets/oc/qtis_pack_64_RH7.20180125/"
local home=pathJoin(qtispackdir,"BIBEXT_1/GDAL_V1.11.2_full/")

setenv("GDAL_DIR",home)
prepend_path("LD_LIBRARY_PATH",pathJoin(home,"/lib"))

setenv("GDAL_DATA",pathJoin(home,"/share/gdal/"))
prepend_path("PATH",pathJoin(home,"/bin"))
--export PYTHONPATH=$gdal_root/lib/python2.7/site-packages/:$PYTHONPATH


-- Action du modulefile  
setenv("GDALHOME",home)
prepend_path("PATH",pathJoin(home,"/bin"))
prepend_path("LD_LIBRARY_PATH",pathJoin(home,"/lib"))




-- Dépendance du modulefile  
conflict("gdal")






 


