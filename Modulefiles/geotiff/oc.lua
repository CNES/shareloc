-- -*- lua -*-  
-- Aide du module accessible avec la commande module help 
help(   
[[ 
gdal_plugins des outils communs  
]]) 
 
-- Information du modulefile 
local os_disponible = "rh7" 
local nom           = "gdal_plugins" 
local version       = "20180125" 
local installation  = "29/07/2019"
local informations  = system_information() 
local rhos          = informations.os 


-- Variable du modulefile  
local qtispackdir = "/softs/projets/oc/qtis_pack_64_RH7.20180125/"
local home=pathJoin(qtispackdir,"BIBEXT_1/GEOTIFF_V1.3.0")

prepend_path("LD_LIBRARY_PATH",pathJoin(home,"/lib"))
prepend_path("PATH",pathJoin(home,"/bin"))



