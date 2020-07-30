ml euclidium


#loc directe
LIG=100.5
COL=50.5
ALT=100
echo "direct localization LIG=${LIG} COL=${COL} ALT=${ALT}"
locdirecte_externe  -grille -path ./grilles_gld_xH/ -file P1BP--2017030824934340CP -ci ${LIG} ${COL} -alt ${ALT} -v -rep GRS80:G-D/:H-M


LIG=50.5
COL=10.0
ALT=10
echo "direct localization LIG=${LIG} COL=${COL} ALT=${ALT}"
locdirecte_externe  -grille -path ./grilles_gld_xH/ -file P1BP--2017030824934340CP -ci ${LIG} ${COL} -alt ${ALT} -v -rep GRS80:G-D/:H-M



