"""A SUPPRIMER A LA FIN DU SPRINT

# Fichier servant à tester le bon fonctionnement des
# bindings de DTMIntersection. Cependant, il se peut qu'il
# ne soit pas necessaire de binder cette classe (utilisée que dans le cpp)"""

import rpc_c

rpc_c.DTMIntersection([1 for i in range(20)])
