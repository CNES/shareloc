# -*- coding: utf-8 -*-
"""
Created on Tue Jun 09 18:44:35 2020

@author: gresloud
"""
import numpy as np
from shareloc.readwrite import  lit_hd_bsq
from shareloc.math_utils import interpol_bilin


#-------------------------------------------------------------------------------
class gld_xH:
    """ multi H direct localization grid handling class 
    """
    def __init__(self,fichier_gld_bsq, format_gri = 'bsq'):
        self.fichier_in = fichier_gld_bsq
        self.format    = format_gri
        self.lig0      = None
        self.col0      = None
        self.nblig     = None
        self.nbcol     = None
        self.paslig    = None
        self.pascol    = None
        self.repter    = None
        self.nbalt     = None
        self.index_alt = {}
        self.gld_lon   = None
        self.gld_lat   = None
        self.alts_down = []
        self.ligmax    = None
        self.colmax    = None
        self.charge()

    def charge(self):
        """
        Lit les grilleset les entetes et definit 2 "cubes" de donnees
        - gld_lon : [alt,lig,col]
        - gld_lat : [alt,lig,col]
        Les grilles bsq sont rangees de H0 a Hx, avec alt montante
        Nous on les range en ordre descendant

        """
        if self.format=='bsq':

            dico_a_lire = {'nblig':('LINES',int),'nbcol':('COLUMNS',int),'bpp':('BITS PER PIXEL',int),\
            'nbalt':('NB ALT',int),'pascol':('PAS COL',float),'paslig':('PAS LIG',float),\
            'col0':('COL0',float),'lig0':('LIG0',float),\
            'repter':('REFERENTIEL TERRESTRE',str)}

            nom_hd      = self.fichier_in[:-4]+'1.hd'
            dico_hd = lit_hd_bsq(nom_hd,dico_a_lire)

            for var in dico_hd:
                setattr(self,var,dico_hd[var])

            """renvoie une structure 3D [i_alt][l,c]"""
            gld_lon = np.zeros((self.nbalt,self.nblig,self.nbcol))
            gld_lat = np.zeros((self.nbalt,self.nblig,self.nbcol))

            codage = float

            for i in range(self.nbalt):
                k = self.nbalt-i
                nom_gri_lon = self.fichier_in[:-4]+str(k)+'.c1'
                nom_gri_lat = self.fichier_in[:-4]+str(k)+'.c2'
                nom_hd      = self.fichier_in[:-4]+str(k)+'.hd'

                gld_lon[i,:,:] = np.fromfile(nom_gri_lon,dtype=codage).reshape((self.nblig,self.nbcol))
                gld_lat[i,:,:] = np.fromfile(nom_gri_lat,dtype=codage).reshape((self.nblig,self.nbcol))

                dico_hd = lit_hd_bsq(nom_hd,{'index':('ALT INDEX',int),'alt':('ALTITUDE',float)})
                self.index_alt[dico_hd['index']] = dico_hd['alt']

            self.gld_lon = gld_lon
            self.gld_lat = gld_lat
            self.alts_down = [self.index_alt[_] for _ in range(int(self.nbalt-1),-1,-1)]
            self.ligmax = self.lig0 + self.paslig*(self.nblig-1)
            self.colmax = self.col0 + self.pascol*(self.nbcol-1)
        else:

            print("format de mnt non reconnu")


    def checkCubeMNT(self,Visee,mnt):
        #Visee: (n,3):
        PointB = None
        dH3D = None
        (ui,vi,zi,hi) = ([],[],[],[])
        ViseeMNT = mnt.TersToMnts(Visee)
        # -----------------------------------------------------------------------
        # Nombre d'intersections valides trouvees
        nbi = 0
        # -----------------------------------------------------------------------
        # On boucle sur les plans du cube MNT
        for f in range(6):
            # -----------------------------------------------------------------------
            # Initialisation du sommet de la visee
            sB = ViseeMNT[0,:]
            # -----------------------------------------------------------------------
            # Initialisation de la position par / au plan
            #print mnt.plans[f,:-1]
            #posB = (mnt.plans[f,:-1]*sB).sum() - mnt.d[f]
            #print posB
            posB = mnt.eq_plan(f,sB)
            # -----------------------------------------------------------------------
            # On boucle sur les segments de la visee et on controle
            # si on traverse ou pas la face courante f du cube MNT
            for p in range(self.nbalt):
                # -----------------------------------------------------------------------
                # Transfert du point B dans le point A
                posA = posB
                sA = sB.copy()
                # -----------------------------------------------------------------------
                # Reinit du point B
                sB = ViseeMNT[p,:] #dg on itere sur les differents points de la visee
                #print sB
                # -----------------------------------------------------------------------
                # Initialisation de la position par / au plan
                posB = mnt.eq_plan(f,sB)
                #print 'posAposB',posA,posB
                # -----------------------------------------------------------------------
                # Test d'intersection : posA et posB de signes opposes
                if (posA * posB <= 0):
                    # -----------------------------------------------------------------------
                    if (not(posA) and  not(posB)):
                        # Trop de solutions !! (A et B sont sur le plan) #comment on le gere ??
                        continue
                    # -----------------------------------------------------------------------
                    elif not(posA):
                        # A est solution (il est sur le plan)
                        ui.append( sA[0] )
                        vi.append( sA[1] )
                        zi.append( sA[2] )
                        hi.append( p - 1)
                    # -----------------------------------------------------------------------
                    elif not(posB):
                        # B est solution (il est sur le plan)
                        ui.append( sB[0] )
                        vi.append( sB[1] )
                        zi.append( sB[2] )
                        hi.append( p )
                    # -----------------------------------------------------------------------
                    else:
                        # -----------------------------------------------------------------------
                        # A et B sont de part et d'autre du plan
                        # Coefficients d'interpolation de l'intersection
                        # entre A et B
                        cA =  posB / (posB - posA)
                        cB = -posA / (posB - posA)
                        # Affectation ou interpolation
                        # NB : pour eviter les pb lors du test
                        #      <estSurCube> (voir + loin)
                        # . coordonn?e <u> (ligne)
                        # -----------------------------------------------------------------------
                        if (f < 2):
                            ui.append( mnt.d[f] )
                        # -----------------------------------------------------------------------
                        else:
                            ui.append( cA * sA[0] + cB * sB[0])
                        # -----------------------------------------------------------------------
                        # . coordonnee <v> (colonne)
                        if (f > 1) and (f < 4):
                            vi.append( mnt.d[f] )
                        else:
                            vi.append( cA * sA[1] + cB * sB[1] )
                        # -----------------------------------------------------------------------
                        # . coordonn?e <z> (altitude)
                        if (f > 3):
                            zi.append( mnt.d[f] )
                        # -----------------------------------------------------------------------
                        else:
                            zi.append( cA * sA[2] + cB * sB[2] )
                        # . coordonn?e <h> (abscisse visee)
                        hi.append( p - cA ) #index non entier de l'intersection
                    # -----------------------------------------------------------------------
                    # Incrementation du nombre d'intersections trouvees
                    nbi+=1
                    # -----------------------------------------------------------------------
                    # Passage a la face du cube suivante
                    break

        # -----------------------------------------------------------------------
        # Tri des points le long de la visee (il y en a au moins deux)
        #on les range par ordre decroissant en altitude
        for p in range(nbi):
            for q in range(p+1,nbi):
                if (hi[q] < hi[p]):
                    dtmp  = ui[p]
                    ui[p] = ui[q]
                    ui[q] = dtmp
                    dtmp  = vi[p]
                    vi[p] = vi[q]
                    vi[q] = dtmp
                    dtmp  = zi[p]
                    zi[p] = zi[q]
                    zi[q] = dtmp
                    dtmp  = hi[p]
                    hi[p] = hi[q]
                    hi[q] = dtmp

        # -----------------------------------------------------------------------
        # Filtrage des points non situes sur le cube
        p = 0
        while (p < nbi):
            #test a l'interieur du cube
            estSurCube = (ui[p] >= mnt.d[0]) and (ui[p] <= mnt.d[1]) and \
                         (vi[p] >= mnt.d[2]) and (vi[p] <= mnt.d[3]) and \
                         (zi[p] >= mnt.d[4]) and (zi[p] <= mnt.d[5])
            if not(estSurCube):
                # On translate tous les points suivants (on ecrase ce point non valide)
                for q in range(p+1,nbi):
                    ui[q - 1] = ui[q]
                    vi[q - 1] = vi[q]
                    zi[q - 1] = zi[q]
                    hi[q - 1] = hi[q]
                nbi-=1
            else:
                p+=1
        # -----------------------------------------------------------------------
        # Pas de solution si 0 ou 1 seul point trouve (on a tangente le cube)
        if (nbi < 2):
            bTrouve = False
            return (True, bTrouve, PointB, dH3D)
        # -----------------------------------------------------------------------
        # Il ne reste que 2 points donc on traverse le cube
        # LAIG-FA-MAJA-2168-CNES: plus de filtrage sur les point identiques. Il peut y avoir un nombre depoints > 2
        # Initialisation du point courant
        # Coordonnees MNT
        PointMnt = np.zeros(3)
        PointMnt[0] = ui[0]
        PointMnt[1] = vi[0]
        PointMnt[2] = zi[0]
        #PointMnt est la premiere intersection avec le cube (lig, col)
        # -----------------------------------------------------------------------
        # h dans gld 3D
        dH3D = hi[0]
        #dH3D correspond a l'index (non entier) d'interpolation en h
        # -----------------------------------------------------------------------
        # Coordonnees terrain
        PointB = mnt.MntToTer(PointMnt)
        #PointB est le point Terreain (lon,lat)
        # -----------------------------------------------------------------------
        # Fin, retour
        bTrouve = True
        return (True, bTrouve, PointB,dH3D)


    #----------------------------------------------------------------
    def intersection(self,Visee, PointB, dH3D, mnt):
        """
        fonction d'intersection mnt
        (Visee, H3D, PointB, dH3D, PointR)
        """
        ViseeMNT   = mnt.TersToMnts(Visee)
        PointB_MNT = mnt.TerToMnt(PointB)
        PointR     = np.zeros(3)
        (npl,_)     = Visee.shape
        H3D = range(npl,-1,-1)

        p1 = PointB_MNT.copy() #[p1[0],p1[1],p1[2]]

        dH3D_p1 = dH3D

        nu = mnt.nl
        nv = mnt.nc

        #1 - Initilialisation et tests prealables
        #   1.1 - Test si le sommet est au-dessu du MNT
        #       - Calcul de l'altitude du MNT ? la position du sommet
        h1 = mnt.MakeAlti(p1[0], p1[1])
        #       - Calcul de l'ecart d'altitude au MNT
        d1 = p1[2] - h1

        #       - Test si le nouveau point haut est au dessus du MNT
        if (d1 < 0):
        #       - Point situ? en dessous du MNT
        #          . ceci signifie que la vis?e rentre dans le MNT par le c?t?
        #          . donc en dessous, pas de solution
            bTrouve = False
            return (True,bTrouve,PointR)


        #   1.2 - Initialisation du rang du premier sommet de la visee
        i0 = int(np.floor(dH3D_p1))

        #   1.3 - Initialisation du point de depart (dans p2)
        p2      = PointB_MNT.copy()
        dH3D_p2 = dH3D

        #2. - Boucle sur les plans de grille
        while (i0 < ViseeMNT.size - 1):
            #2.1 - Initialisation du sommet courant de la visee
            u0 = ViseeMNT[i0][0]
            v0 = ViseeMNT[i0][1]
            z0 = ViseeMNT[i0][2]
            z1 = ViseeMNT[i0 + 1][2]

            #2.2 - Initialisation de la visee MNT
            VM = ViseeMNT[i0 + 1] - ViseeMNT[i0]

            #2.3 - Test si visee verticale
            if (VM[0]==0 and VM[1]==0):
                #2.3.1 - La vis?e est verticale :
                #    - Calcul de l'altitude du MNT ? la position du sommet
                h1 = self.MakeAlti(u0, v0)

                #    Test si le plan suivant est en dessous du MNT
                if (ViseeMNT[i0 + 1][2] <= h1):
                    #Init point de sortie
                    p1[0] = u0
                    p1[1] = v0
                    p1[2] = h1
                    bTrouve = True
                    mnt.MntToTer(p1, PointR)
                    return (True,bTrouve,PointR)
                else:
                    #Positionnement sur le sommet suivant
                    i0+=1
            else:
                #2.3.2 - La visee n'est pas verticale :
                #         elle va donc survoler le MNT
                #         . on peut donc poursuivre
                #         . reste ? d?montrer que la vis?e se crashe sur le MNT...
                #
                # Initialisation du point de d?part
                # Initialisation de son abscisse sur la vis?e
                a2 = dH3D_p2 - i0

                # Correction d'un bug FA DG 10, a la reinitialisation a2 peut valoir 1
                # Ce qui a pour consequence de sauter au segment suivant de la visee
                # Ainsi, on peut perdre des points sur une tranche d'altitude
                # Pour etre sur de scanner le segment de visee, on reinitialise
                # le point au depart du segment, soit a2 = 0
                if (a2 >= 1.):
                    a2 = 0.

                # Initialisation de la premi?re maille MNT intersect?e
                #  - Initialisation des indices de la maille
                uc = int(np.floor(p2[0]))
                vc = int(np.floor(p2[1]))

                # NB :    pr?caution avant de d?marrer :
                #        . on se met du bon cote de la maille
                #        . en principe, on ne doit pas sortir du MNT
                # On rentre par le bas, la maille MNT est la precedente
                if ((p2[0] == uc) and (VM[0] < 0)):
                    uc-=1

                # On rentre par la gauche, la maille MNT est la precedente
                if ((p2[1] == vc) and (VM[1] < 0)):
                    vc-=1

                # LDD - On est deja en dehors des limites, on s'arrete
                if (not((a2 < 1) and (uc > -1) and (uc < (nu - 1)) and (vc > -1) and (vc < (nv - 1)) ) and (a2 < 1)):
                    bTrouve = False
                    (True,bTrouve,PointR)

                # Boucle de recherche iterative de la maille intersectee
                while ((a2 < 1) and (uc > -1) and (uc < (nu - 1)) and (vc > -1) and (vc < (nv - 1))):
                    # - Altitudes min et max de la maille
                    hi = mnt.Zmin_cell[uc,vc]
                    hs = mnt.Zmax_cell[uc,vc]

                    # - Transfert : le point bas devient le point haut
                     #a1 = a2;
                    # p1 devient p2
                    p1 = p2.copy()
                    dH3D_p1 = dH3D_p2.copy()

                    # 4.2 - D?termination d'un nouveau point bas
                    #      - Test d'orientation de la vis?e
                    if (not(VM[0])):
                        # 4.2.1 - La vis?e est orientee pile poil est-ouest
                        #   p2[0] = p1[0] ; // inutile, est deja initialise
                        #       - Test d'orientation de la visee
                        if (VM[1] < 0):
                            # 4.2.1.1 - La vis?e part plein ouest
                            p2[1] = vc
                            vc-=1
                        else:
                            # 4.2.1.2 - La vis?e part plein est
                            vc+=1
                            p2[1] = vc

                        a2    = (p2[1] - v0) / VM[1]
                        p2[2] = z0 + a2 * VM[2]

                    elif (not(VM[1])):
                        # 4.2.2 - La vis?e est orient?e nord-sud
                        #  p2[1] = p1[1] ;
                        #       - Test d'orientation de la visee
                        if (VM[0] < 0):
                            # 4.2.2.1 - La vis?e part plein nord
                            p2[0] = uc
                            uc-=1
                        else:
                            # 4.2.2.2 - La vis?e part plein sud
                            uc+=1
                            p2[0] = uc

                        a2    = (p2[0] - u0) / VM[0]
                        p2[2] = z0 + a2 * VM[2]

                    else:
                        # 4.2.3 - La vis?e est quelconque
                        #            - D?termination du cot? de sortie
                        if ((VM[0] < 0) and (VM[0] <= VM[1]) and (VM[0] <= -VM[1])):
                            # 4.2.3.1 - Vis?e principalement orient?e nord
                            #             - Intersection avec le c?t? nord
                            a2    = (uc - u0) / VM[0]
                            p2[1] = v0 + a2 * VM[1]

                            if ((p2[1] > vc) and (p2[1] < (vc + 1))):
                                # La vis?e sort par le nord
                                p2[0] = uc
                                p2[2] = z0 + a2 * VM[2]
                                uc-=1

                            elif (p2[1] < vc):
                                # La vis?e sort par l'ouest
                                a2 = (vc - v0) / VM[1]
                                p2[0] = u0 + a2 * VM[0]
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                vc-=1

                            elif (p2[1] > (vc + 1)):
                                # La vis?e sort par l'est
                                vc+=1
                                a2 = (vc - v0) / VM[1]
                                p2[0] = u0 + a2 * VM[0]
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]

                            elif (p2[1] == vc):
                                # La vis?e sort par le coin nord-ouest
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                uc-=1
                                vc-=1

                            elif (p2[1] == (vc + 1)):
                                # La vis?e sort par le coin nord-est
                                p2[0] = uc
                                vc+=1
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                uc-=1

                        elif ((VM[1] > 0) and (VM[1] >= VM[0]) and (VM[1] >= -VM[0])):
                            # 4.2.3.2 - Vis?e principalement orient?e est
                            #         - Intersection avec le c?t? est
                            a2    = (vc + 1 - v0) / VM[1]
                            p2[0] = u0 + a2 * VM[0]

                            if ((p2[0] > uc) and (p2[0] < (uc + 1))):
                                #  La vis?e sort par l'est
                                    vc+=1
                                    p2[1] = vc
                                    p2[2] = z0 + a2 * VM[2]

                            elif (p2[0] < uc):
                                # La vis?e sort par le nord
                                p2[0] = uc
                                a2 = (uc - u0) / VM[0]
                                p2[1] = v0 + a2 * VM[1]
                                p2[2] = z0 + a2 * VM[2]
                                uc-=1

                            elif (p2[0] > (uc + 1)):

                                # La vis?e sort par le sud
                                uc+=1
                                p2[0] = uc
                                a2 = (uc - u0) / VM[0]
                                p2[1] = v0 + a2 * VM[1]
                                p2[2] = z0 + a2 * VM[2]

                            elif (p2[0] == uc):
                                # La vis?e sort par le coin nord-est
                                vc+=1
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                uc-=1

                            elif (p2[0] == (uc + 1)):
                                # La vis?e sort par le coin sud-est
                                uc+=1
                                vc+=1
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]


                        elif ((VM[0] > 0) and (VM[0] >= VM[1]) and (VM[0] >= -VM[1])):
                            # 4.2.3.3 - Vis?e principalement orient?e sud
                            #         - Intersection avec le c?t? sud
                            a2    = (uc + 1 - u0) / VM[0]
                            p2[1] = v0 + a2 * VM[1]

                            if ((p2[1] > vc) and (p2[1] < (vc + 1))):
                                # La vis?e sort par le sud
                                uc+=1
                                p2[0] = uc
                                p2[2] = z0 + a2 * VM[2]

                            elif (p2[1] < vc):
                                # La vis?e sort par l'ouest
                                a2 = (vc - v0) / VM[1]
                                p2[0] = u0 + a2 * VM[0]
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                vc-=1

                            elif (p2[1] > vc + 1):
                                # La vis?e sort par l'est
                                vc+=1
                                a2 = (vc - v0) / VM[1]
                                p2[0] = u0 + a2 * VM[0]
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]

                            elif (p2[1] == vc):
                                # La vis?e sort par le coin sud-ouest
                                uc+=1
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                vc-=1

                            elif (p2[1] == vc + 1):
                                # La vis?e sort par le coin sud-est
                                uc+=1
                                vc+=1
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]


                        elif ((VM[1] < 0) and (VM[1] <= VM[0]) and (VM[1] <= -VM[0])):
                            #  4.2.3.4 - Vis?e principalement orient?e ouest
                            #          - Intersection avec le c?t? ouest
                            a2    = (vc - v0) / VM[1]
                            p2[0] = u0 + a2 * VM[0]

                            if ((p2[0] > uc) and (p2[0] < uc + 1)):
                                # La vis?e sort par l'ouest
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                vc-=1

                            elif (p2[0] < uc):
                                # La vis?e sort par le nord
                                p2[0] = uc
                                a2 = (uc - u0) / VM[0]
                                p2[1] = v0 + a2 * VM[1]
                                p2[2] = z0 + a2 * VM[2]
                                uc-=1

                            elif (p2[0] > (uc + 1)):
                                # La vis?e sort par le sud
                                uc+=1
                                p2[0] = uc
                                a2 = (uc - u0) / VM[0]
                                p2[1] = v0 + a2 * VM[1]
                                p2[2] = z0 + a2 * VM[2]

                            elif (p2[0] == uc):
                                # La vis?e sort par le coin nord-ouest
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                uc-=1
                                vc-=1

                            elif (p2[0] == (uc + 1)):
                                # La vis?e sort par le coin sud-ouest
                                uc+=1
                                p2[0] = uc
                                p2[1] = vc
                                p2[2] = z0 + a2 * VM[2]
                                vc-=1


                    # LDD - Verification des bornes min et max de la "couche"
                    bIntersect = False

                    if (p2[2] > z0): # On est remonte trop haut, et ca c'est pas bon !!!
                        bIntersect = not(((p1[2]>hs) and (z0>hs)) or ((p1[2]<hi) and (z0<hi)))

                    elif (p2[2] < z1): # On est descendu trop bas, et ca c'est pas bon non plus !!! (m?me si c'est d?j? plus logique)
                        bIntersect = not(((p1[2]>hs) and (z1>hs)) or ((p1[2]<hi) and (z1<hi)))

                    else:
                        bIntersect = not(((p1[2]>hs) and (p2[2]>hs)) or ((p1[2]<hi) and (p2[2]<hi)))

                    # 5. Test d'intersection de la vis?e avec le cube
                    if (bIntersect):
                        # Il y a intersection entre la vis?e et le cube
                        # 5.1 - Altitudes du MNT
                        h1 = mnt.MakeAlti(p1[0], p1[1])
                        h2 = mnt.MakeAlti(p2[0], p2[1])

                        # 5.2 - Diff?rences d'altitude avec le MNT
                        d1 = p1[2] - h1
                        d2 = p2[2] - h2

                        # 5.3 - Test d'intersection avec le MNT
                        if (d1 * d2 <= 0):
                            # Il y a intersection entre la vis?e et le MNT
                            # 5.3.1 - Calcul de la solution approch?e
                            d2 = 2 * mnt.TOL_Z # Init de d2 > TOL_Z
                            ua = p2[0]
                            va = p2[1]
                            za = h2

                            while (abs(d2) > mnt.TOL_Z):
                                # 5.3.1.1 - Coefficient d'interpolation lin?aire de h
                                ch = (p1[2] - h1) / ((h2 - h1) - (p2[2] - p1[2]))

                                # 5.3.1.2 - Position du point interpole
                                ua = p1[0] + ch * (p2[0] - p1[0])
                                va = p1[1] + ch * (p2[1] - p1[1])
                                za = p1[2] + ch * (p2[2] - p1[2])

                                # 5.3.1.3 - Altitude du point interpole
                                zv = mnt.MakeAlti(ua, va)

                                # 5.3.1.4 - Ecart d'altitude au point interpole
                                d2 = zv - za

                                # 5.3.1.5 - Mise a jour
                                if (d2 < 0):
                                    # Mise a jour du point haut
                                    p1[0] = ua
                                    p1[1] = va
                                    p1[2] = za
                                    h1 = zv

                                else:
                                    # Mise ? jour du point bas
                                    p2[0] = ua
                                    p2[1] = va
                                    p2[2] = za
                                    h2 = zv


                            #// Fin, retour
                            p1[0] = ua
                            p1[1] = va
                            p1[2] = za

                            bTrouve = True
                            PointR = mnt.MntToTer(p1)
                            return (True,bTrouve,PointR)

                # Fin boucle sur les mailles

                # Test si on est toujours dans le cube MNT
                if (a2 >= 1):
                    # Changement de plan
                    i0+=1

                    # Chargement dans p2 du nouveau sommet
                    p2 = ViseeMNT[i0].copy()
                    dH3D_p2 = H3D[i0].copy()

                else:

                    # LDD - On a boucl? sur les mailles, on n'a rien trouv? et on n'a pas atteint le plan suivant
                    # Ca veut dire qu'on sort de l'emprise, pas la peine de continuer
                    bTrouve = False
                    return (True,bTrouve,PointR)
            # Fin cas general (visee non verticale)

        # Fin boucle sur les sommets

        # Fin, retour
        bTrouve = False
        return (True,PointR)

    def fct_locdir_h(self,lig,col,alt):
        """fonction de localisation a altitude constante"""
        #faire une controle sur lig / col !!!!
        # 0.5 < lig < ligmax
        (kh,kb) = self.renvoie_indices_grilles_alt(alt)
        altbas  = self.alts_down[kb]
        althaut = self.alts_down[kh]
        dh = (alt - altbas)/(althaut - altbas)
        mats =  [self.gld_lon[kh:kb+1,:,:],self.gld_lat[kh:kb+1,:,:]]
        P     = np.zeros(3)
        P[2] = alt
        dl = (lig - self.lig0)/self.paslig
        dc = (col - self.col0)/self.pascol
        [vlon,vlat] = interpol_bilin(mats,self.nblig,self.nbcol,dl,dc)
        P[0] = (dh*vlon[0] +(1-dh)*vlon[1])
        P[1] = (dh*vlat[0] +(1-dh)*vlat[1])
        return P

    def fct_locdir_mnt(self,lig,col, mnt):
        """fonction de localisation sur mnt"""
        visee = np.zeros((3,self.nbalt))
        vislonlat = self.fct_interp_visee_unitaire_gld(lig,col)
        visee[0,:] = vislonlat[0]
        visee[1,:] = vislonlat[1]
        visee[2,:] = self.alts_down
        v = visee.T
        (code1, code2, PointB, dH3D) = self.checkCubeMNT(v,mnt)
        (code3,code4,Point_mnt) = self.intersection(v, PointB, dH3D,mnt)
        return Point_mnt

    def fct_locdir_mntopt(self,lig,col, mnt):
        """fonction de localisation sur mnt"""
        visee = np.zeros((3,self.nbalt))
        vislonlat = self.fct_interp_visee_unitaire_gld(lig,col)
        visee[0,:] = vislonlat[0]
        visee[1,:] = vislonlat[1]
        visee[2,:] = self.alts_down
        v = visee.T
        (code, code2, PointB, dH3D) = self.checkCubeMNT(v,mnt)
        #(code,code4,Point_mnt) = self.intersection(v, PointB, dH3D,mnt)
        return code

    def fct_interp_visee_unitaire_gld(self,lig,col):
        dl = (lig - self.lig0)/self.paslig
        dc = (col - self.col0)/self.pascol
        mats =  [self.gld_lon,self.gld_lat]
        res = interpol_bilin(mats,self.nblig,self.nbcol,dl,dc)
        return res

    def fct_interp_gld(self,nblig,nbcol,nbalt=None):
        """renvoie un gld_xH ??"""
        if not nbalt:
            nbalt = self.nbalt
            list_alts = self.alts_down
        else:
            list_alts = np.linspace(self.alts_down[0],self.alts_down[-1],nbalt)

        gld_lon = np.zeros((nbalt,nblig,nbcol))
        gld_lat = np.zeros((nbalt,nblig,nbcol))
        """genere un cube de visee interpole de nlig/ncol visee"""
        #lig_max = self.lig0 + self.paslig * (self.nblig-1)
        #col_max = self.col0 + self.pascol * (self.nbcol-1)

        paslig = (self.ligmax - self.lig0)/(nblig-1)
        pascol = (self.colmax - self.col0)/(nbcol-1)

        for k,alt in enumerate(list_alts):

            res = self.fct_gld_h(self.lig0,self.col0,paslig,pascol,nblig,nbcol,alt)
            gld_lon[k] = res[0]
            gld_lat[k] = res[1]
        return gld_lon,gld_lat


    def fct_gld_mnt(self,lig0,col0,paslig,pascol,nblig,nbcol,mnt):
        """
        fonction de calcul de grille de loc directe sur MNT
        """
        gldmnt = np.zeros((3,nblig,nbcol))
        visee = np.zeros((3,self.nbalt))
        for i in range(nblig):
            for j in range(nbcol):
                col = col0 + pascol*j
                lig = lig0 + paslig*i
                vislonlat = self.fct_interp_visee_unitaire_gld(lig,col)
                visee[0,:] = vislonlat[0]
                visee[1,:] = vislonlat[1]
                visee[2,:] = self.alts_down
                v = visee.T
                (code1, code2, PointB, dH3D) = self.checkCubeMNT(v,mnt)
                (code3,code4,PointR) = self.intersection(v, PointB, dH3D,mnt)
                gldmnt[:,i,j] = PointR
        return gldmnt

    def renvoie_indices_grilles_alt(self,alt):
        """renvoie les indices sup et inf des plans d'altitude encadrant une altitude
        Les grilles gld ne sont pas sensees etre a des hauteurs regulieres"""
        if alt > self.alts_down[0] :
            (indicehaut,indicebas) = (0,0)
        elif alt < self.alts_down[-1]:
            (indicehaut,indicebas) = (self.nbalt-1,self.nbalt-1)
        else:
            i = 0
            while i < self.nbalt and self.alts_down[i] >= alt:
                i+=1
            if i== self.nbalt: #pour gerer alt min
                i = self.nbalt - 1
            indicebas = i     #indice grille bas
            indicehaut = i-1   #indice grille haut
        return (indicehaut,indicebas)


    def fct_gld_h(self,lig0,col0,paslig,pascol,nblig,nbcol,alt):
        """
        fonction de calcul de grille de loc directe a altitude constante
        => a faire: controler les extrapolations
        """
        (kh,kb) = self.renvoie_indices_grilles_alt(alt)
        gldalt  = np.zeros((3,nblig,nbcol))
        altbas  = self.alts_down[kb]
        althaut = self.alts_down[kh]

        dh = (alt - altbas)/(althaut - altbas)
        mats =  [self.gld_lon[kh:kb+1,:,:],self.gld_lat[kh:kb+1,:,:]]
        P     = np.zeros(3)
        P[2] = alt
        for i in range(nblig):
            lig = lig0 + paslig*i
            dl = (lig - self.lig0)/self.paslig
            for j in range(nbcol):
                col = col0 + pascol*j
                dc = (col - self.col0)/self.pascol
                [vlon,vlat] = interpol_bilin(mats,self.nblig,self.nbcol,dl,dc)
                P[0] = (dh*vlon[0] +(1-dh)*vlon[1])
                P[1] = (dh*vlat[0] +(1-dh)*vlat[1])
                gldalt[:,i,j] = P
        return gldalt

    def init_pred_loc_inv(self,nblig_pred=3,nbcol_pred=3):
        """Initialise les polynomes de prediction de localisation inverse de toutes les barrettes de la pdv
        Il s'agit de 4 polynomes calcules sur des grilles 5x5 a hmin et hmax sous
        la forme:
        col_min = a0 + a1*lon + a2*lat + a3*lon**2 + a4*lat**2 + a5*lon*lat
        lig_min = b0 + b1*lon + b2*lat + b3*lon**2 + b4*lat**2 + b5*lon*lat
        col_max = a0 + a1*lon + a2*lat + a3*lon**2 + a4*lat**2 + a5*lon*lat
        lig_max = b0 + b1*lon + b2*lat + b3*lon**2 + b4*lat**2 + b5*lon*lat
        Les coefficients sont determines par moindres carres sur des grandeurs normalisees -1 et 1
        """
        nb_alt     = 2
        nb_coeff   = 6
        nb_mes     = nbcol_pred*nblig_pred

        col_norm    = np.linspace(-1.0,1.0,nbcol_pred)
        lig_norm    = np.linspace(-1.0,1.0,nblig_pred)
        gcol_norm, glig_norm = np.meshgrid(col_norm,lig_norm)
        glon,glat = self.fct_interp_gld(nblig_pred,nbcol_pred,nb_alt)

        #normalisation des variables
        (glon_min,glon_max) = (glon.min(),glon.max())
        (glat_min,glat_max) = (glat.min(),glat.max())
        lon_ofset = (glon_max + glon_min)/2.0
        lon_scale = (glon_max - glon_min)/2.0
        lat_ofset = (glat_max + glat_min)/2.0
        lat_scale = (glat_max - glat_min)/2.0

        glon_norm = (glon - lon_ofset) / lon_scale
        glat_norm = (glat - lat_ofset) / lat_scale


        col_ofset = (self.colmax + self.col0)/2.0
        col_scale = (self.colmax - self.col0)/2.0
        lig_ofset = (self.ligmax + self.lig0)/2.0
        lig_scale = (self.ligmax - self.lig0)/2.0

        glon2     = glon_norm*glon_norm
        glonlat   = glon_norm*glat_norm
        glat2     = glat_norm*glat_norm

        Amin    = np.zeros((nb_mes,nb_coeff))
        Amax    = np.zeros((nb_mes,nb_coeff))
        Bcol    = np.zeros((nb_mes,1))
        Blig    = np.zeros((nb_mes,1))

        #resolution des moindres carres
        imes = 0
        for ilig in range(nblig_pred):
            for icol in range(nbcol_pred):

                Bcol[imes]       = gcol_norm[ilig,icol]
                Blig[imes]       = glig_norm[ilig,icol]
                Amin[imes,0]     = 1.0
                Amax[imes,0]     = 1.0
                Amin[imes,1]     = glon_norm[1,ilig,icol]
                Amax[imes,1]     = glon_norm[0,ilig,icol]
                Amin[imes,2]     = glat_norm[1,ilig,icol]
                Amax[imes,2]     = glat_norm[0,ilig,icol]
                Amin[imes,3]     = glon2[1,ilig,icol]
                Amax[imes,3]     = glon2[0,ilig,icol]
                Amin[imes,4]     = glat2[1,ilig,icol]
                Amax[imes,4]     = glat2[0,ilig,icol]
                Amin[imes,5]     = glonlat[1,ilig,icol]
                Amax[imes,5]     = glonlat[0,ilig,icol]
                imes +=1

        #Calcul des coeffcients
        matAmin       = np.matrix(Amin)
        matAmax       = np.matrix(Amax)
        tAAmin        = matAmin.T*matAmin
        tAAmax        = matAmax.T*matAmax
        tAAmin_inv    = tAAmin.I
        tAAmax_inv    = tAAmax.I

        coef_col_min = tAAmin_inv*matAmin.T*Bcol
        coef_lig_min = tAAmin_inv*matAmin.T*Blig
        coef_col_max = tAAmax_inv*matAmax.T*Bcol
        coef_lig_max = tAAmax_inv*matAmax.T*Blig

        setattr(self,'pred_col_min',coef_col_min.A1)
        setattr(self,'pred_lig_min',coef_lig_min.A1)
        setattr(self,'pred_col_max',coef_col_max.A1)
        setattr(self,'pred_lig_max',coef_lig_max.A1)
        setattr(self,'pred_ofset_scale_lon',  [lon_ofset , lon_scale])
        setattr(self,'pred_ofset_scale_lat',  [lat_ofset , lat_scale] )
        setattr(self,'pred_ofset_scale_lig',  [lig_ofset , lig_scale] )
        setattr(self,'pred_ofset_scale_col',  [col_ofset , col_scale] )
        code_return = 0
        return code_return

    #-------------------------------------------------------------------------------------------------------------------------
    def fct_locinv_pred(self,lon,lat,alt=0):
        """Cette fonction evalue le polynome predicteur en un point et fait
        l'interporpolation altimetrique"""
        seuil_extrapol = 20.0
        extrapol = False
        altmin = self.alts_down[-1]
        altmax = self.alts_down[0]

        #normalisation
        lon_n = (lon - self.pred_ofset_scale_lon[0]) / self.pred_ofset_scale_lon[1]
        lat_n = (lat - self.pred_ofset_scale_lat[0]) / self.pred_ofset_scale_lat[1]
        if abs(lon_n) > (1+seuil_extrapol/100.0):
            #print "Attention, en extrapolation de 20% en longitude:",lon_n
            extrapol = True
        if abs(lat_n) > (1+seuil_extrapol/100.0):
            #print "Attention, en extrapolation de 20% en latitude:",lat_n
            extrapol = True

        #application polynome
        vect_sol = np.array([1,lon_n,lat_n,lon_n**2,lat_n**2,lon_n*lat_n])
        col_min = ((self.pred_col_min*vect_sol).sum() * self.pred_ofset_scale_col[1]) + self.pred_ofset_scale_col[0]
        lig_min = ((self.pred_lig_min*vect_sol).sum() * self.pred_ofset_scale_lig[1]) + self.pred_ofset_scale_lig[0]
        col_max = ((self.pred_col_max*vect_sol).sum() * self.pred_ofset_scale_col[1]) + self.pred_ofset_scale_col[0]
        lig_max = ((self.pred_lig_max*vect_sol).sum() * self.pred_ofset_scale_lig[1]) + self.pred_ofset_scale_lig[0]

        hx = (alt-altmin)/(altmax-altmin)
        col = (1-hx)*col_min + hx*col_max
        lig = (1-hx)*lig_min + hx*lig_max

        if lig > self.ligmax: lig = self.ligmax
        if lig < self.lig0:   lig = self.lig0
        if col > self.colmax: col = self.colmax
        if col < self.col0:   col = self.col0
        return (lig,col,extrapol)

    #-------------------------------------------------------------------------------------------------------------------------
    def fct_loc_inv_mat_dp(self,lig,col,alt=0):
        """calcul de la matrice de derivee partielles permettant de passer d'une
        variation sol a une variation en coordonnees images.
        Renvoie la matrice M telle que
        [dcol,dlig]T = M x [dlon,dlat]T
        dlon/dlas sont en microrad
        Ce calcul se fait sur les noeuds de la grille initiale
        Cette matrice est necessaire pour l'inversion de la loc directe lors du processus iteratif de loc inverse"""

        dl = (lig - self.lig0)/self.paslig
        dc = (col - self.col0)/self.pascol
        il = int(np.floor(dl))
        ic = int(np.floor(dc))
        (kh,kb) = self.renvoie_indices_grilles_alt(alt)

        lon_h00 = self.gld_lon[kh,il  ,ic]
        lon_h01 = self.gld_lon[kh,il  ,ic+1]
        lon_h10 = self.gld_lon[kh,il+1,ic]

        lon_b00 = self.gld_lon[kb,il  ,ic]
        lon_b01 = self.gld_lon[kb,il  ,ic+1]
        lon_b10 = self.gld_lon[kb,il+1,ic]

        lat_h00 = self.gld_lat[kh,il  ,ic]
        lat_h01 = self.gld_lat[kh,il  ,ic+1]
        lat_h10 = self.gld_lat[kh,il+1,ic]

        lat_b00 = self.gld_lat[kb,il  ,ic]
        lat_b01 = self.gld_lat[kb,il  ,ic+1]
        lat_b10 = self.gld_lat[kb,il+1,ic]

        dlon_ch = np.deg2rad(lon_h01 - lon_h00)/self.pascol
        dlon_cb = np.deg2rad(lon_b01 - lon_b00)/self.pascol
        dlon_lh = np.deg2rad(lon_h10 - lon_h00)/self.paslig
        dlon_lb = np.deg2rad(lon_b10 - lon_b00)/self.paslig

        dlat_ch = np.deg2rad(lat_h01 - lat_h00)/self.pascol
        dlat_cb = np.deg2rad(lat_b01 - lat_b00)/self.pascol
        dlat_lh = np.deg2rad(lat_h10 - lat_h00)/self.paslig
        dlat_lb = np.deg2rad(lat_b10 - lat_b00)/self.paslig

        hx = (alt - self.alts_down[kb])/(self.alts_down[kh]-self.alts_down[kb])

        dlon_c = ((1-hx)*dlon_cb + (hx)*dlon_ch)*1e6
        dlat_c = ((1-hx)*dlat_cb + (hx)*dlat_ch)*1e6
        dlon_l = ((1-hx)*dlon_lb + (hx)*dlon_lh)*1e6
        dlat_l = ((1-hx)*dlat_lb + (hx)*dlat_lh)*1e6
        det = dlon_c*dlat_l - dlon_l*dlat_c
        if abs(det) > 0.000000000001:
            Matdp = np.matrix([[dlat_l,-dlon_l],[-dlat_c,dlon_c]])/det
        else:
            print("determinant nul")
        return Matdp
    #-------------------------------------------------------------------------------------------------------------------------
    def fct_locinv(self,P,nb_iterations = 15):
        """Fonction de localisation inverse
        - calcul de la prediction col_0,lig_0
        - calcul de la loc directe lon_0, lat_0
        Puis processus iteratif:
        - calcul de l'erreur au sol dlon et dlat
        - calcul de la correction dcol dlig correspondante
        - calcul de la loc directe lon_i, lat_i

        P de taille [lon, lat, alt]
        si P de taille 2, alt = 0
        """
        (lon,lat) = (P[0],P[1])
        try:
            alt = P[2]
        except:
            alt = 0

        deg2mrad = np.deg2rad(1.0)*1e6
        k=0
        coslon = np.cos(np.deg2rad(lat))
        Rtx = 1e-12*6378000**2
        (lig_i,col_i,extrapol) = self.fct_locinv_pred(lon,lat,alt)
        erreur_m2    = 10.0
        point_valide = 0
        if not extrapol:
            #Processus iteratif
            #while erreur > seuil:1mm
            while (erreur_m2 > 1e-6) and (k < nb_iterations):
                #print k,lig_i,col_i
                P = self.fct_locdir_h(lig_i,col_i,alt)
                dlon_microrad = (P[0] - lon)*deg2mrad
                dlat_microrad = (P[1] - lat)*deg2mrad
                erreur_m2 = Rtx*(dlat_microrad**2+(dlon_microrad*coslon)**2)
                dsol = np.matrix([dlon_microrad,dlat_microrad]).T
                mat_dp = self.fct_loc_inv_mat_dp(lig_i,col_i,alt)
                dimg = mat_dp*dsol
                col_i += -dimg[0,0]
                lig_i += -dimg[1,0]
                k +=1
                point_valide = 1
        return (lig_i,col_i,point_valide)

    #-------------------------------------------------------------------------------

def fct_coloc(gld_xH_src, gld_hX_dst, mnt, \
        l0_src, c0_src, paslig_src, pascol_src, nblig_src, nbcol_src):

    gricoloc = np.zeros((3,nblig_src,nbcol_src))
    for l in range(nblig_src):
        lig = l0_src + paslig_src*l
        for c in range(nbcol_src):
            col = c0_src + pascol_src*c

            Psol = gld_xH_src.fct_locdir_mnt(lig,col, mnt)
            Pdst = gld_hX_dst.fct_locinv(Psol)
            gricoloc[:,l,c] = Pdst
    return gricoloc

