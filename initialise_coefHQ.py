# -*- coding: utf-8 -*-
"""
Projet d'optimisation de M1 : fit d'un nuage de points par une hyperquadrique

Fonction pour l'initialisation des paramètres de l'hyperquadrique  

@author: Thomas Dietenbeck
"""

#%%
import numpy as np

def initialise_coefHQ(x_pts, y_pts, N) :
    
    """ Initialise les parametres de l'HQ 
    
        Un jeu de paramètres est déterminé à partir d'un polygone entourant le nuage de points.
        NB : ne fonctionne correctement que pour N >= 3
        
        Parameters
        ----------
        x_pts, y_pts: tables de float, coordonnees du nuage de points a fitter
        N: entier, nombre de termes de l'HQ à fitter
        
        Returns
        -------
        param_ini : table de float de dimension (N,4), paramètres pour initialiser le fit de l'HQ
        
    """
    # 1) Calcul du centre de gravite et des distances max et min
    xG, yG = np.mean(x_pts), np.mean(y_pts)
    dMax = 2 * np.sqrt( np.max( (x_pts - xG)**2 + (y_pts - yG)**2 ) )
    pad = 0.1 * dMax    # Marge pour s'assurer que les droites incluent tous les points

    i, piN = np.arange(N), np.pi / N
    R = ( dMax / 2 + pad ) / np.cos( piN )  # Rayon du cercle circonscrit au polygone
    D = 1.5 * ( dMax + pad )                # Distance entre les droites d'une paire
    # Definition des N sommets d'un polygone
    ptsPoly = np.zeros( (2, N + 1) )
    ptsPoly[0, :N] = xG + R * np.cos( (2 * i + 1) * piN )
    ptsPoly[1, :N] = yG + R * np.sin( (2 * i + 1) * piN )
    ptsPoly[0, -1], ptsPoly[1, -1] = ptsPoly[0, 0], ptsPoly[1, 0]
    # Definition de N points sur les lignes paralleles
    ptsPar = np.zeros( (2, N) )
    ptsPar[0, :] = xG - (D - R) * np.cos( 2 * (i+1) * piN )
    ptsPar[1, :] = yG - (D - R) * np.sin( 2 * (i+1) * piN )

    # Calcul des coefficients
    Det = ptsPoly[0, i] * ptsPoly[1, i+1]   - ptsPoly[0, i] * ptsPar[1, i] +    \
            ptsPar[0, i] * ptsPoly[1, i]    - ptsPar[0, i] * ptsPoly[1, i+1] +  \
            ptsPoly[0, i+1] * ptsPar[1, i]  - ptsPoly[0, i+1] * ptsPoly[1, i]

    param_ini = np.zeros((N, 4))
    param_ini[i, 0] = 2 * ( ptsPoly[1, i+1] - ptsPoly[1, i] ) / Det
    param_ini[i, 1] = 2 * ( ptsPoly[0, i] - ptsPoly[0, i+1] ) / Det
    param_ini[i, 2] = ( ptsPoly[0, i+1] * ptsPar[1, i]   - ptsPoly[0, i] * ptsPoly[1, i+1] +   \
                         ptsPar[0, i]    * ptsPoly[1, i]  - ptsPoly[0, i] * ptsPar[1, i] +      \
                         ptsPoly[0, i+1] * ptsPoly[1, i]  - ptsPar[0, i]  * ptsPoly[1, i+1] ) / Det
    param_ini[i, 3] = 4
    
    return param_ini
