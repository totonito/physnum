# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
#
# Ceci est un code d'Analyse des resultats obtenues avec le code 'ExerciceY_XXXX'.
# Le code d'analyse est base sur le standard python3 (on recommende python3.6 ou superieur).
# Pour installer les librairies on recommende d'utiliser pip3. Si pip3 n'est pas installe,
# il est possible de suivre la procedure d'un des links suivants:
#   https://linuxconfig.org/how-to-install-pip-on-ubuntu-18-04-bionic-beaver
#   https://linuxize.com/post/how-to-install-pip-on-ubuntu-18.04/
# 
# Ensuite, il faut installer les librairies: 
#   numpy 
#   matplotlib
#   scipy 
#   os
# methode d'installation conseille: utiliser la ligne de commande: 
#   pip3 install --user *nome-librairie*
# dans un terminal linux
#
# Pour utiliser le code d'analyse, il faut que son source (ce ficher) soit 
# dans le repertoire contenant le binaire 'ExerciceY_XXXX' et les 
# fichers d'output '.out'. Pour l'executer, il faut utiliser la ligne de commande
# suivantes dans le terminal linux:
#   python3 Analyse.py
#
# ----------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------
# Modules --------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import PlotResults # importer le module pour plotter les resultats

# ----------------------------------------------------------------------------------------------
# Parametres -----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

# Parametres pour lire et plotter les resultats d'une simulation -------------------------------

# Nome du fichier d'output a analyser
doFigures = True # si vrai les figures pour un fichier fileName sont produites 
delimiter = ' '
fileNames =['output.out']

# dictionaire contenant les parametres pour faire une figure
figure1 = {'NRows':1,'NCols':3} # plot la dynamique
figure1['NPlots'] = [1,1,1]
figure1['indexes'] = [[[0,2]],[[0,3]],[[2,3]]]
figure1['titles'] = ['Position','Vitesse','Phase']
figure1['xLables'] = ['t (s)','t (s)','x (m)']
figure1['yLables'] = ['x (m)','v (m/s)','v (m/s)']
figure1['normalisedPlot'] = [False,False,False]
figure1['axis'] = ['auto','auto','auto']
figure1['legend'] = [[],[],[]]
figure1['legendLocation'] = [3,4,3]
figure1['AnalyticalSolutions'] = []

# dictionaire contenant les parametres pour faire une figure
figure2 = {'NRows':1,'NCols':1} # plot la dynamique
figure2['NPlots'] = [1]
figure2['indexes'] = [[[0,4]]]
figure2['titles'] = ['Energie mechanique']
figure2['xLables'] = ['t (s)']
figure2['yLables'] = ['Emec (J)']
figure2['normalisedPlot'] = [False]
figure2['axis'] = ['auto']
figure2['legend'] = [[]]
figure2['legendLocation'] = [1]
figure2['AnalyticalSolutions'] = []

plt.rcParams["figure.figsize"] = (10,8)


# dictionaire contenant les parametres pour faire tous figures
figureParameters = {'plotData':[figure1, figure2]}
figureParameters['fontSize'] = 10
figureParameters['lineStyle'] = ['-','-','-','-','-','-','-','-','-','-','-','-']
figureParameters['lineWidth'] = 3
figureParameters['marker'] = ['.','.','.','.','.','.','.','.','.','.','.','.']
figureParameters['markerSize'] = 0
figureParameters['plotType'] =  ' ' #'loglog'


# ----------------------------------------------------------------------------------------------
# Lire et plotter les resultats de la simulation -----------------------------------------------
# ----------------------------------------------------------------------------------------------

# verifier s'il faut faire les figures
if(doFigures):
  # Chargement des donnees dans un numpy array
  # initialiser la class pour faire des plots
  plotResult = PlotResults.PlotResults(fileNameList=fileNames,figureParameters=figureParameters,\
  valuesDelimiter=delimiter)
  # faire des plots simples
  plotResult.SimplePlotFigures()

# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

