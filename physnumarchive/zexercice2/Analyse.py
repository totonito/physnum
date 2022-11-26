# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
#
# Ceci est un code d'Analyse des resultats obtenues avec le code 'Exercice2'.
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
# dans le repertoire contenant le binaire 'Exercice2.exe' et les 
# fichers d'output '.out'. Pour l'executer, il faut utiliser la ligne de commande
# suivantes dans le terminal linux:
#   python3 Analyse.py
#
# ----------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------
# Modules -----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
import PlotResults # importer le module pour plotter les resultats
import SchemeAnalysisSuite

# ----------------------------------------------------------------------------------------------
# Parametres -----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

# Parametres pour lire et plotter les resultats d'une simulation -------------------------------

# Nome du fichier d'output a analyser
doFigures = True # si vrai les figures pour un fichier fileName sont produites 
fileName = 'output.out'
delimiter = ' '
# dictionaire contenant les parametres pour faire une figure
figure1 = {'NRows':2,'NCols':2} # plot la dynamique
figure1['NPlots'] = [1,2,2,1]
figure1['indexes'] = [[[1,2]],[[0,1],[0,2]],[[4,5]],[[0,4],[0,5]]]
figure1['titles'] = ['Position X-Y du proton','Temps-position du ballon',\
'Vitesse X-Y du ballon','Temps-vitesse du ballon']
figure1['xLables'] = ['x (m)','t (s)','vx (m/s)','t (s)']
figure1['yLables'] = ['z (m)','position (m)','vz (m/s)','vitesse (m/s)']
figure1['normalisedPlot'] = [False,False,False,False]
figure1['axis'] = ['equal','auto','equal','auto']
figure1['legend'] = [[],['x','y'],[],['vx','vy']]
figure1['legendLocation'] = [0,2,0,2]

figure2 = {'NRows':1,'NCols':1} # la conservation de l'energie
figure2['NPlots'] = [1]
figure2['indexes'] = [[[0,8]]]
figure2['titles'] = ['Energie: erreur de conservation']
figure2['xLables'] = ['t (s)']
figure2['yLables'] = ['erreur %']
figure2['normalisedPlot'] = [True]
figure2['axis'] = ['auto']

# dictionaire contenant les parametres pour faire tous figures
figureParameters = {'plotData':[figure1,figure2]}
figureParameters['fontSize'] = 12
figureParameters['lineStyle'] = '-'
figureParameters['lineWidth'] = 3
figureParameters['marker'] = 'x'
figureParameters['markerSize'] = 0
figureParameters['plotType'] = ' '#'loglog'

# Parametres pour executer le test de convergence ----------------------------------------------

# Variables definies par l'utilisateur
doConvergenceTest = True # si vrai le test de convergence est execute
# Parametres du fichiers d'input
inputs = {'tfin':8.1991e-08,'k':0.0,\
'B0':4.0,'B1':0.0,'B2':0.0,\
'Ex':0.0,'Ey':0.0,'Ez':0.0,\
'x0':0.e0,'y0':0.e0,'z0':0.e0,\
'vx0':0.e0,'vy0':5.e5,'vz0':0.e5,\
'L':1e-2,'q':1.6022e-19,'m':1.6726e-27,\
'schema':'E','alpha':0.5,'nsteps':2000,\
'output':'output','sampling':1,'tol':1.e-5,'maxit':1000}
# Parametre pour le test de convergence
convergenceParametres = dict()
convergenceParametres['fileDelimiter'] = ' '
convergenceParametres['binaryFileName'] = 'Exercice2_solution.exe'
convergenceParametres['inputFileName'] = 'configuration.in.2'
convergenceParametres['numberOfPointForRegression'] = 3
convergenceParametres['meshGenerator'] = 'power10'
convergenceParametres['meshValues'] = [2,5]
convergenceParametres['keyForConvergence'] = 'nsteps'
convergenceParametres['keyForOutputFileName'] = 'output'
convergenceParametres['interpolationValues'] = [75]
convergenceParametres['interpolationIndexes'] = [0]
convergenceParametres['indexList'] = [1,2,4,5]
convergenceParametres['analyticalSolution'] = [2.996820045514575e-12,8.843781816117237e-08,\
33.886182532065142,4.999999988517266e+05] #x(tfinale),y(tfinale),vx(tfinale),vy(tfinale)

convergenceParametres['inputParameters'] = inputs

# parametres pour generer la figure avec l'etude de convergence
figure1 = {'NRows':2,'NCols':2} # plot la dynamique
figure1['NPlots'] = [1,1,1,1]
figure1['indexes'] = [[[0,1]],[[0,2]],[[0,3]],[[0,4]]]
figure1['titles'] = ['Convergence en X','Convergence en Y',\
'Convergence en VX','Convergence en VY']
figure1['xLables'] = ['nombre de pas de temps','nombre de pas de temps',\
'nombre de pas de temps','nombre de pas de temps']
figure1['yLables'] = ['erreur','erreur','erreur','erreur']
figure1['normalisedPlot'] = [False,False,False,False]
figure1['axis'] = ['auto','auto','auto','auto']
figure1['legend'] = []
figure1['legendLocation'] = [1,1,1,1]

# dictionaire contenant les parametres pour faire tous figures
convergenceFigureParameters = {'plotData':[figure1]}
convergenceFigureParameters['fontSize'] = 12
convergenceFigureParameters['lineStyle'] = '-'
convergenceFigureParameters['lineWidth'] = 3
convergenceFigureParameters['marker'] = 'x'
convergenceFigureParameters['markerSize'] = 1
convergenceFigureParameters['plotType'] = 'loglog'


# ----------------------------------------------------------------------------------------------
# Lire et plotter les resultats de la simulation -----------------------------------------------
# ----------------------------------------------------------------------------------------------

# verifier s'il faut faire les figures
if(doFigures):
  # Chargement des donnees dans un numpy array
  #   0) temps
  #   1) position x
  #   2) position y
  #   3) position z
  #   4) vitesse x
  #   5) vitesse y
  #   6) vitesse y
  #   7) moment magnetique
  #   8) energy
  # initialiser la class pour faire des plots
  plotResult = PlotResults.PlotResults(fileName=fileName,figureParameters=figureParameters,\
  valuesDelimiter=delimiter)
  # faire des plots simples
  plotResult.SimplePlotFigures()


# ----------------------------------------------------------------------------------------------
# Executer le test de convergence --------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

# verifier qu'il faut faire le test de convergence
if(doConvergenceTest):

  # initialise un objet SchemeAnalysisSuite
  schemeAnalysis = SchemeAnalysisSuite.SchemeAnalysisSuite(\
  convergenceParameters=convergenceParametres,\
  figureParameters=convergenceFigureParameters)

  # executer le test de convergence
  schemeAnalysis.doConvergenceTest()

# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

