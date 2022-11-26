# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
#
# Ceci est un module du code d'Analyse des resultats des simulations qui a comme objectifs
# le generation de figures et l'interpolation des resultats.
# Le code d'analyse est base sur le standard python3 (on recommende python3.6 ou superieur).
# Pour installer les librairies on recommende d'utiliser pip3. Si pip3 n'est pas installe,
# il est possible de suivre la procedure d'un des links suivants:
#   https://linuxconfig.org/how-to-install-pip-on-ubuntu-18-04-bionic-beaver
#   https://linuxize.com/post/how-to-install-pip-on-ubuntu-18.04/
# 
# Ensuite, il faut installer les librairies: 
#   numpy 
#   matplotlib
# methode d'installation conseille: utiliser la ligne de commande: 
#   pip3 install --user *nome-librairie*
# dans un terminal linux
#
# Pour utiliser le code d'analyse, il faut que son source (ce ficher) soit 
# dans le repertoire contenant le main Analyse.py. 
#
# ----------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------
# Librairies -----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
import numpy as np # importer numpy pour tout le module

# definire la class PlotResultats
class PlotResults():

  # Definir le constructeur de la class PlotResultats
  # inputs:
  #   fileName:		       (string) nome du fichier avec les valeurs
  #   values:		       (array)(Nx,Ny) tableau des valeus
  #   figureParameters:	       (dict) parametres pour faire des figures
  #   interpolationParameters: (dict) parametres pour les interpolations 
  #   valuesDelimiter:	       (string) signe pour delimiter les colonnes
  #					dans le fichier: default = ' '
  # outputs:
  def __init__(self,fileName='',values=np.array([]),\
    figureParameters={},interpolationParameters={},valuesDelimiter=' '):

    # verifier s'il faut lire ou sauvegarder les resultats
    if(len(fileName)): 
      # charger le tableau des resultats dans un numpy array
      self.__values = np.loadtxt(fileName,dtype=np.float64,delimiter=valuesDelimiter)
    else:
      # sauvegarder les donnees
      self.__values = values

    # sauvegarder les parametres pour generer des figures, parametres:
    #   plotData:       list(dict) contient tous informations pour faire des plots
    #                          chaque cles est une figure differents. Chaque cle
    #                                contient un dictionaire avec les donnees:
    #                                NRows:   (int) numero de lignes de subplots 
    #				     NCols:   (int) numero de colonnes de subplots
    #				     NPlots:  (list(int)) numero de lignes dans le meme
    #						    subplot pour chaque subplot
    #				     indexes: (list(list(list(int))))(NSubplots,NPlots,2)
    #						    liste d'indexes des colonnes 
    #						    a imprimer. La premier liste definit
    #						    la deuxieme les lignes sur le meme subplot
    #						    la troisieme les index sous la forme
    #						    [xId,yId] where xId is the index 
    #						    of des abscisses et yId celui des ordonnees
    #				     title:   (list(string))(NSubplots) liste de titre pour
    #						    chauque subplot
    #				     xLables: (list(string))(NSubplots) liste de string containant
    #						    nome des x-axes pour chaque subplot
    #				     yLables: (list(string))(NSubplots) liste de string containant
    #						    nome des y-axes pour chaque subplot
    #				     normalisedPlot: (list(string))(NSubplots) liste de bool
    #						    si vrais le plot est normalise
    #				     axis: (list(string))(NSubplots) axis aspect ratio
    #				     legend: (list(list(string)))(NSubplots,NPlots)
    # 				     		    add a legend to subplot
    #				     legendLocation: (list(string))(NSubplots,NPlots)
    # 				     		     location of the legend
    #
    #   fontSize:       (integer) figure font size
    #   lineStyle:      (string) line stile for plotting
    #   lineWidth:      (integer) width of the line to plot
    #   marker:         (string) type of marker to use
    #   markerSize:     (integer) size of the marker to use for plotting
    #   plotType:       (str) type of plot to be done, available:
    #			      1) linear, 2) loglog
    self.figureParameters = figureParameters

    # sauvegarder les parametres pour faire des interpolations, parametres:
    #   typeInterpolation: (string) le type d'interpolation a utiliser
    #				    available: lineaire
    #   condition:	   (dict) type of condition for finding interpolation index
    #				  available: key: nearValue -> liste de valeurs
    #				  d'interpolation
    #   indexes:    (list(list(int)))(NInterpolations,2) list of indexes to be interpolated
    #				     they should have the form 
    #				     [xId,yId] where xId is the index
    #				     of des abscisses et yId celui des ordonnees
    #   
    self.interpolationParameters = interpolationParameters

    # liste d'index des points pour effecturer les interpolations forme list(list)
    self.__indexInterpolationList = [[]]

    # valeurs interpolees: liste abscisse ordonees des valeurs interpolees
    self.__interpolatedValues = [[]]

# ----------------------------------------------------------------------------------------------

    # destructeur de la class
    def __del__(self):

      # imprimer le nome de la class
      print(self.__class__.__name__ ,'is deleted')

# ----------------------------------------------------------------------------------------------

  # cette methode rendre le tableau des valeurs
  # inputs:
  # outputs:
  #   values: (array)(Nx,Ny) result table
  def getValues(self):

    # rendre le tableau des valeurs
    return self.__values

# ----------------------------------------------------------------------------------------------

  # cette methode copie le tableau des valeurs
  # inputs:
  #   values: (array)(Nx,Ny) result table
  # outputs:
  def setValues(self,results):

    # copier le tableau des valeurs
    self.__values = np.copy(values)

# ----------------------------------------------------------------------------------------------

  # Normaliser une colonne du tableau des valeurs
  # inputs:
  #   index: (int) index de la colonne a normaliser
  # outputs:
  #   normaliseValeur: (array) colonne des valeurs normalisee
  def normaliseFirstValue(self,index):

    # return the normalised value given a index
    return 100.e0*(self.__values[:,index]-\
    self.__values[0,index])/self.__values[0,index]

# ----------------------------------------------------------------------------------------------
  # cette methode genere des figures depuis les valeurs d'entrees
  # et le parametres.
  def SimplePlotFigures(self):

    # importer les utiles pour generer des figures
    import matplotlib.pyplot as plt # importer pyplot depuis matplotlib

    # initialiser les propriete des figures
    plt.rcParams.update({'font.size': self.figureParameters['fontSize']})
    # boucle sur les figures
    for figureId,figureData in enumerate(self.figureParameters['plotData']):
      # intialiser une nouvelle figure avec subplots
      fig1,axList = plt.subplots(figureData['NRows'],figureData['NCols'], \
      sharey=False, sharex=False,num=figureId+1)
      # initialiser l'index des lignes
      rowId = 0
      # initialiser l'index des colonnes
      colsId = 0
      # boucle sur les indexes
      for plotId,plotList in enumerate(figureData['indexes']):
        # pointer vers l'ax
        if((figureData['NRows']==1) and (figureData['NCols']==1)):
          ax = axList
        elif((figureData['NRows']==1) or (figureData['NCols']==1)):
          ax = axList[plotId]
        else:
          ax = axList[rowId,colsId]
        # selectioner le type de plot
        plotFunction = self.selectPlotType(ax)
        # boucle sur le numero les plots
        for valueIndex in plotList:
          # verifier s'il faut normaliser les donnees 
          if(figureData['normalisedPlot'][plotId]):          
            # faire la figure
            plotFunction(self.__values[:,valueIndex[0]],\
            self.normaliseFirstValue(valueIndex[1]),\
            linestyle=self.figureParameters['lineStyle'],\
            linewidth=self.figureParameters['lineWidth'],\
            marker=self.figureParameters['marker'],\
            markersize=self.figureParameters['markerSize'])
          else:
            # faire la figure
            plotFunction(self.__values[:,valueIndex[0]],\
            self.__values[:,valueIndex[1]],\
            linestyle=self.figureParameters['lineStyle'],\
            linewidth=self.figureParameters['lineWidth'],\
            marker=self.figureParameters['marker'],\
            markersize=self.figureParameters['markerSize'])
        # selctioner les axes
        ax.set_aspect(figureData['axis'][plotId], adjustable='box')
        # ajouter le label des x
        ax.set_xlabel(figureData['xLables'][plotId],\
        fontsize=self.figureParameters['fontSize'])
        # ajouter le label des y
        ax.set_ylabel(figureData['yLables'][plotId],\
        fontsize=self.figureParameters['fontSize'])
        # ajouter un titre
        ax.set_title(figureData['titles'][plotId],\
        fontsize=self.figureParameters['fontSize'])
        # ajouter une grille
        ax.grid()
        # verifier que la liste des legendes existe
        if('legend' in figureData):
          # verifier que la list de la legende n'est pas vide
          if(len(figureData['legend'][plotId])!=0):
            # ajouter la legende
            ax.legend(figureData['legend'][plotId],\
            loc=figureData['legendLocation'][plotId],\
            fontsize=self.figureParameters['fontSize'])
         
        # mis a jour index colonnes
        colsId = colsId+1
        # verifier la valeur de l'index des colonnes
        if(colsId>=figureData['NCols']):
          # mis a jour index lignes
          rowId = rowId+1
          # initialiser index colonnes
          colsId = 0

    # plotter les figures
    plt.show()

# ----------------------------------------------------------------------------------------------

  # cette methode selection le type de plot
  # inputs:
  #   ax: (axis) axis handler
  # outputs:
  #   plot: (func) function for plotting
  def selectPlotType(self,ax):

    # selectioner le type de plot
    if(self.figureParameters['plotType']=='loglog'):
      # retourner une fonction de plot loglog
      return ax.loglog
    else:
      # retourner une fonction de plot lineaire
      return ax.plot

# ----------------------------------------------------------------------------------------------

