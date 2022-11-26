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
  #   fileNameList:	       (list(string)) list des nomes du fichier avec les valeurs
  #   values:		       (list(array))(Nfiles,Nx,Ny) tableau des valeus
  #   figureParameters:	       (dict) parametres pour faire des figures
  #   interpolationParameters: (dict) parametres pour les interpolations 
  #   valuesDelimiter:	       (string) signe pour delimiter les colonnes
  #					dans le fichier: default = ' '
  # outputs:
  def __init__(self,fileNameList=[],values=np.array([]),\
    figureParameters={},interpolationParameters={},valuesDelimiter=' '):

    # initialiser la liste des valeurs
    self.__values = []
    # verifier s'il faut lire ou sauvegarder les resultats
    if(len(fileNameList)!=0):
      # boucle sur les nomes des fichiers a lire
      for fileName in fileNameList:  
        # charger le tableau des resultats dans un numpy array
        self.__values.append(np.loadtxt(fileName,dtype=np.float64,delimiter=valuesDelimiter))
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
    #   fontSize:            (integer) figure font size
    #   lineStyle:           (list(string)) list of line stile for plotting
    #					    one for each value array
    #   lineWidth:           (integer) width of the line to plot
    #   marker:              (list(string)) list of types of marker to use:
    #					    one for each value array 
    #   markerSize:          (integer) size of the marker to use for plotting
    #   plotType:            (str) type of plot to be done, available:
    #			           1) linear, 2) loglog
    #   AnalyticalSolutions: (function)(NSubplots) fonction for solution comparison
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
  def setValues(self,values):

    # copier le tableau des valeurs
    self.__values = np.copy(values)

# ----------------------------------------------------------------------------------------------

  # Normaliser une colonne du tableau des valeurs
  # inputs:
  #   value: (array) array des valeurs a normaliser
  #   index: (int) index de la colonne a normaliser
  # outputs:
  #   normaliseValeur: (array) colonne des valeurs normalisee
  def normaliseFirstValue(self,value,index):

    # return the normalised value given a index
    return 100.e0*(value[:,index]-\
    value[0,index])/value[0,index]

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
        # select axis of subplot
        ax = self.selectAxisHandler(rowId,colsId,plotId,axList,figureData)
        # selectioner le type de plot
        plotFunction = self.selectPlotType(ax)
        # boucle sur le numero les plots
        for valueIndex in plotList:
          # boucle sur tous tableaux des figures
          for valueId,value in enumerate(self.__values):
            # verifier s'il faut normaliser les donnees 
            if(figureData['normalisedPlot'][plotId]):          
              # faire la figure
              plotFunction(value[:,valueIndex[0]],\
              self.normaliseFirstValue(value,valueIndex[1]),\
              linestyle=self.figureParameters['lineStyle'][valueId],\
              linewidth=self.figureParameters['lineWidth'],\
              marker=self.figureParameters['marker'][valueId],\
              markersize=self.figureParameters['markerSize'])
              # check if there is an analytical function for comparison
              if(len(figureData['AnalyticalSolutions'])==\
              figureData['NCols']*figureData['NRows']):
                # calculer la solution analytique
                analyticalSolution = \
                figureData['AnalyticalSolutions'][plotId](value[:,valueIndex[0]])
                # faire la figure
                plotFunction(value[:,valueIndex[0]],\
                analyticalSolution/valueIndex[1],\
                linestyle=self.figureParameters['lineStyle'][valueId],\
                linewidth=self.figureParameters['lineWidth'],\
                marker=self.figureParameters['marker'][valueId],\
                markersize=self.figureParameters['markerSize'])
            else:
              # faire la figure
              plotFunction(value[:,valueIndex[0]],\
              value[:,valueIndex[1]],\
              linestyle=self.figureParameters['lineStyle'][valueId],\
              linewidth=self.figureParameters['lineWidth'],\
              marker=self.figureParameters['marker'][valueId],\
              markersize=self.figureParameters['markerSize'])
              # check if there is an analytical function for comparison
              if(len(figureData['AnalyticalSolutions'])==\
              figureData['NCols']*figureData['NRows']):
                # calculer la solution analytique
                analyticalSolution = \
                figureData['AnalyticalSolutions'][plotId](value[:,valueIndex[0]])
                # faire la figure
                plotFunction(value[:,valueIndex[0]],
                analyticalSolution,\
                linestyle=self.figureParameters['lineStyle'][valueId],\
                linewidth=self.figureParameters['lineWidth'],\
                marker=self.figureParameters['marker'][valueId],\
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
  # cette methode genere des figures 2D a different iterations
  # depuis un flattened array. La structure des données 
  # doit etre: 0: colonne des iterations, 1: colonne index
  # du 1iere axe, 2: colonne 2ieme axe, 3-> valeurs 
  def plotSimple2DTimeFlattenedFigure(self):

    # load tools
    from ArrayTools import extractArrayFromValueReshape

    # importer les utiles pour generer des figures
    import matplotlib.pyplot as plt # importer pyplot depuis matplotlib
    from matplotlib import ticker # ticker for increasing number of ticks

    # initialiser les propriete des figures
    plt.rcParams.update({'font.size': self.figureParameters['fontSize']})
    # boucle sur les index a plotter
    for figureId,figureData in enumerate(self.figureParameters['plotData']):
      # intialiser une nouvelle figure avec subplots
      fig1,axList = plt.subplots(figureData['NRows'],figureData['NCols'], \
      sharey=False, sharex=False,num=figureId+1)
      # initialiser l'index des lignes
      rowId = 0
      # initialiser l'index des colonnes
      colsId = 0
      # seulement le 1iere valeur du self.__values est utilisee
      value = self.__values[0]
      # loop on the values to extract
      for extractValueId,extractValue in enumerate(figureData['extractValue']):
        # select axis of subplot
        ax = self.selectAxisHandler(rowId,colsId,extractValueId,axList,figureData)
        # selectioner le type de plot
        plotFunction = self.selectPlotType(ax)
        # extraire l'array 2D en fonction du 1iere index
        plotTable,xAxis,yAxis = extractArrayFromValueReshape(extractValue,\
        figureData['indexes'],figureData['tableSize'],self.__values[0])
        # fair le plot avec 3 inputs
        image = plotFunction(xAxis,yAxis,plotTable,\
        levels=self.figureParameters['contour_levels'])
        # ajouter le label des x
        ax.set_xlabel(figureData['xLables'],\
        fontsize=self.figureParameters['fontSize'])
        # ajouter le label des y
        ax.set_ylabel(figureData['yLables'],\
        fontsize=self.figureParameters['fontSize'])
        # ajouter un titre
        ax.set_title(\
        ''.join([figureData['titles'],' at value: ',str(extractValue)]),\
        fontsize=self.figureParameters['fontSize'])
        # ajouter colorbar	 
        fig1.colorbar(image,ax=ax)
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

  # cette methode genere des figures avec vecteurs 2D depuis des tableau de taille
  # (N+1)x(M+1). Les abscisses sont donnees par la ligne (1:M,0) par contre les 
  # ordonnees par la colonne (0,1:N)
  def plotSimple2D(self):

    # importer les utiles pour generer des figures
    import matplotlib.pyplot as plt # importer pyplot depuis matplotlib

    # initialiser les propriete des figures
    plt.rcParams.update({'font.size': self.figureParameters['fontSize']})
  
    # boucle sur les index a plotter
    for figureId,figureData in enumerate(self.figureParameters['plotData']):
      # seulement le 1iere valeur du self.__values est utilisee
      # intialiser une nouvelle figure avec subplots
      fig1,ax = plt.subplots(1,1, \
      sharey=False, sharex=False,num=figureId+1)
      # selectioner le type de plot
      plotFunction = self.selectPlotType(ax)
      value = self.__values[0]
      # extraire les abscisses
      xAxis = value[0,1:]
      # extraire les ordonees
      yAxis = value[1:,0]
      # create 2D image
      image = plotFunction(xAxis,yAxis,value[1:,1:],\
      levels=self.figureParameters['contour_levels'])
      # ajouter le label des x
      ax.set_xlabel(figureData['xLables'],\
      fontsize=self.figureParameters['fontSize'])
      # ajouter le label des y
      ax.set_ylabel(figureData['yLables'],\
      fontsize=self.figureParameters['fontSize'])
      # ajouter colorbar	 
      fig1.colorbar(image,ax=ax)
      # ajouter un titre
      ax.set_title(figureData['titles'],fontsize=self.figureParameters['fontSize'])

    # plotter le figures
    plt.show()

# ----------------------------------------------------------------------------------------------
  # cette methode genere des figures avec vecteurs 2D a different iterations
  # depuis un flattened array. La structure des données 
  # doit etre: 0: colonne des iterations, 1: colonne index
  # du 1iere axe, 2: colonne 2ieme axe, 3: valeurs
  # 1iere composante vecteur 4: 2ieme composante vecteur
  def plotSimple2DTimeFlattenedVectorFigure(self):
    # load tools
    from ArrayTools import extractArrayFromValueReshape

    # importer les utiles pour generer des figures
    import matplotlib.pyplot as plt # importer pyplot depuis matplotlib

    # initialiser les propriete des figures
    plt.rcParams.update({'font.size': self.figureParameters['fontSize']})
    # boucle sur les index a plotter
    for figureId,figureData in enumerate(self.figureParameters['plotData']):
      # intialiser une nouvelle figure avec subplots
      fig1,axList = plt.subplots(figureData['NRows'],figureData['NCols'], \
      sharey=False, sharex=False,num=figureId+1)
      # initialiser l'index des lignes
      rowId = 0
      # initialiser l'index des colonnes
      colsId = 0
      # seulement le 1iere valeur du self.__values est utilisee
      value = self.__values[0]
      for extractValueId,extractValue in enumerate(figureData['extractValue']):
        # select quiver as axis of subplot
        ax = self.selectAxisHandler(rowId,colsId,extractValueId,axList,figureData)
        # selectioner le type de plot
        plotFunction = self.selectPlotType(ax)
        # x vector list
        indexList = figureData['indexes'][0:4]
        # extraire l'array 2D en fonction du 1iere index
        XVectorArray,xAxis,yAxis = extractArrayFromValueReshape(extractValue,\
        indexList,figureData['tableSize'],self.__values[0])
        # y vector list
        indexList = figureData['indexes'][0:3]
        indexList.append(figureData['indexes'][4])
        # extraire l'array 2D en fonction du 2ieme index
        YVectorArray,gxAxis,yAxis = extractArrayFromValueReshape(extractValue,\
        indexList,figureData['tableSize'],self.__values[0])
        # colormap vector list
        indexList = figureData['indexes'][0:3]
        indexList.append(figureData['indexes'][5])
        # extraire l'array 2D colormap en fonction du 3ieme index
        CVectorArray,gxAxis,yAxis = extractArrayFromValueReshape(extractValue,\
        indexList,figureData['tableSize'],self.__values[0])
	# faire le plot
        image=ax.quiver(xAxis,yAxis,XVectorArray,YVectorArray,CVectorArray)
        # ajouter le label des x
        ax.set_xlabel(figureData['xLables'],\
        fontsize=self.figureParameters['fontSize'])
        # ajouter le label des y
        ax.set_ylabel(figureData['yLables'],\
        fontsize=self.figureParameters['fontSize'])
        # ajouter colorbar	 
        fig1.colorbar(image,ax=ax)
        # ajouter un titre
        ax.set_title(\
        ''.join([figureData['titles'],' at value: ',str(extractValue)]),\
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
    elif(self.figureParameters['plotType']=='semilogy'):
      # retourner une fonction de plot semilogy
      return ax.semilogy
    elif(self.figureParameters['plotType']=='semilogx'):
      # retourner une fonction de plot semilogx
      return ax.semilogx
    elif(self.figureParameters['plotType']=='contour'):
      # retourner une fonction pour plotter les contours
      return ax.contour
    elif(self.figureParameters['plotType']=='contourf'):
      # retourner une fonction pour plotter les contours
      # avec filling
      return ax.contourf
    else:
      # retourner une fonction de plot lineaire
      return ax.plot

# ----------------------------------------------------------------------------------------------

  # This method select the axis handler
  # inputs:
  #   rowId:	  (integer32) subplot row index
  #   colsId:	  (integer32) subplot column index
  #   plotId:	  (integer32) subplot index
  #   axList:	  (list(axis)) list of subplot axis
  #   figureData: (figureParameters) figure parameters structure
  # output
  #   ax:	  (axis) axis handler
  def selectAxisHandler(self,rowId,colsId,plotId,axList,figureData):

    # pointer vers l'ax
    if((figureData['NRows']==1) and (figureData['NCols']==1)):
      ax = axList
    elif((figureData['NRows']==1) or (figureData['NCols']==1)):
      ax = axList[plotId]
    else:
      ax = axList[rowId,colsId]

    # return ax handler
    return ax

# ----------------------------------------------------------------------------------------------

