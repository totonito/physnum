# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
#
# Ceci est un module du code d'Analyse des schemas numeriques qui a comme objectifs
# d'etudier leur convergence et proprietes de conservation.
# Le code d'analyse est base sur le standard python3 (on recommende python3.6 ou superieur).
# Pour installer les librairies on recommende d'utiliser pip3. Si pip3 n'est pas installe,
# il est possible de suivre la procedure d'un des links suivants:
#   https://linuxconfig.org/how-to-install-pip-on-ubuntu-18-04-bionic-beaver
#   https://linuxize.com/post/how-to-install-pip-on-ubuntu-18.04/
# 
# Ensuite, il faut installer les librairies: 
#   numpy 
#   scipy
#   os
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

# creer la class SchemeAnalysisSuite
class SchemeAnalysisSuite():

  # constructeur de la class
  def __init__(self,convergenceParameters={},figureParameters={}):

    # sauvegarder les parametres de convergence
    # Ce dictionaire contient tout les parametres pour
    # effectuer un etude de convergence ou conservation
    # parametres:
    #   fileDelimeter:		    (str) character qui delimite les colonnes
    #					  du fichier d'output
    #   binaryFileName:		    (str) nom de l'executable
    #   inputFileName:		    (str) nom du fichier d'input
    #   numberOfPointForRegression: (int) nombre de points pour effecture
    #					  une regressione lineaire
    #   meshGenerator:		    (str) nombre de points pour generer un
    #					  maillage, available:
    #					  0) input
    #					  1) linear
    #					  2) power10
    #   meshValues:		    (list(float64))(NValues) valeurs pour
    #					  generer le maillage: si input
    #					  est utilise meshValues doit 
    #					  contenir le maillage entier
    #					  autrement seulement les valeurs
    #					  minimale et maximale sont demandees
    #					  le nombre de points est facultatif
    #   keyForConvergence:	    (str) clee du parametre de convergence
    #					  (nom du parametre utilise pour
    #					  l'etude de convergence)
    #   keyForOutputFileName:	    (str) clee du parametre contenant le nom
    #					  du fichier d'output
    #   interpolationValues:	    (list(float64)) abscisse(s) d'interpolation
    #					  pour tous les valeurs
    #   interpolationIndexes:	    (list(int)) indexe(s) des abscisse(s)
    #					  d'interpolation pour tous les valeurs
    #   indexList:		    (list(int))(NIndexes) liste des indexes des
    #					  ordonees utilisees a etudier
    #   analyticalSolution:	    (list(float64)) liste des solutions analytique
    #					  avec lesquelles calculer l'erreur 
    #					  de convergence
    #   inputParameters:	    (dict) dictionaire avec les parametres
    #					   a ecrire dans le ficher d'input
    self.convergenceParameters = convergenceParameters

    # sauvegarder les parametres pour faire les figures
    # attention: la methode de doConvergenceTest demande
    # un figure avec NIndexes subplots. Les parametres
    # de la figure de convergence doit etre le premier
    # dictionaire dans figureParameters['plotData']
    # pour plus d'information: regarder le module PlotResults
    self.figureParameters = figureParameters

    # liste des nomes du ficher d'output
    self._outputFileNameList = []
    # array de valeurs utilisees pour etudier les schema
    # comme par, exemple, les valeur de pas de temps pour
    # l'etude de convergence. Taille: (NMesh) 
    self.__mesh = np.array([])
    # array des resultats. Taille: (NIndexes,NMesh)
    self.__results = np.array([])
    # array des erreurs
    self.__errors = np.array([])
    # array des regressions lineaires des erreurs
    # 0) slopes 1) intersection
    self.__errorsLinearRegression = \
    np.array((2,len(self.convergenceParameters['indexList'])),dtype=np.float64)
    # liste contenant les coordonnees des nodes d'interpolation
    # [[coordonnees node 1],[coordonnees node 2],...]
    self.__interpolationNodeCoordinates = []
    # liste contenant les indexes des nodes d'interpolation
    # [[indexes node 1],[indexes node 2],...]
    self.__interpolationNodeIndexes = []


# ----------------------------------------------------------------------------------------------

  # destructeur de la class
  def __del__(self):

      # imprimer le nome de la class
      print(self.__class__.__name__ ,'is deleted')

# ----------------------------------------------------------------------------------------------

  # cette methode retourne la mesh
  # output:
  #   mesh: (array)(NMesh) array contenant le maillage
  def getMesh(self):
  
    # retourner le maillage
    return self.__mesh

# ----------------------------------------------------------------------------------------------

  # cette methode permets de changer le maillage
  # inputs:
  #   mesh: (array)(NMesh) array contenant le maillage
  def setMesh(self,mesh):

    # copier le maillage
    self.__mesh = np.copy(mesh)

# ----------------------------------------------------------------------------------------------

  # cette methode retourne la resultats
  # output:
  #   results: (array)(NIndexes,NMesh) array contenant les resultats
  def getMesh(self):
  
    # retourner le maillage
    return self.__results

# ----------------------------------------------------------------------------------------------

  # Cette methode construit une maillage
  def generateMesh(self):

    # verifier la procedure pour genere une mesh
    if(self.convergenceParameters['meshGenerator']=='input'):
      # copier le maillage
      self.__mesh = np.copy(self.convergenceParameters['meshValues'])
    else:
      # verifier si le nombre de points est donne
      if(len(self.convergenceParameters['meshValues'])==2):
        # ajuter le nombre de points
        self.convergenceParameters['meshValues'].append(\
        self.convergenceParameters['meshValues'][1]-\
        self.convergenceParameters['meshValues'][0]+1)
        # calculer le maillage lineaire
        self.__mesh = np.linspace(self.convergenceParameters['meshValues'][0],\
        self.convergenceParameters['meshValues'][1],\
        num=self.convergenceParameters['meshValues'][2],\
        endpoint=True,dtype=np.float64)
        # verifier le type de maillage
        if(self.convergenceParameters['meshGenerator']=='power10'):
          # construire un maillage exposant10
          self.__mesh = np.power(10.0,self.__mesh,dtype=np.float64)

# ----------------------------------------------------------------------------------------------

  # cette methode permet de executer le binaire et 
  # lire les donnees.
  # inputs:
  #   meshValue: (value) valeur de l'element du maillage
  #			 a ecrire dans le fichier d'input
  def runExecutableAndReadData(self,meshValue):

    # importer la librerie os
    import os
  
    # surecrire la valeur de convergence pour le fichier # d'input
    self.convergenceParameters['inputParameters'][\
    self.convergenceParameters['keyForConvergence']]=meshValue
    # generer un nome pour le fichier d'output
    self.convergenceParameters['inputParameters'][\
    self.convergenceParameters['keyForOutputFileName']]=\
    ''.join(['output_',\
    self.convergenceParameters['keyForConvergence'],str(meshValue),'.out'])
    # ajouter le nome du fichier d'output a la list
    self._outputFileNameList.append(\
    self.convergenceParameters['inputParameters'][\
    self.convergenceParameters['keyForOutputFileName']])
    # open input file
    with open(self.convergenceParameters['inputFileName'],'w') as inputFile:
      # boucle sur les clees du fichier d'input
      for key,value in self.convergenceParameters['inputParameters'].items():
        # ecrire valeurs dans le fichier d'input
        inputFile.write(''.join([key,' = ',str(value),'\n']))
    # executer le binaire
    os.system(''.join(['./',self.convergenceParameters['binaryFileName'],\
    ' ',self.convergenceParameters['inputFileName']]))
    # lire les resultats
    self.__results = np.loadtxt(\
    self.convergenceParameters['inputParameters'][\
    self.convergenceParameters['keyForOutputFileName']],\
    dtype=np.float64,delimiter=self.convergenceParameters['fileDelimiter'])  

# ----------------------------------------------------------------------------------------------

  # cette methode selectionne la methode pour chercher les nodes
  # d'interpolation et celle d'interpolation en fonction du nombre 
  # de coordonnees.
  # outputs:
  #   interpolator: (func) methode d'interpolation
  #   nodeFinder:   (func) methode pour trouver les nodes d'interpolation
  def selectInterpolationAndNodeFinderMethods(self):

    interpolator=None # initialiser la methode d'interpolation
    nodeFinder=None   # initialiser la methode pour trouver les nodes
    # verifier la taille de la list des index d'interpolation
    if(len(self.convergenceParameters['interpolationIndexes'])==1):
      # copier la methode d'interpolation
      interpolator = self.interpolator1D
      # copier la methode pour trouver les nodes d'interpolation
      nodeFinder = self.findInterpolationNodes1D
    else:
      # imprimer un message d'erreur
      print('Error: interpolation scheme not implemented')

    # retourner les deux methodes
    return interpolator,nodeFinder

# ----------------------------------------------------------------------------------------------

  # cette methode permet d'interpoler valeurs en 1D
  # inputs: 
  #   valueIndex: (int) index du valeur qu'il faut interpoler
  # outputs:
  #   value: (float64) interpolated value
  def interpolator1D(self,valueIndex):

    # effectuer l'interpolation
    return np.interp(self.convergenceParameters['interpolationValues'][0],\
    [self.__interpolationNodeCoordinates[0][0],self.__interpolationNodeCoordinates[1][0]],
    [self.__results[self.__interpolationNodeIndexes[0][0],valueIndex],\
    self.__results[self.__interpolationNodeIndexes[1][0],valueIndex]])

# ----------------------------------------------------------------------------------------------

  # cette methode recherche le point du maillage (1D) 
  # dans l'array des solutions et interpole linearment
  # la solution. Si la valeur du maillage n'est pas
  # trouve le dernier resultat est retourne
  def findInterpolationNodes1D(self):

    # chercher l'index le plus proche du valeur donnee
    nearestNodeIndex = (np.abs(\
    self.__results[:,self.convergenceParameters['interpolationIndexes'][0]]-\
    self.convergenceParameters['interpolationValues'][0])).argmin()
    # extraire le valeur du node le plus proche
    nearestInterpolationValue = self.__results[\
    nearestNodeIndex,self.convergenceParameters['interpolationIndexes'][0]]
    # verifier s'il est le node de gauche
    if(nearestInterpolationValue<=self.convergenceParameters['interpolationValues'][0]):
      # copier l'nearestNodeIndex comme premier node et nearestNodeIndex+1 comme deuxieme
      self.__interpolationNodeIndexes=[[nearestNodeIndex],\
      [np.amin([nearestNodeIndex+1,self.__results.shape[0]-1])]]
    else:
      # copier l'nearestNodeIndex-1 comme premier node et nearestNodeIndex comme deuxieme
      self.__interpolationNodeIndexes=[[np.amax([nearestNodeIndex-1,0])],\
      [nearestNodeIndex]]
    # boucle sur les valeurs au nodes
    for nodeCoordinates in self.__interpolationNodeIndexes:
      # sauvegarder les valeurs au nodes
      self.__interpolationNodeCoordinates.append([\
      self.__results[nodeCoordinates[0],self.convergenceParameters['interpolationIndexes'][0]]])

# ----------------------------------------------------------------------------------------------

  # cette methode retourne la fonction de 
  # regressione lineaire desideree en fonction
  # du type du maillage
  # output:
  #   regression: (funct) regression lineaire
  def selectLinearRegression(self):

    # importer la librerie avec la regression lineaire 
    import scipy.stats.mstats as mstats

    # select la regression lineare
    if(self.convergenceParameters['meshGenerator']=='power10'):
      # utiliser la regression du logarithm base 10
      return self.linearRegressionLog10
    else:
      # utiliser la regression lineaire standard
      return mstats.linregress

# ----------------------------------------------------------------------------------------------

  # cette methode overload la regression lineare avec
  # celui du logarithme
  # input:
  #   a: (value) premier valeur
  #   b: (value) deuxieme valeur
  # output:
  #   meme que la refression lineaire
  def linearRegressionLog10(self,a,b):

    # importer la librerie avec la regression lineaire 
    import scipy.stats.mstats as mstats

    # retourner la regression base 10
    return mstats.linregress(np.log10(np.abs(a)),np.log10(np.abs(b)))

# ----------------------------------------------------------------------------------------------

  # cette methode execute l'etude de convergence pour une serie de valeurs
  def doConvergenceTest(self):

    # importer le module pour faire des figures
    import PlotResults
  
    # calculer le maillage pour le test de convergence
    self.generateMesh()
    # initialiser l'array des erreurs
    self.__errors = np.zeros((self.__mesh.shape[0],\
    len(self.convergenceParameters['indexList'])),dtype=np.float64)
    # select linear regression method
    linearRegression = self.selectLinearRegression()

    # initialiser les methodes d'interpolation
    # et pour pour trouver les nodes d'interpolation
    interpolator,findNode = self.selectInterpolationAndNodeFinderMethods()

    # boucle sur les elements du maillage
    for elementId,element in enumerate(self.__mesh):
      # executer une nouvelle simulation
      self.runExecutableAndReadData(element)
      # extraire les nodes d'interpolation
      findNode() 
      # boucle sur les indexes des valeurs de convergence
      for indexId,index in enumerate(self.convergenceParameters['indexList']):
        # calculer l'erreur en interpolant la solution
        self.__errors[elementId,indexId] = np.abs(interpolator(index)-\
        self.convergenceParameters['analyticalSolution'][indexId])

    # boucle sur les index
    for indexId,index in enumerate(self.convergenceParameters['indexList']):
      # calculer la regression lineaire sur les derniers x-points
      slope,intersection,r_value,p_value,std_err = \
      linearRegression(self.__mesh[self.__mesh.shape[0]-\
      self.convergenceParameters['numberOfPointForRegression']:-1],\
      self.__errors[self.__mesh.shape[0]-\
      self.convergenceParameters['numberOfPointForRegression']:-1,indexId])
      # imprimer la pente
      print('Convergence test index: ',index,' slope: ',slope)
      # verifier la presence de la figure pour la convergence
      if(len(self.figureParameters['plotData'])!=0):
        # sauvegarder la pente comme legende
        self.figureParameters['plotData'][0]['legend'].append(\
        [''.join(['slope: ',str(np.around(slope,decimals=2))])])

    # verifier la presence de la figure pour la convergence
    if(len(self.figureParameters['plotData'])!=0):
      # creer le dictionaire des parametres
      parameters = self.figureParameters.copy()
      # remove all other plot data
      parameters['plotData']=[self.figureParameters['plotData'][0]]
      # initialiser array pour plot
      plotValues = np.zeros((self.__errors.shape[0],\
      self.__errors.shape[1]+1),dtype=np.float64)
      # copier la mesh
      plotValues[:,0] = self.__mesh
      # coper les erreurs
      plotValues[:,1:] = self.__errors
      # initialiser un objet pour generer des figures
      plotResult = PlotResults.PlotResults(values=plotValues,\
      figureParameters=parameters)
      # fare la figure
      plotResult.SimplePlotFigures()
  
# ----------------------------------------------------------------------------------------------

