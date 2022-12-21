#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
                          // Fichier .tpp car inclut fonctions template
#include <numeric>
using namespace std; // ouvrir un namespace avec la librerie c++ de base


template<typename T> T scalarProduct(valarray<T> const& array1,
valarray<T> const& array2){
  // compute and return the norm2 of a valarray
  return (array1*array2).sum();
}

template<typename T> T norm2(valarray<T> const& array){
  // compute and return the norm2 of a valarray
  return sqrt((array*array).sum());
} 

/* definir a fonction template pour calculer le produit vecteur
   entre 2 valarray de dimension 3
   inputs:
     array1, array2: (valarray<T>)(N) vecteurs de taille N
   outputs:
     produitVecteur: (T) produit vectoriel array1 x aray2 
*/
template<typename T> valarray<T> produitVecteur(valarray<T> const& array1,\
valarray<T> const& array2){

  valarray<T> array3=valarray<T>(3);

  array3[0] = array1[1]*array2[2] - array1[2]*array2[1]; // premier composante
  array3[1] = array1[2]*array2[0] - array1[0]*array2[2]; // deuxieme composante
  array3[2] = array1[0]*array2[1] - array1[1]*array2[0]; // troisieme composante

  return array3;
} 


class Engine
{

protected:
  // definition des constantes
  const double pi = 3.1415926535897932384626433832795028841971e0;
  const double G = 6.674*pow(10,-11);
  const double g = 9.81;
  // definition des variables
  double gamma; 		// parametre gamma 
  unsigned int nsteps;  // Nombre de pas
  double rho_0; 	 	// parametre rho à l'altitude 0
  double P_0;			// pression au sol
  double epsilon;		// décalage pour l'intégration en descente
  double z0;			// position initiale  
  unsigned int sampling;// Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;    // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile; // Pointeur vers le fichier de sortie

  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
     write: (bool) ecriture de tous les sampling si faux
  */  
  

  void printOut(bool write)
  { 
    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      *outputFile << z << ' ' << rho << endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }
  
  // Iteration temporelle, a definir au niveau des classes filles
  virtual void step()=0;


  // donnes internes
  double z, dz;
  double rho;

  double drho(double r) const
  {
	  double K = P_0*pow(rho_0, -gamma);
	  //cout << (-g/(K*gamma))*pow(r,2-gamma) << endl;
	  return (-g/(K*gamma))*pow(r,2-gamma);
	  
  }
  
double rho_ana(double z)
{
	double K = P_0*pow(rho_0, -gamma);
	return pow( pow(rho_0,gamma-1) - g*z*(gamma-1) / (K*gamma), 1/(gamma-1) );
}
public:

  Engine(ConfigFile configFile)
  {
    // variable locale
    
    // Stockage des parametres de simulation dans les attributs de la classe
    dz       = configFile.get<double>("dz",dz);
    nsteps   = configFile.get<double>("nsteps",nsteps);
    rho_0    = configFile.get<double>("rho_0",rho_0);
    P_0      = configFile.get<double>("P_0",P_0);
	gamma	 = configFile.get<double>("gamma",gamma);
	epsilon  = configFile.get<double>("epsilon",epsilon);
    sampling = configFile.get<unsigned int>("sampling",sampling); // lire le parametre de sampling 

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output","output.out").c_str()); 
    outputFile->precision(14); // Les nombres seront ecrits avec 15 decimales
  };

  // Destructeur virtuel
  virtual ~Engine()
  {
    outputFile->close();
    delete outputFile;
  };

  // Simulation complete
  virtual void run() = 0;

}; 




// Extension de la class Engine implementant l'integrateur Runge-Kutta 4
class Enginemonter: public Engine
{
public:
  Enginemonter(ConfigFile configFile): Engine(configFile) {}


  void step()
  {
	  double k1(0.0),k2(0.0),k3(0.0),k4(0.0);
	  k1 = dz * drho(rho);
	  k2 = dz * drho(rho + 0.5 * k1);
	  k3 = dz * drho(rho + 0.5 * k2);
	  k4 = dz * drho(rho + k3);
	  rho += (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
    z+=dz;
    
  }
  virtual void run() override
  {
    z = 0;   // initialiser la position
    rho = rho_0;//initialiser rho
    last = 0; // initialise le parametre d'ecriture
    printOut(true); // ecrire premier pas de temps

    while(rho > rho_0*1e-6)
    {
      step();  // faire la mise a jour de la simulation 
      printOut(false); // ecrire pas de temps actuel
    }
    printOut(true); // ecrire dernier pas de temps
  }
  
};

class Enginedescendre: public Engine
{
public:
  Enginedescendre(ConfigFile configFile): Engine(configFile) {}


  void step()
  {
	  double k1(0.0),k2(0.0),k3(0.0),k4(0.0);
	  
	  k1 = -dz * drho(rho);
	  k2 = -dz * drho(rho + 0.5 * k1);
	  k3 = -dz * drho(rho + 0.5 * k2);
	  k4 = -dz * drho(rho + k3);
	  rho += (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
	  z -= dz;
  }
    virtual void run() override
  {
    z = (gamma*P_0/((gamma-1) * rho_0 * g)) - epsilon ;   // initialiser la position
    dz = z/nsteps;
    rho = rho_ana(z);//initialiser rho
    last = 0; // initialise le parametre d'ecriture
    printOut(true); // ecrire premier pas de temps

    while( z > 1.0e-6 )
    {
      step();  // faire la mise a jour de la simulation 
      printOut(false); // ecrire pas de temps actuel
    }
    printOut(true); // ecrire dernier pas de temps
  }
};


// programme
int main(int argc, char* argv[])
{
  string inputPath("configuration.in.example"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);


  string schema(configFile.get<string>("schema"));

  Engine* engine; // definer la class pour la simulation
  // choisir quel schema utiliser
  if(schema == "monter" || schema == "M")
  {
    // initialiser une simulation avec schema Euler
    engine = new Enginemonter(configFile);
  }
  else if(schema == "descendre" || schema == "D")
  {
    engine = new Enginedescendre(configFile);
  }
  else
  {
    cerr << "Schema inconnu" << endl;
    return -1;
  }

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation 
  cout << "Fin de la simulation." << endl;
  return 0;
}
