#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
                          // Fichier .tpp car inclut fonctions template
#include <numeric>
using namespace std; // ouvrir un namespace avec la librerie c++ de base


template<typename T> T scalarProduct(valarray<T> const& array1,\
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

private:
  // definition des constantes
  const double pi = 3.1415926535897932384626433832795028841971e0;
  const double G = 6.674*pow(10,-11);
  const double g = 9.81;
  // definition des variables
  double tfin;          // Temps final
  unsigned int nsteps;  // Nombre de pas de temps
  double rho_0; 	 	// parametre rho à l'altitude 0
  double P_0;			// pression au sol
  double gamma;			// parametre adiabadicite
  double m_t;           // masse de la terre
  double m_v;			// masse du vaisseau
  double C_x;			// coefficient de trainee
  double d;				// diametre du vaisseau
  valarray<double> x0=valarray<double>(0.e0,3); // vecteur contenant la position initiale 
  valarray<double> v0=valarray<double>(0.e0,3); // vecteur contenant la vitesse initiale 
  
  unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
     write: (bool) ecriture de tous les sampling si faux
  */  
  
  double E_mec(ce que tu as besoin, je pense rien normalement)// TODO calculer l'energie mecanique
  {}
  
  double Pnc(ce que tu as besoin)// TODO calculer la puissance des forces non conservatives
  {}
  
  double rho (ce que tu as besoin)// TODO calculer rho
  {}
  
  void printOut(bool write)
  { 
    
	double E_mec = Emec(ce que tu as besoin);
	double Pnc = Pnc (ce que tu as besoin);
	double rho = rho(ce que tu as besoin);
	valarray<double> acc =valarray<double>(0.e0,3);
	acceleration(acc);
	double norme_acc = norm2(acc);
    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      *outputFile << t << " " << x[0] << " " << x[1] << " " \
      << v[0] << " " << v[1] << " " \
      << E_mec<< " " << Pnc << " " << rho << " " \
      << norme_acc << " " << dt << endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

  // Iteration temporelle, a definir au niveau des classes filles
  virtual void step()=0;




protected:

  // donnes internes
  double t,dt;
  valarray<double> x =valarray<double>(3);
  valarray<double> v =valarray<double>(3);
  
  // TODO
  /* Calcul de l'acceleration totale
     output:
       a: (valarray<double>)(3) vecteur acceleration
  */
  void acceleration(valarray<double>& a) const
  { 
    a[0]      = 0.0; // composante x acceleration
    a[1]      = 0.0; // composante y acceleration
  }


public:

  /* Constructeur de la classe Engine
     inputs:
       configFile: (ConfigFile) handler du fichier d'input
  */
  Engine(ConfigFile configFile)
  {
    // variable locale
    
    // Stockage des parametres de simulation dans les attributs de la classe
    tfin     = configFile.get<double>("tfin",tfin);	       
    nsteps   = configFile.get<unsigned int>("nsteps",nsteps); 
    x0[0]    = configFile.get<double>("x0",x0[0]);	  
    x0[1]    = configFile.get<double>("y0",x0[1]);        
    v0[0]    = configFile.get<double>("vx0",v0[0]);	    
    v0[1]    = configFile.get<double>("vy0",v0[1]);       
    rho_0    = configFile.get<double>("rho_0",rho_0);
    P_0      = configFile.get<double>("P_0",P_0);
	gamma	 = configFile.get<double>("gamma",gamma);
	m_t		 = configFile.get<double>("m_t",m_t);
	m_v		 = configFile.get<double>("m_v",m_v);
	C_x		 = configFile.get<double>("C_x",C_x);
	d		 = configFile.get<double>("d",d);
    sampling = configFile.get<unsigned int>("sampling",sampling); // lire le parametre de sampling

    dt = tfin / nsteps; // calculer le time step

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output","output.out").c_str()); 
    outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales
  };

  // Destructeur virtuel
  virtual ~Engine()
  {
    outputFile->close();
    delete outputFile;
  };

  // Simulation complete
  void run()
  {
    t = 0.e0; // initialiser le temps
    x = x0;   // initialiser la position
    v = v0;   // initialiser la vitesse
    last = 0; // initialise le parametre d'ecriture
    printOut(true); // ecrire premier pas de temps

    for(unsigned int i(0); i<nsteps; ++i) // boucle sur les pas de temps
    {
      step();  // faire la mise a jour de la simulation 
      printOut(false); // ecrire pas de temps actuel
    }
    printOut(true); // ecrire dernier pas de temps
  };

}; 




// Extension de la class Engine implementant l'integrateur Runge-Kutta 4
class EngineRungeKutta4: public Engine
{
public:
  EngineRungeKutta4(ConfigFile configFile): Engine(configFile) {}


  void step()
  {
    valarray<double> a =valarray<double>(0.e0,3); //TODO écrire RK4
    
    valarray<double> x_ = x; // position
    valarray<double> v_ = v; // vitesse

    x+=0.0;
    v+=0.0;

    t+=dt;
    
  }
};

class EngineRungeKutta4adaptatif: public Engine
{
public:
  EngineRungeKutta4adaptatif(ConfigFile configFile): Engine(configFile) {}


  void step()
  {
    valarray<double> a =valarray<double>(0.e0,3); //TODO écrire RK4 avec pas de temps adaptatif
    
    valarray<double> x_ = x; // position
    valarray<double> v_ = v; // vitesse

    x+=0.0;
    v+=0.0;

    t+=dt;
    
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
  if(schema == "RungeKutta4" || schema == "RK4")
  {
    // initialiser une simulation avec schema Euler
    engine = new EngineRungeKutta4(configFile);
  }
  else if(schema == "adaptatif" || schema == "A" || schema == "RK4A")
  {
    engine = new EngineRungeKutta4adaptatif(configFile);
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
