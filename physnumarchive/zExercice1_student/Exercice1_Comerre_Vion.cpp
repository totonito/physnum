#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "ConfigFile.h" // Fichier .tpp car inclut un template

using namespace std;

class Engine
{

private:
 
  // definition des variables
  double t,dt; //<- temps  et pas de temps
  double x; // position de la fusee
  double v; // vitesse de la fusee
  //Parametres physiques
  double mF;     // mass de le carburant[kg]
  double D_f;  // debit massique de propulseurs [kg/s]
  double v_e ;// vitesse  d’ ejection [m/s]
  double tfin   ; // temps de combustion  [s]
  double rho_0 ;   // densité de l'air au sol [kg/m3]
  double Cx    ;  // coefficient  aerodynamique[kg/m3]
  double Surf     ;    // surface de section  [m2]
  double lambda ;  //factoure decadiment densite [m]
  double g0  ;  // acceleration de la pesanteur [m/s2]

  unsigned long long int nsteps; // Nombre d'iterations
  unsigned int sampling; // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last; // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile; // Pointeur vers le fichier de sortie

  // Ecriture des diagnostics
  // inputs:
  //   write: (bool) force l'écriture des données si vrai
  void printOut(bool write)
  {
    // Ecriture tous les [sampling] pas de temps ou si 'write' est vrai
    // sequence des sorties: 1) temps, 2) mass 3)position (m) 3) vitesse (m/s)
    // 4) énergie mécanique de la fusée (J) 5)puissance mécanique de la fusée
    if((!write && last>=sampling) || (write && last!=1))
    {
      *outputFile << t << " " << mF << " " << x << " " << v \
       << " " << energy(mF,x,v) << " " << power(mF,x,v) << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }
  // TODO: calculer la force de frottement sur la fusée
  double Fvisc(double xLocal,double vLocal){
    double Fa(0.0);
    Fa = (1.0/2.0) * rho_0 * exp(-xLocal/lambda) * Cx * Surf * vLocal * vLocal;
    return Fa;
  }

  // TODO: calculer la force de propulsion sur la fusée
  double Fprop()
  {
    double Fp;
    Fp = D_f * v_e;
    return Fp;
  }
  
  // TODO: calculer l'acceleration F(xLocal,vLocal)/m de la fusée
  // inputs:
  //   xLocal: (double) position de la fusée
  //   vLocal: (double) vitesse de la fusée
  // outputs:
  //   acceleration: (double) acceleration de la fusée [m/s^2]
  double acceleration(double mLocal, double xLocal,double vLocal)
  {
    double a(0.0);
    a = -g0 + (1.0 / mF) * ( -Fvisc(xLocal,vLocal) + Fprop() );
    return a;
  }

  // TODO: calculer l'energie mecanique de la fusée
  // inputs:
  //   xLocal: (double) position de la fusée
  //   vLocal: (double) vitesse de la fusée
  // outputs:
  //   energy: (double) energie mecanique de la fusée [J]
  double energy(double mLocal, double xLocal, double vLocal){
    double kinetic(0.0);
    double potential(0.0);
    potential = mLocal * g0 * xLocal;
    kinetic = 0.5 * mLocal * vLocal * vLocal;
    return kinetic + potential;
  }

  // TODO: calculer la puissance mecanique de la fusée
  // inputs:
  //   mLocal: (double) masse de le carburant
  //   xLocal: (double) position de la fusée
  //   vLocal: (double) vitesse de la fusée
  // outputs:
  //   energy: (double) puissance mecanique de la fusée [W]
  double power(double mLocal, double xLocal, double vLocal){
    double p_friction(0.0);
    double p_propulsion(0.0);
    double p_massvariation(0.0);
    p_friction =  - vLocal * Fvisc(xLocal,vLocal);
    p_propulsion = Fprop() * vLocal;
    p_massvariation = - D_f * ( g0 * xLocal + 0.5 * vLocal * vLocal);
        return p_propulsion + p_friction +  p_massvariation;
  }

  // Iteration temporelle
  void step()
  {
    // TODO: Mettre a jour (mF,x,v) avec le schema d'Euler
    // Utiliser la méthode acceleration(x,v) définie ci-dessus
    
    double vdot = acceleration(mF,x,v);

    mF -= dt * D_f ;   // mettre a jour la mass
    x  += dt * v;   // mettre a jour la position
    v  += vdot * dt;   // mettre a jour la vitesse

  }


public:

  // Constructeur
  Engine(int argc, char* argv[]) {

    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice1 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice1 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    // Stockage des parametres de simulation dans les attributs de la classe
    tfin     = configFile.get<double>("tfin");                   // lire temps final
    nsteps   = configFile.get<unsigned long long int>("nsteps"); // lire nombre pas de temps
    dt       = tfin / (double(nsteps)); 	                 // calculer pas de temps
    mF       = configFile.get<double>("mF");    	         // lire la mass de la fusée
    v_e      = configFile.get<double>("v_e");	                 
    D_f      = configFile.get<double>("D_f");                     
    rho_0    = configFile.get<double>("rho_0");                     
    Cx       = configFile.get<double>("Cx");                     
    Surf     = configFile.get<double>("Surf");                     
    lambda   = configFile.get<double>("lambda");                     
    g0       = configFile.get<double>("g0");                     
    sampling = configFile.get<unsigned int>("sampling");        

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales
  };

  // Destructeur
  ~Engine()
  {
    outputFile->close(); // fermer le ficher des sorties
    delete outputFile;   // touer la class outputFile
  };

  // Simulation complete
  void run()
  {
    unsigned long long int i=0; // declare and initialise the index
    t = 0.e0; // initialiser temps
    x = 0.e0;   // initialiser position
    v = 0.e0;   // initialiser vitesse
    last = 0; // initialise ecriture
    printOut(true); // ecrire donnees initialies
    //------ boucle sur les pas de temps ----------
    // on arrête la simulation si la fusée sort de l'intervalle entre la terre et la lune
    while(i<nsteps){
      step();  // integrer la dynamique
      t += dt; // mettre a jour le temps
      printOut(false); // false pour imprimer a chaque pas de temps
      i += 1; // increase the index 
    }
    printOut(true); // imprimer le dernier pas de temps
  };

};

// programme
int main(int argc, char* argv[])
{
  Engine engine(argc, argv); // construire la object engine
  engine.run(); // executer la simulation
  cout << "Fin de la simulation." << endl;
  return 0;
}
