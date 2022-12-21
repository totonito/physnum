#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <algorithm>      // librerie algorithmes de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
                          // Fichier .tpp car inclut fonctions template
#include <numeric>
using namespace std; // ouvrir un namespace avec la librerie c++ de base
typedef valarray<double> vecteur; // definir un type vecteur


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
  const double G = 6.674e-11;
  const double g = 9.81;
  const double RT = 6378100e0;
  // definition des variables
  double gamma; // si gamma est modifié, la fonction rho doit être changée
  double tfin;          // Temps final
  unsigned int nsteps;  // Nombre de pas de temps
  double rho_0; 	 	// parametre rho à l'altitude 0
  double P_0;			// pression au sol
  double m_t;           // masse de la terre
  double m_a;			// masse du vaisseau Appolo
  double C_x;			// coefficient de trainee
  double d;				// diametre du vaisseau
  // double K = P_0*pow(rho_0,-1*gamma);
  double dz;
  
  valarray<double> x0 = valarray<double>(0.e0,2); // vecteur contenant la position initiale 
  valarray<double> v0 = valarray<double>(0.e0,2); // vecteur contenant la vitesse initiale 
  
  unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
     write: (bool) ecriture de tous les sampling si faux
  */  
  

  void printOut(bool write)
  { 
    
	double Emec_ = Emec();
	double Pnc_ = Pnc ();
  double rho_ = rho_ana(norm2(x));

	valarray<double> acc =valarray<double>(0.e0,2);
	acceleration(acc);
	double norme_acc = norm2(acc);
	
    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      *outputFile << t << " " << x[0] << " " << x[1] << " " \
      << v[0] << " " << v[1] << " " \
      << Emec_ << " " << Pnc_ << " " << rho_ << " " \
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

  valarray<double> x = valarray<double>(0.e0,2);
  valarray<double> v = valarray<double>(0.e0,2);

  double rho; // sujet a modifications
  // double z(sqrt(pow(x[0],2.0) + pow(y[0],2))-6378100); // sujet a modifications
  	  

  
  // TODO
  /* Calcul de l'acceleration totale
     output:
       a: (valarray<double>)(3) vecteur acceleration
  */
  
    double Emec() const
  {
	return 0.5 * m_a * norm2(v) * norm2(v) - ( m_a * m_t * G) / norm2 (x);
  }
  
  double Pnc() const
  {
	  return (m_a * m_t * G) * (x[0]*v[0] + x[1]*v[1])/pow(scalarProduct(x,x),1.5);
  }
  

  double rho_ana(const double& z)
  {
    if (abs(P_0) < 1e-10 || abs(rho_0) < 1e-10 || z-RT > 600000)
      return 0.0;
    else
    {
      double K = P_0*pow(rho_0, -gamma);
      return pow(pow(rho_0,gamma-1.0) - g*z*(gamma-1.0) / (K*gamma), 1.0/(gamma-1.0)); 
    }
  }

  void acceleration(valarray<double>& a) const
  { 
    a[0]      = (-1.0/8.0) * C_x * norm2(v) * pi * d * d * rho * v[0] / m_a  \
                 - (G * m_t * x[0]) / pow(norm2(x),3.0) ; // composante x acceleration
    a[1]      = (-1.0/8.0) * C_x * norm2(v) * pi * d * d * rho * v[1] / m_a  \
                 - (G * m_t * x[1]) / pow(norm2(x),3.0) ;  // composante y acceleration

  }

  valarray<double> acceleration_rk4(const valarray<double>& x_1, const valarray<double>& v_1)
  {
    valarray<double> a = valarray<double>(0.e0,3);
    a[0]      = (-1.0/8.0) * C_x * norm2(v_1) * pi * d * d * rho_ana(norm2(x_1)) * v_1[0] / m_a  \
                 - (G * m_t * x_1[0]) / pow(norm2(x_1),3.0) ; // composante x acceleration
    a[1]      = (-1.0/8.0) * C_x * norm2(v_1) * pi * d * d * rho_ana(norm2(x_1)) * v_1[1] / m_a  \
                 - (G * m_t * x_1[1]) / pow(norm2(x_1),3.0) ;  // composante y acceleration
    return a;

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
	  m_t		 = configFile.get<double>("m_t",m_t);
	  m_a		 = configFile.get<double>("m_a",m_a);
	  C_x		 = configFile.get<double>("C_x",C_x);
	  d		 = configFile.get<double>("d",d);
	  gamma	 = configFile.get<double>("gamma",gamma);
    sampling = configFile.get<unsigned int>("sampling",sampling); // lire le parametre de sampling
    dz = configFile.get<double>("dz",dz);
    string schema(configFile.get<string>("schema"));


  if(schema == "RK4")
    dt = tfin / nsteps;
    else if (schema == "RK4A")
    {
      dt = configFile.get<double>("dt",dt);
      nsteps = tfin / dt;
    }
    

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

    // for(unsigned int i(0); i<nsteps; ++i) // boucle sur les pas de temps
    // {
    //   step();  // faire la mise a jour de la simulation 
    //   printOut(false); // ecrire pas de temps actuel
    // }

    while (t < tfin) // boucle sur le temps final
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

  void step() // rk4

  {
    valarray<double> x_ = x; // position
    valarray<double> v_ = v; // vitesse
    valarray<double> k1x = dt * v;
    valarray<double> k1v = dt * acceleration_rk4(x + k1x * 0.5, v + 0.5 * acceleration_rk4(x, v));
    valarray<double> k2x = dt * (v + k1v * 0.5);
    valarray<double> k2v = dt * acceleration_rk4(x + k2x * 0.5, v + k1v * 0.5);
    valarray<double> k3x = dt * (v + k2v * 0.5);
    valarray<double> k3v = dt * acceleration_rk4(x + k3x, v + k2v * 0.5);
    valarray<double> k4x = dt * (v + k3v);
    valarray<double> k4v = dt * acceleration_rk4(x + k4x, v + k3v);

    // Update the position and velocity using the RK4 formula
    x = x_ + (k1x + k2x * 2.0 + k3x * 2.0 + k4x) / 6.0;
    v = v_ + (k1v + k2v * 2.0 + k3v * 2.0 + k4v) / 6.0;
    t += dt; 

  }
};

class EngineRungeKutta4adaptatif: public Engine
{
public:
  EngineRungeKutta4adaptatif(ConfigFile configFile): Engine(configFile) {}


  valarray<double> rk4(const valarray<double>& state, const double& dt_1, const double& div = 1.0)
  {
    // Split the state vector into position and velocity
    valarray<double> x_2 = state[slice(0, 2, 1)];
    valarray<double> v_2 = state[slice(2, 2, 1)];
    double dt_2 = dt_1 / div;

    // Compute the intermediate values for the Runge-Kutta algorithm
    valarray<double> k1x = dt_2 * v_2;
    valarray<double> k1v = dt_2 * acceleration_rk4(x_2 + k1x * 0.5, v + 0.5 * acceleration_rk4(x_2, v_2));
    valarray<double> k2x = dt_2 * (v_2 + k1v * 0.5);
    valarray<double> k2v = dt_2 * acceleration_rk4(x_2 + k2x * 0.5, v_2 + k1v * 0.5);
    valarray<double> k3x = dt_2 * (v_2 + k2v * 0.5);
    valarray<double> k3v = dt_2 * acceleration_rk4(x_2 + k3x, v_2 + k2v * 0.5);
    valarray<double> k4x = dt_2 * (v_2 + k3v);
    valarray<double> k4v = dt_2 * acceleration_rk4(x_2 + k4x, v_2 + k3v);

    x_2 += (k1x + k2x * 2.0 + k3x * 2.0 + k4x) / 6.0;
    v_2 += (k1v + k2v * 2.0 + k3v * 2.0 + k4v) / 6.0;

    // Concatenate the updated position and velocity into the state vector
    valarray<double> retour(4);
    retour[slice(0, 2, 1)] = x_2;
    retour[slice(2, 2, 1)] = v_2;
    return retour;
  }

  void step()
  {
    
    valarray<double> a =valarray<double>(0.e0,3); //TODO écrire RK4 avec pas de temps adaptatif
    // Initializing the state vectors
    valarray<double> state_0(4);
    state_0[slice(0, 2, 1)] = x;
    state_0[slice(2, 2, 1)] = v;

    valarray<double> state_1(4);
    state_1[slice(0, 2, 1)] = x;
    state_1[slice(2, 2, 1)] = v;
    
    // Initializing the time step
    double dt_1 = 0.0;
    double t_fin = 259200.0;
    dt_1 = min(dt, t_fin - t);

    // Initializing the tolerance
    double tol = 1.e-6;
    unsigned int i = 0;
    
    vecteur y_0 = rk4(state_0,dt_1);
    vecteur y_1 = rk4(state_1,dt_1, 2.0);
    vecteur y_1tilde = rk4(y_1,dt_1, 2.0);

    valarray<double> x_0 = y_0[slice(0, 2, 1)];
    valarray<double> x_1 = y_1tilde[slice(0, 2, 1)];
    valarray<double> x_diff = x_1 - x_0;
    double diff = norm2(x_diff);
    if(diff>tol)
    {
      while (diff>tol)
      {
        dt_1 = dt_1 * 0.999 * pow((tol/diff),(1.0/5.0));
        std::cout << "dt_1 = " << dt_1 << std::endl;

        y_0 = rk4(state_0,dt_1);
        y_1 = rk4(state_1,dt_1, 2.0);
        y_1tilde = rk4(y_1,dt_1, 2.0);
        x_0 = state_0[slice(0, 2, 1)];
        x_1 = state_1[slice(0, 2, 1)];
        x_diff = x_1 - x_0;
        diff = sqrt(x_diff[0]*x_diff[0] + x_diff[1]*x_diff[1]);
      }
      x = y_1tilde[slice(0, 2, 1)];
      v = y_1tilde[slice(2, 2, 1)];
      t = t + dt_1;
      ++i;
    }
    else
    {
      x = y_1tilde[slice(0, 2, 1)];
      v = y_1tilde[slice(2, 2, 1)];
      t = t + dt_1;
      ++i;
    }
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
