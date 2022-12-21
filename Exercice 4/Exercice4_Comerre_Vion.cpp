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
  const double G = 6.674e-11;
  const double g = 9.81;
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
  
  valarray<double> x0=valarray<double>(0.e0,3); // vecteur contenant la position initiale 
  valarray<double> v0=valarray<double>(0.e0,3); // vecteur contenant la vitesse initiale 
  
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
	valarray<double> acc =valarray<double>(0.e0,3);
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
  valarray<double> x = valarray<double>(2);
  valarray<double> v = valarray<double>(2);

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
  
  double rho_ana(double z)
  {
    if (abs(P_0) < 1e-10 || abs(rho_0) < 1e-10)
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
    valarray<double> a =valarray<double>(0.e0,3); //TODO écrire RK4 avec rho qui évolue aussi
    valarray<double> k1_x,k1_v,k2_x,k2_v,k3_x,k3_v,k4_x,k4_v (0.e0,3);
    valarray<double> x_ = x; // position
    valarray<double> v_ = v; // vitesse
    // double rho__ = rho_ana(x);     // rho, réfléchir si on intègre dans l'autre sens, il peut ne pas être ici si c'est pas possible
    acceleration(a);
	k1_x = dt * v;
	k1_v = dt * a;

	x = x_ + k1_x/2.0;
	v = v_ + 0.5 * k1_v;
  rho = rho_ana(norm2(x));
  acceleration(a); 
  k2_x = dt * v;
  k2_v = dt * a;

  x = x_ + 0.5 * k2_x;
  v = v_ + 0.5 * k2_v;
  rho = rho_ana(norm2(x));
  acceleration(a);
  k3_x = dt * v;
  k3_v = dt * a;

  x = x_ + k3_x;
  v += v_ + k3_v;
  rho = rho_ana(norm2(x));
  acceleration(a);
  k4_x = dt * v;
  k4_v = dt * a;

  x = x_ + (1.0/6.0) * (k1_x + 2.0*k2_x + 2.0*k3_x + k4_x);
  v = v_ + (1.0/6.0) * (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v);
  rho = rho_ana(norm2(x));
  acceleration(a);
  t=t+dt;
  // x = x_ + dt * v;
  // v = v_ + dt * a;
  // t=t+dt;
  // acceleration(a);
  }
};

class EngineRungeKutta4adaptatif: public Engine
{
public:
  EngineRungeKutta4adaptatif(ConfigFile configFile): Engine(configFile) {}

  void rk4(valarray<double>& x_1, valarray<double>& v_1, double t_1, double dt_1, double rho_1, valarray<double>& a_1, double div=1.0)
  {
  valarray<double> a = a_1; //TODO écrire RK4 avec rho qui évolue aussi
  valarray<double> k1_x,k1_v,k2_x,k2_v,k3_x,k3_v,k4_x,k4_v (0.e0,3);
  valarray<double> x_ = x_1; // position
  valarray<double> v_ = v_1; // vitesse
  dt_1 = dt_1 / div;
	k1_x = dt_1 * v_1;
	k1_v = dt_1 * a_1;

	x_1 = x_ + 0.5 * k1_x;
	v_1 = v_ + 0.5 * k1_v;
  rho_1 = rho_ana(norm2(x_1));
  acceleration(a); 
  k2_x = dt_1 * v_1;
  k2_v = dt_1 * a_1;

  x_1 = x_ + 0.5 * k2_x;
  v_1 = v_ + 0.5 * k2_v;
  rho_1 = rho_ana(norm2(x_1));
  acceleration(a);
  k3_x = dt_1 * v_1;
  k3_v = dt_1 * a;

  x_1 = x_ + k3_x;
  v_1 = v_ + k3_v;
  rho_1 = rho_ana(norm2(x_1));
  acceleration(a);
  k4_x = dt_1 * v_1;
  k4_v = dt_1 * a_1;

  x_1 = x_ + (1.0/6.0) * (k1_x + 2.0*k2_x + 2.0*k3_x + k4_x);
  v_1 = v_ + (1.0/6.0) * (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v);
  rho_1 = rho_ana(norm2(x_1));
  acceleration(a_1);
  }

  void step()
  {
    valarray<double> a =valarray<double>(0.e0,3); //TODO écrire RK4 avec pas de temps adaptatif
    
    valarray<double> x_0 = x; // position 0
    valarray<double> x_1 = x; // position 1
    
    valarray<double> v_0 = v; // vitesse 0
    valarray<double> v_1 = v; // vitesse 1
    
    valarray<double> a_0 = a; // accélération 0
    valarray<double> a_1 = a; // accélération 1

    double dt_1 = dt;
    

    //rungekutta 4 with adaptative time step
    rk4(x_0,v_0,t,dt_1,rho,a_0);
    rk4(x_1,v_1,t,dt_1,rho,a_1,2.0);
    rk4(x_1,v_1,t,dt_1,rho,a_1,2.0);
    valarray<double> x_diff = x_1-x_0;
    double diff = norm2(x_diff);
    double tol = 1.e0;
    unsigned int i = 0;
    while (diff>tol)
    {
      x_0 = x;
      v_0 = v;
      a_0 = a;
      x_1 = x;
      v_1 = v;
      a_1 = a;
      dt_1 = dt_1 * 0.999 *pow((tol/diff),(1.0/(1.0+i)));
      rk4(x_0,v_0,t,dt_1,rho,a_0);
      rk4(x_1,v_1,t,dt_1,rho,a_1,2.0);
      rk4(x_1,v_1,t,dt_1,rho,a_1,2.0);
      x_diff = x_1-x_0;
      diff = norm2(x_diff);
      ++i;
    }
    x = x_1;
    v = v_1;
    t = t + dt_1;
    rho = rho_ana(norm2(x));
    acceleration(a);
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
