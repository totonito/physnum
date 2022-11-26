#include <iostream>       // basic input output streamsih"iorf uizejbfzfjb zhkfj
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
                          // Fichier .tpp car inclut fonctions template
#include <numeric>
using namespace std; // ouvrir un namespace avec la librerie c++ de base

/* definir a fonction template pour calculer le produit interne
   entre deux valarray
   inputs:
     array1: (valarray<T>)(N) vecteur de taille N
     array2: (valarray<T>)(N) vecteur de taille N
   outputs:
   * 
     produitInterne: (T) produit entre les deux vecteurs
*/ 

template<typename T> T scalarProduct(valarray<T> const& array1,\
valarray<T> const& array2){
  // compute and return the norm2 of a valarray
  return (array1*array2).sum();
} 

template<typename T> T norm2(valarray<T> const& array){
  // compute and return the norm2 of a valarray
  return sqrt((array*array).sum());
} 

template<typename T> T norm2carre(valarray<T> const& array){
  // compute and return the square of the norm to reduce approximation made by the operator sqrt
  return (array*array).sum();
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
  // initialiser le nouveau valarray
  valarray<T> array3=valarray<T>(3);
  // calculer le produit vecteur
  array3[0] = array1[1]*array2[2] - array1[2]*array2[1]; // premier composante
  array3[1] = array1[2]*array2[0] - array1[0]*array2[2]; // deuxieme composante
  array3[2] = array1[0]*array2[1] - array1[1]*array2[0]; // troisieme composante
  // retourner le array3
  return array3;
} 



/* La class Engine est le moteur principale de ce code. Il contient 
   les methodes de base pour lire / initialiser les inputs, 
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{

protected:
  // definition des constantes
  const double pi=3.1415926535897932384626433832795028841971e0;
  // definition des variables
  double tfin;          // Temps final
  unsigned int nsteps;  // Nombre de pas de temps
  double kappa; 	 // parameter kappa
  double mass;          // mass du proton
  double B0;            // parameter B0
  double B1;            // parameter B1
  double B2;            // parameter B2
  double L;             // parameter L
  double q;             // charge du proton 

  valarray<double> E =valarray<double>(0.e0,3); // vecteur contenant le champ ́electrique 
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
  // TODO calculer l'energie mecanique
    double Energy = (0.5) * mass * scalarProduct(v,v)  - q * x[1] * E[1]; 

    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      magnetic_moment(); // moment magnetique du proton
      update_magneticfield();
      valarray<double> vp = (scalarProduct(v,B)/norm2carre(B))*B;
      
      *outputFile << t << " " << x[0] << " " << x[1] << " " \
      << x[2] << " " << v[0] << " " << v[1] << " " << v[2] << " " \
      << mu << " "<< Energy<< " " << scalarProduct(v,B)/norm2(B) << " " << vp[0] <<  " " << vp[1] << " "<< vp[2] << endl; // write output on file
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
  double t,dt;  // Temps courant pas de temps
  double mu; // moment magnetique

  // Ci-dessous, on définit deux vecteurs séparés, x et v, pour le vecteur-position
  // et le vecteur-vitesse, respectivement.
  valarray<double> x =valarray<double>(3); // Position actuelle du proton 
  valarray<double> v =valarray<double>(3); // Vitesse actuelle du proton
  
  // On pourrait aussi définir un seul vecteur y, de taille 6, qui regrouperait (x,v)
  // valarray<double> y =valarray<double>(6); // (Position,vitesse) actuelle du proton

  
  valarray<double> B =valarray<double>(3); // vecteur du champ magnetique 

  
  // TODO
  /* Calcul de l'acceleration totale
     output:
       a: (valarray<double>)(3) vecteur acceleration
  */
  void acceleration(valarray<double>& a) const
  {
    valarray<double> v_cross_B =valarray<double>(0.e0,3); //  
    v_cross_B = produitVecteur(v,B);
    double B_x(0);
    double B_z(0);
	B_x = B2*x[0]*kappa* sin(kappa*x[2]);
	B_z = B0 + B1 * x[0]/L + B2 * cos(kappa * x[2]);
    a[0]      = (q/mass) * v[1] * B_z ; // composante x acceleration
    a[1]      =  ( E[1] + v[2]*B_x - v[0]*B_z )  * (q/mass); // composante y acceleration
    a[2]      =  v[1] * B[0] * (-q/mass); // composante z acceleration
  }

  void update_magneticfield() 
  {
    B[0] = B2*x[0]*kappa* sin(kappa*x[2]);         
    B[1] = 0;
    B[2] = B0 + B1 * x[0]/L + B2 * cos(kappa * x[2]);
  }

  void magnetic_moment() 
  {
    double module_B;
    double v2 = scalarProduct(v,v);
    double v_par;
    update_magneticfield();
    module_B = norm2(B);
    v_par = scalarProduct(v,B)/module_B;
    mu = mass*(v2 - v_par*v_par)/(2*module_B);
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
    tfin     = configFile.get<double>("tfin",tfin);	        // lire la temps totale de simulation
    nsteps   = configFile.get<unsigned int>("nsteps",nsteps); // lire la nombre de pas de temps
    x0[0]    = configFile.get<double>("x0",x0[0]);	    // lire composante x position initiale
    x0[1]    = configFile.get<double>("y0",x0[1]);        // lire composante y position initiale
    x0[2]    = configFile.get<double>("z0",x0[2]);        // lire composante z position initiale
    v0[0]    = configFile.get<double>("vx0",v0[0]);	    // lire composante x vitesse initiale
    v0[1]    = configFile.get<double>("vy0",v0[1]);       // lire composante y vitesse initiale
    v0[2]    = configFile.get<double>("vz0",v0[2]);       // lire composante z vitesse initiale
    B0       = configFile.get<double>("B0",B0);           // lire B0 parameter
    B1       = configFile.get<double>("B1",B1);           // lire B1 parameter
    B2       = configFile.get<double>("B2",B2);           // lire B2 parameter
    L        = configFile.get<double>("L",L);             // lire L parameter
    q        = configFile.get<double>("q",q);             // lire la charge électrique du proton
    mass     = configFile.get<double>("m",mass);          // lire la masse du proton
    kappa    = configFile.get<double>("k",kappa);         // lire le parametre kappa
    E[0]     = configFile.get<double>("Ex",E[0]);         // lire composante x champ ́electrique
    E[1]     = configFile.get<double>("Ey",E[1]);         // lire composante y champ ́electrique
    E[2]     = configFile.get<double>("Ez",E[2]);         // lire composante z champ ́electrique

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

// Extension de la class Engine implementant l'integrateur Boris-Buneman
class EngineBorisBuneman: public Engine
{
	public:
	EngineBorisBuneman(ConfigFile configFile): Engine(configFile) {}
	
	void step()
	{
		double omega(0);
		valarray <double> champ_unitaire = valarray<double>(0.e0,3);
		valarray <double> x_moins = valarray<double>(0.e0,3);
		valarray <double> v_moins = valarray<double>(0.e0,3);
		valarray <double> v_plus = valarray<double>(0.e0,3);
		x_moins = x + v * dt * 0.5;
		v_moins = v + (q/mass) * E * dt * 0.5;
		x = x_moins;
		update_magneticfield();
		omega = q * norm2(B) /mass;
		champ_unitaire =  B / norm2(B); 
		v_plus = v_moins + (omega * dt) / (1 + pow(omega * dt *0.5,2)) * 
		(produitVecteur(v_moins,champ_unitaire) + omega * dt * 0.5 * 
		produitVecteur(produitVecteur(v_moins,champ_unitaire),champ_unitaire)) ;
		v = v_plus + (q/mass)*E * dt *(0.5);
		x = x_moins + v * dt * 0.5;
		t+=dt;
	}
};//fin de la classe BorisBuneman

// Extension de la class Engine implementant l'integrateur Runge-Kutta 2
class EngineRungeKutta2: public Engine
{
public:
  EngineRungeKutta2(ConfigFile configFile): Engine(configFile) {}

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le scheme: Runge-kutta 2
  */
  void step()
  {
    valarray<double> a = valarray<double>(0.e0,3); // 
    // sauvgarder la position et la vitesse au debut de l'intervalle temporel
    valarray<double> x_ = x; // position
    valarray<double> v_ = v; // vitesse
	valarray<double> k_1x (0.e0,3);
	valarray<double> k_1v (0.e0,3);
	valarray<double> k_2x (0.e0,3);
	valarray<double> k_2v (0.e0,3);
	
	acceleration(a);
	k_1x = dt * v;
	k_1v = dt * a;
	t += dt/2;
	x += 0.5 * k_1x;
	v += 0.5 * k_1v;
	acceleration(a);
	k_2x = dt * v ;
	k_2v = dt * a;
	x = x_;
	v = v_;
	x += k_2x;
	v += k_2v;
 	  
    t += dt/2;
    
  }
}; // fin de la classe RungeKutta2



// Extension de la class Engine implementant l'integrateur d'Euler
class EngineEuler: public Engine
{
private:
  unsigned int maxit=1000; // nombre maximale d iterations
  double tol=1.e12;        // tolerance methode iterative
  double alpha=1.0;        // parametre pour le scheme d'Eurler (alpha dans [0,1] -> 1 Esplicit, 0 Implicit, 0.5 semi-implicit)
public:
  
  EngineEuler(ConfigFile configFile): Engine(configFile){
    tol   = configFile.get<double>("tol"); // lire la tolerance pour le method iterative
    maxit = configFile.get<unsigned int>("maxit"); // lire le nombre des iterations maximale 
    alpha = configFile.get<double>("alpha"); // lire le parametre alpha (alpha dans [0,1] -> 1 Eesplicit, 0 Implicit, 0.5 semiimplicit)
  }

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     les schemas: Euler explicite, implicite et semi-implicite
  */
  
  void EulerExplicite(valarray<double>& x, valarray<double>& v , double& t)
  {
	  valarray<double> a = valarray<double>(0.e0,3);
	  acceleration(a);
	  v += a * dt;
	  x += v * dt;
      t += dt;
  }
  
  void EulerImplicite(valarray<double>& x, valarray<double>& v , double& t)
  {
	 valarray<double> x_save = x;
	 valarray<double> v_save = v; 
	 valarray<double> a = valarray<double>(0.e0,3);
	 acceleration(a);
	  valarray<double> a_save = a; 
	 unsigned int it=0;
	 valarray<double> error_x = valarray<double>(0.e0,3);
	 valarray<double> error_v = valarray<double>(0.e0,3);
	 double error = 999e0;
	 
	 while(error > tol and it < maxit)
	 {
		 x = x_save + ( v) * dt;
		 v = v_save + ( a) * dt;
		 update_magneticfield();
		 acceleration(a);
		 it += 1 ;
		 error_x = x - (x_save + ( v) * dt);
		 error_v = v - (v_save + ( a) * dt);
		 error = (sqrt(norm2carre(error_x) + norm2carre(error_v)));
	 }
	 t += dt;
  }
  
  
  
  void EulerSemiImplicite(valarray<double>& x, valarray<double>& v , double& t)
  {
	 valarray<double> x_save = x;
	 valarray<double> v_save = v; 
	 valarray<double> a = valarray<double>(0.e0,3);
	 acceleration(a);
	  valarray<double> a_save = a; 
	 unsigned int it=0;
	 valarray<double> error_x = valarray<double>(0.e0,3);
	 valarray<double> error_v = valarray<double>(0.e0,3);
	 double error = 999e0;
	 
	 while(error > tol and it < maxit)
	 {
		 x = x_save + (v_save * alpha + (1-alpha) * v) * dt;
		 v = v_save + (a_save * alpha + (1-alpha) * a) * dt;
		 update_magneticfield();
		 acceleration(a);
		 it += 1 ;
		 error_x = x - (x_save + (v_save * alpha + (1-alpha) * v) * dt);
		 error_v = v - (v_save + (a_save * alpha + (1-alpha) * a) * dt);
		 error = (sqrt(norm2carre(error_x) + norm2carre(error_v)));
	 }
	 t += dt;
  }
  
   void step()
  {
    EulerSemiImplicite(x,v,t);
  }
  
}; // fin de la classe Euler 





// programme
int main(int argc, char* argv[])
{
  string inputPath("configuration.in.example"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Schema numerique ("Euler"/"E" ou "RungeKutta2"/"RK2" )
  string schema(configFile.get<string>("schema"));

  Engine* engine; // definer la class pour la simulation
  // choisir quel schema utiliser

  if(schema == "B" || schema == "Bonjour!")
  {
	  // initialiser une simulation avec schema Boris-Buneman
    engine = new EngineBorisBuneman(configFile);
  }
  else if(schema == "Euler" || schema == "E")
  {
    // initialiser une simulation avec schema Euler
    engine = new EngineEuler(configFile);
  }
  else if(schema == "RungeKutta2" || schema == "RK2")
  {
    // initialiser une simulation avec schema runge-kutta 2
    engine = new EngineRungeKutta2(configFile);
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
