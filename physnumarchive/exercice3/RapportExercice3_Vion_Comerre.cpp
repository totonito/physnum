#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 

using namespace std;
const double pi=3.1415926535897932384626433832795028841971e0;
class Exercice3
{

private:
  double t, dt, tFin;
  double m, g, L;
  double d, Omega, kappa;
  double theta, thetadot;
  int N_excit, nstep_per;
  int sampling;
  int last;
  ofstream *outputFile;

  void printOut(bool force)
  {
    if((!force && last>=sampling) || (force && last!=1))
    {
      double emec = 0.5 * m * pow(L*thetadot , 2) + m*g*L * (1 - cos(theta)) ; // TODO: Evaluer l'energie mecanique
      double pnc = L*L*m*thetadot*acceleration(theta,thetadot,t) + L * thetadot * m * g * sin(theta); // TODO: Evaluer la puissance des forces non conservatives

      *outputFile << t << " " << theta << " " << thetadot << " " << emec << " " << pnc << endl;
     
      last = 1;
    }
    else
    {
      last++;
    }
  }

  double acceleration(const double theta_, const double thetadot_, const double t_)
  {
    // TODO: Modifier selon l'expression analytique
    return (-kappa / m) * thetadot_ - (g/L) * sin(theta_) + (d/L) * pow(Omega,2) * sin(Omega*t_) * sin (theta_);
  }

  void step()
  {
    // TODO: Modifier  selon l'algorithme  
    double theta_copie (theta);
    double thetadot_demi (0);
    theta += thetadot * dt + 0.5 * acceleration(theta, thetadot,t) * dt * dt;
    thetadot_demi = thetadot + 0.5 *  acceleration(theta_copie, thetadot,t) * dt;
    thetadot += 0.5 * dt * (acceleration(theta_copie , thetadot_demi,t) + acceleration(theta,thetadot_demi,t+dt));
    
    // TEST avec euler explicite de la fonction accéleration
    /*theta += thetadot * dt ;
    thetadot += acceleration (theta , thetadot ,t)*dt;*/
  }


public:
  Exercice3(int argc, char* argv[])
  {
    const double pi=3.1415926535897932384626433832795028841971e0;
    string inputPath("configuration.in.example"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

        // t final (overwritten if N_excit >0)
       // time step (overwritten if nstep_per >0)
    d        = configFile.get<double>("d");         // amplitude forcing term
    Omega    = configFile.get<double>("Omega");     // angular frequency forcing term 
    kappa    = configFile.get<double>("kappa",kappa);     // coefficient for friction
    m        = configFile.get<double>("m");         // mass
    g        = configFile.get<double>("g");         // gravity acceleration
    L        = configFile.get<double>("L");         // length
    theta    = configFile.get<double>("theta0");    // initial condition in theta
    thetadot = configFile.get<double>("thetadot0"); // initial condition in thetadot
    sampling = configFile.get<int>("sampling");     // number of time steps between two writings on file
    N_excit  = configFile.get<int>("N");            // number of periods of excitation
    nstep_per= configFile.get<int>("nstep");        // number of time step per period
	dt       = configFile.get<double>("dt"); 
	tFin     = configFile.get<double>("tFin"); 
    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
    if(N_excit>0){
      tFin = N_excit*(2*pi/Omega);
    }
     if(nstep_per>0){
      dt = 2.0*pi/(nstep_per * Omega);
    }
  }

  ~Exercice3()
  {
    outputFile->close();
    delete outputFile;
  };

  void run()
  {
    t = 0.;
    last = 0;
    printOut(true);
    while( t < tFin-0.5*dt )
    {
      step();
      t += dt;
      printOut(false);
    }
    printOut(true);
  };

};

int main(int argc, char* argv[])
{
  Exercice3 engine(argc, argv);
  engine.run();
  return 0;
}
