cmake_minimum_required(VERSION 3.23)
project(physnum)

set(CMAKE_CXX_STANDARD 14)

include_directories("Exercice 4")

add_executable(physnum
        "Exercice 4/Analyse_Apollo_student.m"
        "Exercice 4/Analyse_rho_descendre_Comerre_Vion.m"
        "Exercice 4/Analyse_rho_monter_Comerre_Vion.m"
        "Exercice 4/ConfigFile.h"
        "Exercice 4/ConfigFile.hpp"
        "Exercice 4/Exercice4_Comerre_Vion.cpp"
        "Exercice 4/Exercice4_Comerre_Vion_rho.cpp"
        physnumarchive/exercice3/RapportExercice3_Vion_Comerre/Questiona_Vion_Comerre.m
        physnumarchive/exercice3/RapportExercice3_Vion_Comerre/Questionb_Vion_Comerre.m
        physnumarchive/exercice3/RapportExercice3_Vion_Comerre/Questionc_Vion_Comerre.m
        physnumarchive/exercice3/RapportExercice3_Vion_Comerre/Questiond2_Vion_Comerre.m
        physnumarchive/exercice3/ConfigFile.h
        physnumarchive/exercice3/ConfigFile.hpp
        physnumarchive/exercice3/ParameterScan_Vion_Comerre.m
        physnumarchive/exercice3/RapportExercice3_Vion_Comerre.cpp
        physnumarchive/exercice3/somme.cpp
        physnumarchive/zExercice1_student/Analyse_Exercice1_Cas260secondes_Vion_Comerre.m
        physnumarchive/zExercice1_student/Analyse_Exercice1_Convergence_frottements_Vion_Comerre.m
        physnumarchive/zExercice1_student/Analyse_Exercice1_Convergence_sans_frottements_Vion_Comerre.m
        physnumarchive/zExercice1_student/ConfigFile.h
        physnumarchive/zExercice1_student/Exercice1_Comerre_Vion.cpp
        physnumarchive/zexercice2/Analyse.m
        physnumarchive/zexercice2/Analyse_Exercice1_Vion_Comerre.m
        physnumarchive/zexercice2/ConfigFile.hpp
        physnumarchive/zexercice2/Exercice2_Comerre_Vion.cpp
        physnumarchive/zexercice2/ParameterScan_2.m)
