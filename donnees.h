/*
 *  donnees.h
 *  
 *
 *  Created by nat on 11/03/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef donnees_h
#define donnees_h

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>

#define R 8.3145
#define PI 3.1416
#define kB 1.3807e-23
#define kdelta 1.e-10
#define ZEPS 1.e-10

//---------------------

class chargement
{
	FILE *fe;
public:
	chargement();
	int ouverture (const char *);
	void lecture();
	void lecture1();
	void fermeture();
	void lignesuivante();
};


class sortie
{
	FILE *fs;
	int lignenb;
public:
	sortie();
	int ouverture (const char *);
	void ecriture(std::vector <double> &Y);
	void ecritureX(std::vector <double> &Y);
	void ecritureU(std::vector <double> &U);
	void ecritureV(std::vector <double> &U, int compt);
	void ecriture_growth_dissolution_FeC(std::vector <double> &t, \
										 std::vector <double> &Rlong, \
										 std::vector <double> &Xca, \
										 std::vector <double> &Xfea, \
										 std::vector <double> &Xcb, \
										 std::vector <double> &Xfeb, \
										 std::vector <double> &XV, \
										 std::vector <double> &Gtrf, \
										 std::vector <double> &Gfricf);
	void fermeture();
};



class segment
{
public:
	segment();
	double date;
	double temperature;
	double pas;
};

#endif

