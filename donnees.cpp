#include "donnees.h"

int inttest=0;


chargement::chargement()
{
  fe=0;
}

sortie::sortie()
{
  fs=0;
  lignenb=0;
}

segment::segment()
{
  date=0;
  temperature=0;
  pas=0;
}

int chargement::ouverture(const char *nomfe)
{
	if(fe=fopen(nomfe,"r")) return 1;
	printf("\nNom de fichier de donnees inconnu : %s\n\n\n",nomfe);
	return 0;
}


void chargement::lignesuivante() 
{
	char c;c=getc(fe);while(c!='\n') c=getc(fe);
}


void chargement::fermeture()
{
	fclose(fe);
}


int sortie::ouverture (const char *nomfs)
{
	if(fs=fopen(nomfs,"w"))
	{
		lignenb=1;
		return 1;
	}
	printf("\nImpossible d'ouvrir le fichier de sortie : %s\n\n\n",nomfs);
	return 0;
}


void sortie::fermeture()
{
	fclose(fs);
}


