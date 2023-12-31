/*
Copyright ESIEE (2009) 

m.couprie@esiee.fr

This software is an image processing library whose purpose is to be
used primarily for research and teaching.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software. You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

/* authors : J. Cousty - L. Najman and M. Couprie */


/* $Id: mcunionfind.h,v 1.3 2006/02/28 07:49:12 michel Exp $ */
/* ============== */
/* types publics  */
/* ============== */

#include "mcutil.h"

typedef struct {
  int32_t Size;
  int32_t *Fth;
  int32_t *Rank;
} Tarjan;

/* ============== */
/* prototypes     */
/* ============== */

/*
extern Tarjan * CreeTarjan(int32_t taille);
extern void TarjanTermine(Tarjan * T);
extern void TarjanInit(Tarjan * T);
extern void TarjanPrint(Tarjan * T);
extern void TarjanMakeSet(Tarjan * T, int32_t x);
extern int32_t TarjanFind(Tarjan * T, int32_t x);
extern int32_t TarjanLink(Tarjan * T, int32_t x, int32_t y);
extern int32_t TarjanLinkSafe(Tarjan * T, int32_t x, int32_t y);
*/

/* ==================================== */
Tarjan * CreeTarjan(int32_t taille)
/* ==================================== */
/*! \fn Tarjan * CreeTarjan(int32_t taille)
    \param taille : nombre total d'��ents
    \return pointeur sur une structure Tarjan
    \brief cr� une structure pour la fusion d'ensemble
    \warning ne fait pas l'initialisation de la structure
*/
#undef F_NAME
#define F_NAME "CreeTarjan"
{
  Tarjan * T = (Tarjan *)calloc(1,sizeof(Tarjan));
  T->Size = taille;
  if (T == NULL)
  {   fprintf(stderr, "%s : malloc failed for T\n", F_NAME);
      return(NULL);
  }
  T->Fth = (int32_t *)calloc(1,taille * sizeof(int32_t));
  if (T->Fth == NULL)
  {   fprintf(stderr, "%s : malloc failed for T->Fth\n", F_NAME);
      return(NULL);
  }
  T->Rank = (int32_t *)calloc(1,taille * sizeof(int32_t));
  if (T->Rank == NULL)
  {   fprintf(stderr, "%s : malloc failed for T->Rank\n", F_NAME);
      return(NULL);
  }
  return T;
} //CreeTarjan()

/* ==================================== */
void TarjanTermine(Tarjan * T)
/* ==================================== */
/*! \fn void TarjanTermine(Tarjan * T)
    \param T: une structure Tarjan
    \brief lib�e la m�oire
*/
{
  free(T->Fth);
  free(T->Rank);
  free(T);
} //TarjanTermine()


/* ==================================== */
void TarjanMakeSet(Tarjan * T, int32_t x)
/* ==================================== */
/*! \fn void TarjanMakeSet(Tarjan * T, int32_t x)
    \param T: une structure Tarjan
    \param x: un ��ent
    \brief ajoute le singleton {x} �la famille d'ensembles
*/
{
  T->Fth[x] = x;
  T->Rank[x] = 0;
} //TarjanMakeSet()

/* ==================================== */
void TarjanInit(Tarjan * T)
/* ==================================== */
/*! \fn void TarjanInit(Tarjan * T)
    \param T: une structure Tarjan
    \brief initialise la structure (cr� les singletons)
*/
{
  int32_t i;
  for (i = 0; i < T->Size; i++) TarjanMakeSet(T, i);
} //TarjanInit()

/* ==================================== */
void TarjanPrint(Tarjan * T)
/* ==================================== */
/*! \fn void TarjanPrint(Tarjan * T)
    \param T: une structure Tarjan
    \brief affiche la structure
*/
{
  int32_t i;
  for (i = 0; i < T->Size; i++)
    printf("%d: Rank = %d ; Fth = %d\n", i, T->Rank[i], T->Fth[i]);
} //TarjanPrint()



/* ==================================== */
int32_t TarjanFind(Tarjan * T, int32_t x)
/* ==================================== */
/*! \fn int32_t TarjanFind(Tarjan * T, int32_t x)
    \param T: une structure Tarjan
    \param x: un ��ent
    \return un repr�entant
    \brief retourne le repr�entant de l'ensemble auquel appartient x
    \warning x doit appartenir �un ensemble de la famille - pas de v�ification
*/
{
  if (T->Fth[x] != x) T->Fth[x] = TarjanFind(T, T->Fth[x]);
  return T->Fth[x];
} //TarjanFind()

/* ==================================== */
int32_t TarjanLink(Tarjan * T, int32_t x, int32_t y)
/* ==================================== */
/*! \fn int32_t TarjanLink(Tarjan * T, int32_t x, int32_t y)
    \param T: une structure Tarjan
    \param x, y: deux repr�entants
    \return un repr�entant
    \brief fusionne les ensembles repr�ent� par x et y et retourne le repr�entant de la fusion
    \warning x et y doivent �re des repr�entants - pas de v�ification
*/
{
  if (T->Rank[x] > T->Rank[y]) { int32_t tmp = x; x = y; y = tmp; }
  if (T->Rank[x] == T->Rank[y]) T->Rank[y] += 1;
  T->Fth[x] = y;
  return y;
} //TarjanLink()

/* ==================================== */
int32_t TarjanLinkSafe(Tarjan * T, int32_t x, int32_t y)
/* ==================================== */
/*! \fn int32_t TarjanLinkSafe(Tarjan * T, int32_t x, int32_t y)
    \param T: une structure Tarjan
    \param x, y: deux ��ents
    \return un repr�entant
    \brief fusionne les ensembles auxquels appartiennent x et y et retourne le repr�entant de la fusion
*/
{
  x = TarjanFind(T, x);
  y = TarjanFind(T, y);
  if (T->Rank[x] > T->Rank[y]) { int32_t tmp = x; x = y; y = tmp; }
  if (T->Rank[x] == T->Rank[y]) T->Rank[y] += 1;
  T->Fth[x] = y;
  return y;
} //TarjanLinkSafe()
