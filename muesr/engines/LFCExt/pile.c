#include "pile.h"



void pile_init(pile * p, unsigned int nElements)
{
	p->nElements = nElements;
	p->ranks = malloc(nElements * sizeof(double));
	p->elements = malloc(nElements * sizeof(struct vec3));
    
    unsigned int i;
    for (i = 0; i < nElements; ++i)
    {
        p->ranks[i] = -1.0;
        p->elements[i] = vec3_zero();
    }
}

void pile_add_element(pile * p, double rank, struct vec3 v)
{
	unsigned int i;
	for ( i = 0; i < p->nElements; i++)
	{
		if (p->ranks[i] == -1.0) {
			p->ranks[i] = rank;
			p->elements[i] = v;
			break;
		}
		if (p->ranks[i] > rank) {
			pile_move_elements_from_position(p, i);
			p->ranks[i] = rank;
			p->elements[i] = v;
			break;
		}
	}
}

void pile_move_elements_from_position(pile * p, unsigned int pos)
{
	// the first -1 is for 0 indexing
	// the second -1 is becouse if only the last element must be moved it is thrashed!
	if (p->nElements < 2) {
		return;
	}
	unsigned int i;
	for (i = (p->nElements-1) ; i-- > pos ;)
	{
		p->ranks[i+1] = p->ranks[i];
		p->elements[i+1] = p->elements[i];
	}
}

void pile_free(pile * p)
{
	free(p->ranks);
	free(p->elements);
}
