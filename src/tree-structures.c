#include "sampling.h"

NODEPTR make_node(double pr) {
  NODEPTR newnode;
  newnode = (struct Node *) R_alloc(1, sizeof(struct Node));
  newnode->prob = pr;
  newnode->update = 0;
  newnode->logmarg = 0.0;
  newnode->one = NULL;
  newnode->zero = NULL;
  return(newnode);
}

typedef int (*compfn)( const void* , const void*);

int sortvars(struct Var *vars, double *prob, int p)
{
  int i, n;

  /* Fill in variable information. */


  for (i = 0; i < p; i++) {
    vars[i].prob = prob[i];
    vars[i].index = i;
  }

  /* Make "list" from "probs".  Involves sorting and flipping and such. */
  n = 0;
  for (i = 0; i < p; i++) {
    if (vars[i].prob < 0.0) {
      REprintf("Warning: Probability %d (%lf) less than zero, setting to zero.\n",
	      i, vars[i].prob);
      vars[i].leaveout = TRUE;
      vars[i].prob = 0.0;
    }
    else if (vars[i].prob == 0.0)
      vars[i].leaveout = TRUE;	/* Must be out. */
    else if (vars[i].prob < .5) {
      vars[i].leaveout = FALSE;
      vars[i].logit = log((1.0-vars[i].prob)/(vars[i].prob));
      vars[i].flip = TRUE;
      n++;
    }
    else if (vars[i].prob < 1.0) {
      vars[i].leaveout = FALSE;
      vars[i].logit = log((vars[i].prob)/(1.0-vars[i].prob));
      vars[i].flip = FALSE;
      n++;
    }
    else if (vars[i].prob == 1.0)
      vars[i].leaveout = TRUE;	/* Must be in. */
    else {
      REprintf("Warning: Probability %d (%lf) more than one, setting to one.\n",
	      i, vars[i].prob);
      vars[i].leaveout = TRUE;
      vars[i].prob = 1.0;
    }
  }

  if (n == 0) {
    error("Probabilities are all 0 or 1 - Quitting!\n");
  }
  /* Ok, vars is set up.  Need to sort to get "list". */
  qsort((char *) vars, p, sizeof(struct Var),(compfn) compare);

  return(n);
}


int compare(struct Var *i, struct Var *j)
{
  if (i->leaveout) return(1);
  if (j->leaveout) return(-1);
  if (i->logit > j->logit)
    return (-1);
  if (i->logit < j->logit)
    return (1);
  return (0);
}


int update_probs(double *probs, struct Var *vars, int m, int k, int p) {
  int i, update;
  double wt, newprob, diff;
  wt = (double) m / (double) k;
  update = 0;
  for (i=0, diff=0.0; i <p; i++){
    diff += (probs[vars[i].index] - vars[i].prob)*(probs[vars[i].index] - vars[i].prob);
  }	
  if (sqrt(diff/ (double) p) >  .025) {
    update = 1;
    for (i = 0;   i < p; i++) {
      if (m < p) {
	newprob = probs[vars[i].index]* wt + vars[i].prob * (1 - wt);
      }
      else  newprob = probs[vars[i].index];
      if (newprob > .975) newprob = .975;
      if (newprob < .025) newprob = .025;
      vars[i].prob = newprob;
    }}
return(update);
}

