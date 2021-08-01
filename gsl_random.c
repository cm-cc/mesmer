/*---------------------------------------- gslwr.c --*/
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*LOOK HERE FOR AVAILABLE GENERATORS: the name has to be put in the arg of gsl_rng_alloc() */
/*http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html#Random-number-generator-algorithms*/
/*http://www.nongnu.org/gsl-shell/doc/random.html*/

static gsl_rng* r;  				/*  r is (almost) never used explicitly and does 
                        			    not need to be seen by the Fortran program.  */

void gslrndinit_(int* s){
  /* some generators: gsl_rng_mt19937 gsl_rng_ranlux gsl_rng_ranlux389 gsl_ranlxs2 */
   r = gsl_rng_alloc(gsl_rng_ranlxd2);		/*  constant gsl_rng_taus is unknown to Fortran  */
   gsl_rng_set(r, (unsigned long int)(*s));	/*  s is cast from int to (unsigned long int)  */
}

void gslrnd_(double* x, int* n){
  int i;
  for(i=0; i<*n; i++) 
    x[i] = gsl_rng_uniform(r);
}

