#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#else
#include <unistd.h>
#endif

#include "iso8859.h"
#include "randtest.h"

#define FALSE 0
#define TRUE  1

#ifdef M_PI
#define PI	 M_PI
#else
#define PI	 3.14159265358979323846
#endif

extern double pochisq(const double ax, const int df);


#ifdef _WIN32	
static int optind = 1;

static int getopt(int argc, char *argv[], char *opts)
{
    static char *opp = NULL;
    int o;
    
    while (opp == NULL) {
        if ((optind >= argc) || (*argv[optind] != '-')) {
	   return -1;
	}
	opp = argv[optind] + 1;
	optind++;
	if (*opp == 0) {
	    opp = NULL;
	}	
    }
    o = *opp++;
    if (*opp == 0) { 
	opp = NULL;
    }
    return strchr(opts, o) == NULL ? '?' : o;
}
#endif

/*  Main program  */

int main(int argc, char *argv[])
{
	int i, oc, opt;
	long ccount[256];	      /* Bins to count occurrences of values */
	long totalc = 0;	      /* Total character count */
	char *samp;
	double montepi, chip,
	       scc, ent, mean, chisq;
	FILE *fp = stdin;
	int counts = FALSE,	      /* Print character counts */
		  fold = FALSE,	      /* Fold upper to lower */
		  binary = FALSE,	  /* Treat input as a bitstream */
		  terse = FALSE;	  /* Terse (CSV format) output */

	if (optind < argc) {
	   if (optind != (argc - 1)) {
              printf("Duplicate file name.\n");
	     // help();
	      return 2;
	   }
           if ((fp = fopen(argv[optind], "rb")) == NULL) {
              printf("Cannot open file %s\n", argv[optind]);
	      return 2;
	   }
	}
	
#ifdef _WIN32

	/** Set the binary mode */

	    _setmode(_fileno(fp), _O_BINARY);
#endif

        samp = binary ? "bit" : "byte";
	memset(ccount, 0, sizeof ccount);

	/* Initialise for calculations */

	rt_init(binary);

	/* Scan input file and count character occurrences */

	while ((oc = fgetc(fp)) != EOF) {
	   unsigned char ocb;

	   if (fold && isISOalpha(oc) && isISOupper(oc)) {
	      oc = toISOlower(oc);
	   }
	   ocb = (unsigned char) oc;
	   totalc += binary ? 8 : 1;
	   if (binary) {
	    int b;
	    unsigned char ob = ocb;

	    for (b = 0; b < 8; b++) {
		ccount[ob & 1]++;
		ob >>= 1;
	    }
	   } else {
	       ccount[ocb]++;	      /* Update counter for this bin */
	   }
	   rt_add(&ocb, 1);
	}
	fclose(fp);

	/* Complete calculation and return sequence metrics */

	rt_end(&ent, &chisq, &mean, &montepi, &scc);

	if (terse) {
           printf("0,File-%ss,Entropy,Chi-square,Mean,Monte-Carlo-Pi,Serial-Correlation\n",
              binary ? "bit" : "byte");
           printf("1,%ld,%f,%f,%f,%f,%f\n",
	      totalc, ent, chisq, mean, montepi, scc);
	}

	/* Calculate probability of observed distribution occurring from the results of the Chi-Square test */

    	chip = pochisq(chisq, (binary ? 1 : 255));

	/* Print bin counts if requested */

	if (counts) {
	   if (terse) {
              printf("2,Value,Occurrences,Fraction\n");
	   } else {
              printf("Value Char Occurrences Fraction\n");
	   }
	   for (i = 0; i < (binary ? 2 : 256); i++) {
	      if (terse) {
                 printf("3,%d,%ld,%f\n", i,
		    ccount[i], ((double) ccount[i] / totalc));
	      } else {
		 if (ccount[i] > 0) {
                    printf("%3d   %c   %10ld   %f\n", i,

		 /* Shows ISO 8859-1  Latin1 characters and delete other codes and space */

                       (!isISOprint(i) || isISOspace(i)) ? ' ' : i,
		       ccount[i], ((double) ccount[i] / totalc));
		 }
	      }
	   }
	   if (!terse) {
              printf("\nTotal:    %10ld   %f\n\n", totalc, 1.0);
	   }
	}

	/* Print calculated results */

	if (!terse) {
           printf("\nEntropy = %f bits per %s.\n", ent, samp);
           printf("File size: %ld %s \n", totalc, samp);
	   	   printf("\nChi square distribution for %ld samples is %1.2f, and randomly\n",
	      totalc, chisq);
	   if (chip < 0.0001) {
              printf("would exceed this value less than 0.01 percent of the times.\n\n");
	   } else if (chip > 0.9999) {
              printf("would exceed this value more than than 99.99 percent of the times.\n\n");
	   } else {
              printf("would exceed this value %1.2f percent of the times.\n\n",
		 chip * 100);
	   }
	   printf("Arithmetic mean value of data %ss is %1.4f (%.1f = random).\n",
	      samp, mean, binary ? 0.5 : 127.5);
           printf("Monte Carlo value for Pi is %1.9f (error %1.2f percent).\n",
	      montepi, 100.0 * (fabs(PI - montepi) / PI));
           printf("Serial correlation coefficient is ");
	   if (scc >= -99999) {
              printf("%1.6f (totally uncorrelated = 0.0).\n", scc);
	   } else {
              printf("undefined (all values equal!).\n");
	   }
	}

	return 0;
}
