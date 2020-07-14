/* Libs */
	#include<igraph.h>
	#include<stdio.h>
	#include<math.h>   
	#include<errno.h>  
	#include<stdlib.h> 
	#include<string.h>
  #include <time.h>

/* Functions */
  void initialstate(int n, int Ns, struct igraph_vector_t pop, struct igraph_vector_t escore);
	void pleaseprint(int n, int Ns, int t, int w, int generations, int wmaxe, igraph_t g,
                      struct igraph_vector_t pop, struct igraph_vector_t escore, 
                      struct igraph_vector_t gc, struct igraph_matrix_t AUX, char input[][9], int peak[], int q);
  void steamroller(int generations, int repeat, struct igraph_matrix_t AUX, struct igraph_matrix_t matriz);
	void create_roulette(int n, struct igraph_vector_t pop, struct igraph_vector_t escore, 
                        struct igraph_vector_t roulette_wheel);
	char* mutation(char[], double bias);

/* The main business */

  int main(int argc, char *argv[]) {
	
  /* Reads the inputs and variable declarations */
    	FILE *ifile = fopen(argv[1], "r");
      FILE *outfile = fopen(argv[2], "w");
      int generations = atoi(argv[3]);
      int repeat = atoi(argv[4]);
      int n = atoi(argv[5]);
      double prob = atof(argv[6]);
      double bias = atof(argv[7]);
      int seed = atoi(argv[8]);
      FILE *pkfile = fopen(argv[9], "r");
    	igraph_t g;
    	igraph_vector_t nahe, escore, ranges, gc, pop, cc, cg, auxpop;
    	igraph_matrix_t output, AUX;
    	int i, j, id, mut, Ns, t, k, v, w, aux, wmaxe, q, a, peak[200];
    	double dice, maxe, maxgc;
      const char E[10] = "score";
      const char W[10] = "genotype";
    	char son[1][9];
      char string[200][9];
      char line[9];

      srand(time(NULL));                          // randomize seed
    	igraph_i_set_attribute_table(&igraph_cattribute_table);
    	igraph_read_graph_gml(&g, ifile); 
    	fclose(ifile);
      srand(seed);

    	Ns = igraph_vcount(&g);
    	char gene[Ns + 1][9];
      n = (int)(n);                                     // population size
      igraph_vector_init(&nahe, 1);
      igraph_vector_init(&pop, n + 1);
      igraph_vector_init(&auxpop, n + 1);
      igraph_vector_init(&cc, 1);
      igraph_vector_init(&cg, 1);
      igraph_vector_init(&escore, 1);
      igraph_vector_init(&ranges, n + 1);
      igraph_vector_init(&gc, 1); 
      igraph_matrix_init(&output, 1, 1);
      igraph_matrix_init(&AUX, 1, 1);
    	igraph_vector_resize(&gc, Ns);
    	igraph_vector_resize(&cc, Ns + 1);
    	igraph_vector_resize(&cg, Ns + 1);
    	igraph_matrix_resize(&output, generations + 1, 10);
    	igraph_matrix_null(&output);
      igraph_matrix_resize(&AUX, generations * repeat, 10);
      igraph_matrix_null(&AUX);

  /* Reads the sequences in the global peak */
      q=0;
      while(!feof(pkfile)) {
        fscanf(pkfile,"%s",line);
        strcpy(string[q], line);
        q+=1;
      }
      fclose(pkfile);

  /* Reads the Escores and sequences of all the vertex in the input network and saves the ones in the peak */
    	for(i = 0; i < Ns; i++) {
      		strcpy(gene[i], VAS(&g,W,i));
          a = 1;
          for(j = 0; j < q; j++) {
            a = strcmp(string[j], gene[i]);
            if (a == 0)
              peak[j] = i;
        }
    	}
     	VANV(&g,E,&escore);
     	maxe = igraph_vector_max(&escore);
     	wmaxe = igraph_vector_which_max(&escore);

  /* Calculates the GC content of all the vertex in the network */
    	for(j = 0; j < Ns; j++){ 
    		for (i = 0, VECTOR(cc)[j] = 0, VECTOR(cg)[j]=0; gene[j][i]; i++){
        	VECTOR(cc)[j] += (gene[j][i] == 'C');
        	VECTOR(cg)[j] += (gene[j][i] == 'G');
        }
       	VECTOR(gc)[j]=(VECTOR(cc)[j]+VECTOR(cg)[j])/8.0;
    	}
      maxgc = igraph_vector_sum(&gc) / Ns;

  /* Loop over all runs */
    	for (w = 0; w < repeat; w++) {
      	t = 0;
      	mut = 0;

  /* Sets and prints the initial conditions */
    	initialstate(n, Ns, pop, escore);

      pleaseprint(n, Ns, t, w, generations, wmaxe, g, pop, escore, gc, AUX, gene, peak, q);

  /* Loop over all generations */
  		for (t = 1; t < generations; t++) {           

	    	/* Function that creates the roulette wheel */
	    		create_roulette(n, pop, escore, ranges); 

        	/* Keeps the record of the parents id for the roulette wheel */
          		igraph_vector_update(&auxpop,&pop);

      		/* Loop over all individuals in the population to select the parents */
				  for(j = 0; j < n; j++) {

        		/* Turn the roulette wheel */
              		dice = rand() / (RAND_MAX + 1.);
        			for(i = 0; i < n; i++) {
        	  			if(dice < VECTOR(ranges)[i]) {
        		  			id = (int)VECTOR(auxpop)[i];
        	  				break;
          				}
          			}
        		/* Mutation process */
       				if(rand() / (RAND_MAX + 1.) < prob) {                //Probability of having a mutation
	            		igraph_neighbors(&g, &nahe, id, IGRAPH_ALL);     //Gets the neighbours of the parent vertex
                  		aux = 0;
            			while(aux == 0){
                    		strcpy(son[0], mutation(gene[id],bias));       //Copies the sequence of the parent with mutation
                    		aux = strcmp(son[0], gene[id]);
                  		}      

        			/* Loop over all neighbours of the parent vertex to check if the mutation belongs to the network */
            			for(k = 0; k < igraph_vector_size(&nahe); k++){
              				v = (int)VECTOR(nahe)[k];
                      		aux = strcmp(son[0], gene[v]);
	              			/* Verification of the mutation viability */
                      		if(aux == 0){
                        		mut += 1;
                        		VECTOR(pop)[j] = v;
                        		break;
                      		}
            			}
            			/* If the mutation is not part of the network, choose another parent */
            			if(mut < 1){
                  			j = j - 1;
                  		}
        			}
        			/* Case without mutation */    
        			else{
          				VECTOR(pop)[j] = id;
          			}
        			mut = 0;
    			}
  		}		

      for(i = 0; i < n; i++) {
      fprintf(outfile,"%f\t", VECTOR(pop)[i]);
    }
    fprintf(outfile,"\n");

    }

  /* Frees memory and that's it */
    	igraph_destroy(&g);
   		igraph_vector_destroy(&nahe);
    	igraph_vector_destroy(&escore);
    	igraph_vector_destroy(&ranges);
    	igraph_vector_destroy(&gc);
    	igraph_vector_destroy(&pop);
    	igraph_vector_destroy(&cc);
    	igraph_vector_destroy(&cg);
    	igraph_matrix_destroy(&output);
    
    	fclose(outfile); 
    
    	return 0;
  
  }

/* This one sets the initial conditions */
  void initialstate(int n, int Ns, struct igraph_vector_t pop, struct igraph_vector_t escore) {

    int i, j = 0, v;
    igraph_vector_t aux, aux2;
    igraph_vector_init(&aux, 1);
    igraph_vector_init(&aux2, Ns);

    igraph_vector_update(&aux, &escore);
    igraph_vector_sort(&aux);

    for (i = 0; i < Ns; i++) {
      if (VECTOR(escore)[i] == VECTOR(aux)[0]) {
        VECTOR(aux2)[j] = i;
        j += 1;
      }
    }

    //v = rand() % j;
    for (i = 0; i < n; i++) {
      VECTOR(pop)[i] = VECTOR(aux2)[0];
    }

    igraph_vector_destroy(&aux);
    igraph_vector_destroy(&aux2);
  }

/* This one is the print function for each generation */

	void pleaseprint(int n, int Ns, int t, int w, int generations, int wmaxe, igraph_t g,
                      struct igraph_vector_t pop, struct igraph_vector_t escore, 
                      struct igraph_vector_t gc, struct igraph_matrix_t AUX, char input[][9], int peak[], int q) {

  	int i, j, v;
    double standEsc = 0, standGCc = 0, enormfactor = 0.0, gcnormfactor = 0.0, fc, fg, ft, fa, div = 0.0, max;
  	igraph_vector_t aux, aux2, aux3, aux4;
  	igraph_vector_init(&aux, n + 1);
  	igraph_vector_fill(&aux, 0.0);
    igraph_vector_init(&aux2, n + 1);
    igraph_vector_fill(&aux2, 0.0);
    igraph_vector_init(&aux3, n + 1);
    igraph_vector_fill(&aux3, 0.0);
    igraph_vector_init(&aux4, n);
    igraph_vector_fill(&aux4, 0.0);

  	/* Loop over all individuals in the population to calculate the averages Escore and GC content */
  		for(i = 0; i < n; i++) {
    		v=(int)VECTOR(pop)[i];
    		enormfactor += VECTOR(escore)[v];
    		gcnormfactor += VECTOR(gc)[v];
    		VECTOR(aux)[i] = VECTOR(escore)[v];
        	VECTOR(aux2)[i] = VECTOR(gc)[v];
  		}
  
      MATRIX(AUX, t + generations*w, 0) = enormfactor/n;
      MATRIX(AUX, t + generations*w, 1) = gcnormfactor/n;

      /* Normalizes the E-score */
      max = igraph_vector_max(&escore);
      MATRIX(AUX, t + generations*w, 0) = MATRIX(AUX, t + generations*w, 0)/max;


    /* Loop over all individuals in the population to find the number of different individuals */
  		for(i = 0; i < n; i++) {
       		for (j = 0; j < i; j++) {
           		if ((int)VECTOR(pop)[i] == (int)VECTOR(pop)[j])
             	break;
        	}
            if (i == j){
             	MATRIX(AUX, t + generations * w, 3) += 1;
            }
  		}

    /* Loop over all individuals in the population to find global peak is reached */ 
      for(i = 0; i < n; i++) {
        for(j = 0; j < q; j++) {
            if (VECTOR(pop)[i] == peak[j]){
              MATRIX(AUX, t + generations * w, 4) += 1;
              break;
            }
        }
      }
      MATRIX(AUX, t + generations * w, 4) /= n;

    igraph_vector_destroy(&aux);
    igraph_vector_destroy(&aux2);
    igraph_vector_destroy(&aux3);
    igraph_vector_destroy(&aux4);
	}

/* This one prints the output */

  void steamroller(int generations, int repeat, struct igraph_matrix_t AUX, struct igraph_matrix_t matriz) {

    int i, j, k;
    double stand = 0.0;

    /* Calculates averages */

      for (j = 0; j < 5; j++) {
        for (i = 0; i < generations; i++) {
          for (k = 0; k < repeat; k++) {
            MATRIX(matriz, i, j) += MATRIX(AUX, i + generations * k, j);
          }
          MATRIX(matriz, i, j) /= (double)repeat;
        }  
      }

    /* Calculates standard deviations */

      for (j = 0; j < 5; j++){
        for(i = 0; i < generations; i++) {
          for (k = 0; k < repeat; k++) {
            stand += pow(MATRIX(AUX, i + generations * k, j) - MATRIX(matriz, i, j), 2);
          }
        MATRIX(matriz, i, j + 5) = sqrt(stand / (double)(repeat));
        stand = 0.0;
        }        
      }

  }

/* This one creates the roulette */
	void create_roulette(int n, struct igraph_vector_t pop, struct igraph_vector_t escore, 
                      struct igraph_vector_t roulette_wheel) {
	
		double normfactor = 0.0;
		igraph_vector_t pop_escore;
    igraph_vector_init(&pop_escore, n + 1);
		int i;

	/* Loop over all individuals in the population to calculate the normalization factor */
		for(i = 0; i < n; i++) {
			VECTOR(pop_escore)[i] = VECTOR(escore)[(int)VECTOR(pop)[i]] + (rand() / (RAND_MAX + 1.))*0.001 - 0.0005;
		}

		normfactor = igraph_vector_sum(&pop_escore);

		igraph_vector_scale(&pop_escore, 1/normfactor);

	/* Loop over all individuals in the population to create the roulette wheel */

		VECTOR(roulette_wheel)[0] = VECTOR(pop_escore)[0];

		for(i = 1; i <= n; ++i) {
			VECTOR(roulette_wheel)[i] = VECTOR(roulette_wheel)[i - 1] + VECTOR(pop_escore)[i];
		}

    igraph_vector_destroy(&pop_escore);
	}

/* This one makes mutations */
	 char *mutation(char input[], double bias) {
    	int i, n;
    	char * str = malloc(9);
    	double random_value, random_two;
    	random_value = rand() / (RAND_MAX + 1.);
      random_two = rand() / (RAND_MAX + 1.);

    	strcpy(str,input);
    	n = rand() % 7 + 1;
    	if(str[n] == 'A'){
    	  if(random_value < bias)
    	    str[n] = 'G';
    	  else{
          if(random_two < 0.5)
            str[n] = 'T';
          else 
            str[n] = 'C';
         } 
        goto END;
    	}
    	else if(str[n] == 'C'){
    	  if(random_value < bias)
          str[n] = 'T';
        else{
          if(random_two < 0.5)
            str[n] = 'G';
          else 
            str[n] = 'A';
         } 
    	  goto END;
    	}
    	else if(str[n] == 'T'){
    	  if(random_value < bias)
          str[n] = 'C';
        else{
          if(random_two < 0.5)
            str[n] = 'A';
          else 
            str[n] = 'G';
         } 
    	  goto END;
    	}
    	else {
    	  if(random_value < bias)
          str[n] = 'A';
        else{
          if(random_two < 0.5)
            str[n] = 'T';
          else 
            str[n] = 'C';
         } 
    	  goto END;
   		}

    	END: return str;
	 }
