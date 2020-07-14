/* Libs */
	#include<igraph.h>
	#include<stdio.h>
	#include<math.h>   
	#include<errno.h>  
	#include<stdlib.h> 
	#include<string.h>
  #include <time.h>

/* The main business */

  int main(int argc, char *argv[]) {
	
  /* Reads the inputs and variable declarations */
      FILE *ifile = fopen(argv[1], "r");
      FILE *outfile = fopen(argv[2], "w");
      double delta = atof(argv[3]);
      FILE *pkfile = fopen(argv[4], "r");

    	igraph_t g;
    	igraph_vector_t escore, nrgeo, aux, aux2, path, vecprint, vecprint2, avbias, vecinos;
      igraph_vector_ptr_t links;
      igraph_vs_t summit;
      igraph_integer_t to, from;
    	int i, j, id, mut, Ns, t, k, v, w, wmaxe, edges, pos, transversions, transitions, len, inicio, l, accessible, stop, cunha;
      int actrans, ac, size, sumlen, x, x2, count, c, a, q, peak[200];
      double dice, maxe, tiav, tvav, dba, bias, acbias, record, final, effective, effective2;
      char string[200][9];
      char line[9];
    	const char E[10] = "score";
    	const char W[10] = "genotype";

    	igraph_i_set_attribute_table(&igraph_cattribute_table);
    	igraph_read_graph_gml(&g, ifile);
    	fclose(ifile);

    	Ns = igraph_vcount(&g);
    	char gene[Ns + 1][9];
    	char reverse[Ns + 1][9];
      char geneaux[Ns + 1][9];
      edges = igraph_ecount(&g);
      igraph_vector_init(&escore, 1);
      igraph_vector_init(&aux, 1);
      igraph_vector_init(&aux2, Ns);

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
          strcpy(reverse[i], gene[i]);
          a = 1;
          for(j = 0; j < q; j++) {
            a = strcmp(string[j], gene[i]);
            if (a == 0)
              peak[j] = i;
        }
    	}

      VANV(&g,E,&escore);
      wmaxe = igraph_vector_which_max(&escore);
      igraph_vs_1(&summit, wmaxe);

  /* Loop to get the reverse comps */
      for(i = 0; i < Ns; i++) {
          for(j = 0; j < 8; j++){
              if (gene[i][j] == 'A')
                reverse[i][7-j] = 'T';
              else if (gene[i][j] == 'T')
                reverse[i][7-j] = 'A';
              else if (gene[i][j] == 'C')
                reverse[i][7-j] = 'G';
              else
                reverse[i][7-j] = 'C';
          }
      }
        tiav = 0;
        tvav = 0;
        t = 0;
        a = 0;
        sumlen = 0;
        x = 0;
        x2 = 0;

    igraph_vector_update(&aux, &escore);
    igraph_vector_sort(&aux);

/* Loop over all sequences to get the ones with the lowest E-score */
    for (i = 0; i < Ns; i++) {
      if (VECTOR(escore)[i] < VECTOR(aux)[(int)(Ns / 10)]) {
        VECTOR(aux2)[t] = i;
        t += 1;
      }
    }
    igraph_vector_destroy(&aux);

    igraph_vector_init(&vecprint, 1);
    igraph_vector_init(&vecprint2, 1);
    igraph_vector_init(&avbias, 1);
    igraph_vector_resize(&avbias,t);
    igraph_vector_init(&vecinos, t);

/* Loop over all initial conditions */
  for (k = 0; k < t; k++) {

    count = 0;
    record = 0;
    transitions = 0;
    igraph_vector_ptr_init(&links, 1);
    igraph_vector_init(&nrgeo, 1);

    inicio = (int)VECTOR(aux2)[k];

    igraph_get_all_shortest_paths(&g, &links, &nrgeo, inicio, summit, IGRAPH_ALL);

    len = (int)igraph_vector_ptr_size(&links);

    sumlen += len;
    igraph_vector_resize(&vecprint,sumlen);
    igraph_vector_resize(&vecprint2,sumlen);

  /* Loop over all paths */
    for (w = 0; w < len; w++) {
      transversions = 0;
      transitions = 0;
      actrans = 0;
      ac = 0;

    /* Copies the shortest path w */  
      igraph_vector_init(&path, 1);
      igraph_vector_copy(&path, VECTOR(links)[w]);

      accessible = 0;
      cunha = 0;
      stop = igraph_vector_size(&path) - 1;
      for (l = 0; l < igraph_vector_size(&path) - 1; l++) {
        from = VECTOR(path)[l];
        to = VECTOR(path)[l + 1];
        if(cunha == 10){
        for(j = 0; j < q; j++) {
            if (to == peak[j]){
              stop = l + 1;
              cunha += 1;
              break;
            }
        }
        }
        if(VECTOR(escore)[to] > VECTOR(escore)[from] - delta) {
          accessible = 1;
          break;
        } 
      }

    /* Loop over all edges */
      for (l = 0; l < stop; l++) {

        from = VECTOR(path)[l];
        to = VECTOR(path)[l + 1];

        dba = fabs(VECTOR(escore)[from] - VECTOR(escore)[to]); 
    	  j = 0;
        /* Finds and records the diferences between the sequences */
    	  for (i = 0; i < 8; i++){
    	   	if (gene[from][i] != gene[to][i]){
          		j += 1;
          		pos = i;
          }
        }
      	/* If there's more than one single difference, check the reverse complement */
      	if (j == 1)
          strcpy(geneaux[to], gene[to]);
        else{
          strcpy(geneaux[to], reverse[to]);
      		j = 0;
      		/* Finds and records the differences between the sequence and the reverse complement */
        	for (i = 0; i < 8; i++){
    			  if (gene[from][i] != geneaux[to][i]){
            			j += 1;
            			pos = i;
        		}
        	}
        }
        if(j == 1){
        	if((gene[from][pos] == 'A')||(gene[from][pos] == 'G')){
        		if ((geneaux[to][pos] == 'A')||(geneaux[to][pos] == 'G')){
        			transitions += 1;
              if(accessible)
                actrans += 1;              

        		}
        		else{
        			transversions += 1;
              if(accessible)
                ac += 1; 
        		}
        	}
        	else{
        		if ((geneaux[to][pos] == 'T')||(geneaux[to][pos] == 'C')){
        			transitions += 1;

              if(accessible)
                actrans += 1;
        		}	
        		else{
        			transversions += 1;
              if(accessible)
                ac += 1; 
        		}
        	}
        }
      }
      effective = actrans + (ac / 2.0);
      effective2 = transitions + (transversions / 2.0);
      VECTOR(vecprint)[x] = (double)transitions/(double)effective2;
      if(accessible){
        VECTOR(vecprint2)[x2] = (double)actrans/(double)effective;
        x2 += 1;
        count += 1;
        record += (double)actrans/(double)effective;
      }
      x += 1;

      igraph_vector_destroy(&path);

     }

    VECTOR(avbias)[k] = (double)record / (double)len;
    igraph_vector_ptr_destroy(&links);
    igraph_vector_destroy(&nrgeo);
                                                                  
  }

  /* Print the thing */

    for (k = 0; k < x; k++) {
      fprintf(outfile,"%f\t", VECTOR(vecprint)[k]);
      if (k < x2)
        fprintf(outfile,"%f\n", VECTOR(vecprint2)[k]);
      else
        fprintf(outfile,"\n");
    }

    c = 0;
    fprintf(outfile,"--------\n");
    for (k = 0; k < igraph_vector_size(&avbias); k++) {
      if(VECTOR(avbias)[k] == VECTOR(avbias)[k]){
        final += VECTOR(avbias)[k];
        c += 1;
      }
    }
    final /= c;
    fprintf(outfile,"%f\n", final);

    igraph_vector_destroy(&vecprint);
    igraph_vector_destroy(&vecprint2);
    igraph_destroy(&g);
    
    return 0;

  }
