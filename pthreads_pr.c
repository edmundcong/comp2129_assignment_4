#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <immintrin.h>

#include "pagerank.h"

typedef struct vect_calc {
  double vectNorm;
  int *p_inlinks;
  int *p_outlinks;
  double *currRank;
  double leading_term;
  double *prevRank;
  double dampener;
  int **inlinks_data;
  int start;
  int end;
  int j;
} vect_struct;

void *calc_worker(void *calc_struct){
  vect_struct calc = *(vect_struct *)calc_struct; //casting
  for (int i = calc.start; i < calc.end; i++) {
  printf("%d\n", calc.start);
    calc.vectNorm = 0;
      if (calc.p_inlinks[i] == 0) {
          for (int j = calc.j; j < calc.p_inlinks[i]; j++) {
            calc.vectNorm = calc.vectNorm + calc.prevRank[calc.inlinks_data[i][j]]
            /(calc.p_outlinks[calc.inlinks_data[i][j]]);
          }
      }
      //assigning currRank[i] to newest t row result
      calc.currRank[i] = calc.leading_term + calc.dampener * calc.vectNorm;
  }

  return NULL;
}

void pagerank(node *list, int npages, int nedges, int nthreads, double dampener) {
    double vectNorm        = EPSILON;
    double squaredEpsilon  = EPSILON*EPSILON;
    double leading_term    = (1.0 - dampener)/npages;
    double difference_vect = 0;
    double zeroth_vect     = 1.0 / npages;
    double *tempRank;

    double *currRank = malloc(npages*sizeof(double));
    double *prevRank = malloc(npages*sizeof(double));

    int * *inlinks_data = (int * *)malloc(npages*sizeof(int *));
    int *p_inlinks      = malloc(npages*sizeof(int));
    int *p_outlinks     = malloc(npages*sizeof(int));
    int inlinks_count = 0;

    pthread_t t_id[nthreads];
    vect_struct *calc = malloc(sizeof(vect_struct)*nthreads);
    int chunk = 0;

    struct node *list_copy = malloc(sizeof(struct node));
    memcpy(list_copy, list, sizeof(node));

    // populate arrays for storing page data
    node *temp_node = list;
    node *temp_inlink;
    node *first_page;

    for (int i = 0; i < npages; i++) {
        p_outlinks[i] = temp_node->page->noutlinks; //get page's inlinks
        if (temp_node->page->inlinks != NULL) {
            first_page = list_copy->page->inlinks;
            //count length of page's total inlinks
            while (list_copy->page->inlinks->next != NULL) {
                list_copy->page->inlinks = list_copy->page->inlinks->next;
                inlinks_count            = inlinks_count + 1;
            }
            list_copy->page->inlinks = first_page;
            inlinks_count            = inlinks_count + 1;
            p_inlinks[i]             = inlinks_count;

            inlinks_data[i] = malloc(p_inlinks[i] * sizeof(int));
            temp_inlink     = temp_node->page->inlinks;

            for (int j = 0; j < p_inlinks[i]; j++) {
                inlinks_data[i][j] = temp_inlink->page->index;
                temp_inlink        = temp_inlink->next;
            }
        } else {
            p_inlinks[i] = 0;
        }
        inlinks_count = 0;

        list_copy = list_copy->next;
        temp_node = temp_node->next;
    }

    // initialise ranks to 1/N (of pages)
    for (int i = 0; i < npages; i++) {
        currRank[i] = prevRank[i] = zeroth_vect;
    }


    //sqrt(vectNorm) > EPSILON same as vectNorm > EPSILON^2
    while (vectNorm > squaredEpsilon) {
        // calculate each page rank
        // for (int i = 0; i < npages; i++) {
        //     vectNorm = 0;
        //     if (p_inlinks[i] 	== 0) {
        //         for (int j = 0; j < p_inlinks[i]; j++) {
        //             vectNorm = vectNorm + prevRank[inlinks_data[i][j]]/(p_outlinks[inlinks_data[i][j]]);
        //         }
        //     }
        //     //assigning currRank[i] to newest t row result
        //     currRank[i] = leading_term + dampener * vectNorm;
        // }

        for (int i = 0; i < nthreads; i++){
          calc[i] = (vect_struct){
          .vectNorm = 0, //result
          .p_inlinks = p_inlinks,
          .p_outlinks = p_outlinks,
          .currRank = currRank,
          .inlinks_data = inlinks_data,
          .leading_term = leading_term,
          .prevRank = prevRank,
          .dampener = dampener,
          .start = i*chunk,
          .end = (i*chunk) + chunk,
          .j = 0,
          };
        }

        for (int i = 0; i < nthreads; i++){
          pthread_create(t_id + 1, NULL, calc_worker, calc+i);
        }

        for (int i = 0; i < nthreads; ++i){
          pthread_join(t_id[i], NULL);
        }

        vectNorm = 0;


        // produce new vector norm
        for (int i = 0; i < npages; i++) {
            difference_vect = (currRank[i] - prevRank[i]);
            vectNorm        = vectNorm + (difference_vect * difference_vect); //squaring difference
        }

        //swapping current and prev (moving down t row)
        tempRank = currRank;
        currRank = prevRank;
        prevRank = tempRank;     //making prev t-1 (prev[0]->prev[1], cur[1]->cur[2])

    }

    temp_node = list;
    for (int i = 0; i < npages; i++) {
        printf("%s %.4lf\n", temp_node->page->name, prevRank[i]);
        temp_node = temp_node->next;
    }
}

/*
 ######################################
 ### DO NOT MODIFY BELOW THIS POINT ###
 ######################################
 */

int main(int argc, char * *argv) {

    /*
     ######################################################
     ### DO NOT MODIFY THE MAIN FUNCTION OR HEADER FILE ###
     ######################################################
     */

    config conf;

    init(&conf, argc, argv);

    node   *list    = conf.list;
    int    npages   = conf.npages;
    int    nedges   = conf.nedges;
    int    nthreads = conf.nthreads;
    double dampener = conf.dampener;

    pagerank(list, npages, nedges, nthreads, dampener);

    release(list);

    return 0;
}
