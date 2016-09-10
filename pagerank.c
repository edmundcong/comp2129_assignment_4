#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <immintrin.h>

#include "pagerank.h"



void pagerank(node *list, int npages, int nedges, int nthreads, double dampener) {
    double vectNorm       = EPSILON;
    double squaredEpsilon = EPSILON*EPSILON;
    double difference_vect;
    double *tempRank;
    double leading_term   = (1.0 - dampener)/npages;
    double initial_p = 1.0 / npages;

    double *current_rank = malloc(npages*sizeof(double));
    double *prev_rank = malloc(npages*sizeof(double));

    int * *inlinks_data = (int * *)malloc(npages*sizeof(int *));
    int *p_inlinks      = malloc(npages*sizeof(int));
    int *p_outlinks     = malloc(npages*sizeof(int));

    int inlinks_count   = 0;

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
        list_copy     = list_copy->next;
        temp_node     = temp_node->next;
    }

    // initialise ranks to 1/N (of pages)
    for (int i = 0; i < npages; i++) {
        current_rank[i] = prev_rank[i] = initial_p;
    }


    for (int i = 0; i < nthreads; i++)

        //sqrt(vectNorm) > EPSILON same as vectNorm > EPSILON^2
        while (vectNorm > squaredEpsilon) {
            // calculate each page rank
            for (int i = 0; i < npages; i++) {
                vectNorm = 0;
                if (p_inlinks[i] != 0) {
                    for (int j = 0; j < p_inlinks[i]; j++) {
                        vectNorm = vectNorm + prev_rank[inlinks_data[i][j]]/(p_outlinks[inlinks_data[i][j]]);
                    }
                }
                //assigning current_rank[i] to newest t row result
                current_rank[i] = leading_term + dampener * vectNorm;
            }

            vectNorm = 0;

            // produce new vector norm
            for (int i = 0; i < npages; i++) {
                difference_vect = (current_rank[i] - prev_rank[i]);
                vectNorm        = vectNorm + (difference_vect * difference_vect); //squaring difference
            }

            //swapping current and prev (moving down t row)
            tempRank = current_rank;
            current_rank = prev_rank;
            prev_rank = tempRank; //making prev t-1 (prev[0]->prev[1], cur[1]->cur[2])

        }

    temp_node = list;
    for (int i = 0; i < npages; i++) {
        printf("%s %.4lf\n", temp_node->page->name, prev_rank[i]);
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
