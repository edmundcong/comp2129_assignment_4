#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <immintrin.h>

#include "pagerank.h"


void pagerank(node *list, int npages, int nedges, int nthreads, double dampener) {
    double p_results[npages][npages]; //stores values for P array
    // double p_converge[npages]; //converge results {-0.036, etc}
    double convergence = EPSILON;
    double con_temp    = 0.0;
    double first_term  = (1.0-dampener)/npages;
    double temp_val    = 0.0; //temp value for p_results;

    int i = 1;
    int j = 0;

    node *list_ptr    = list; //operate on list_ptr
    node *list_inlinks;

    for (int i = 0; i < npages; i++) { //set 0, 0->n values to 1/npages
        p_results[0][i] = 1.0/npages;
    }

    printf("%s: {%s}\n", list_ptr->page->name, list_ptr->page->inlinks->page->name);
    printf("%s->%s\n", list->page->name, list->page->inlinks->page->name);

    while (convergence >= EPSILON) {
        while (list_ptr->next != NULL) {
            //only performs this loop on the first iteration
            while (list_ptr->page->inlinks != NULL) {
                temp_val += (p_results[i-1][j] //calculating each P
                             /list_ptr->page->inlinks->page->noutlinks);
                if (list_ptr->page->inlinks->next != NULL) { //traver list
                    list_ptr->page->inlinks = list_ptr->page->inlinks->next;
                } else {
                    break;
                }
            }
            j++;
            p_results[i][j] = (temp_val*dampener) + first_term;

            list_ptr = list_ptr->next;
            temp_val = 0;

            if (list_ptr->page->inlinks == NULL) p_results[i][j+1] = first_term;
        }


        for (int j = 1; j <= npages; j++) {
            con_temp += (p_results[i][j]-p_results[i-1][j-1])*(p_results[i][j]-p_results[i-1][j-1]);
        }

        convergence = sqrt(con_temp);
        printf("%f\n", convergence);
        j = 0;
        i++;
        // list_ptr = *list_head;


        printf("%s: {%s}\n", list_ptr->page->name, list_ptr->page->inlinks->page->name);
        printf("%s: {%s}\n", list_ptr->next->page->name, list_ptr->page->inlinks->page->name);


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
