#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <immintrin.h>

#include "pagerank.h"


void pagerank(node *list, int npages, int nedges, int nthreads, double dampener) {
    struct node *list_copy = malloc(sizeof(struct node));
    // struct node *temp = malloc(sizeof(struct node));
    node node_array[npages];

    memcpy(list_copy, list, sizeof(node));


    double p_scores = (1.0/npages);
    double p_array[npages*2];
    double temp_result = 0;
    int    i           = 0;

    /*+1.0 to started off with (will be modified after first iteration).
       just so it can pass the initial while loop condition*/
    double vect_norm = EPSILON + 1.0;

    for (int i = 0; i < npages; i++) {
        p_array[i] = p_scores;
    }
    // printf("list: %s\n", list->page->inlinks->page->name);
    // printf("list_copy: %s\n", list_copy->page->inlinks->page->name);
    // printf("list_copy: %s\n", list_copy->next->page->name);
    while (vect_norm > EPSILON) {
        i = 0;
        for (int j = npages; j < npages*2; j++) {
            node_array[i] = *list;
            i++;
            while (list->page->inlinks != NULL) {
                temp_result += p_array[j-npages]/list->page->inlinks->page->noutlinks;
                printf("%s<-%s: %f\n", list->page->name, list->page->inlinks->page->name, temp_result);
                if (list->page->inlinks->next != NULL) {
                    list->page->inlinks = list->page->inlinks->next;
                } else {
                    break;
                }
            }
            if (list->page->inlinks == NULL) {
                p_array[j] = (1-dampener)/npages;
                break;
            }
            p_array[j] = (temp_result*dampener) + (1-dampener)/npages;
            if (list->next != NULL) {
                list = list->next;
            }
            temp_result = 0;
        }
        //need to reset the list so it doesn't stay at the end (D)
        //output should b
        printf("list_copy: %s\n", list_copy->page->inlinks->page->name);
        list = list_copy;

        vect_norm = 0;

        for (int j = 0; j < npages; j++) {
            vect_norm += (p_array[j+npages] - p_array[j])*(p_array[j+npages] - p_array[j]);
            printf("p_ar[%d]: %f, p_ar[%d]: %f, v_norm %f\n", j+npages, p_array[j+npages], j, p_array[j], vect_norm);
            p_array[j] = p_array[j+npages];
        }
        printf("\n");
        vect_norm = sqrt(vect_norm);

        printf("vect_norm %.3lf\n\n", vect_norm);
    }

        for (int i = 0; i < npages; i++) {
            printf("node_array[%d]: %s\n", i, node_array[i].page->name);
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
