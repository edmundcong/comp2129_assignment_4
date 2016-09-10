#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <immintrin.h>

#include "pagerank.h"


void pagerank(node *list, int npages, int nedges, int nthreads, double dampener) {
    double *rankOne = malloc(npages*sizeof(double));
    double *rankTwo = malloc(npages*sizeof(double));

    double vectNorm       = EPSILON;
    double squaredEpsilon = EPSILON*EPSILON;
    double calcDampener   = (1 - dampener)/npages;
    double minusedVectCalc;

    double *currRank = rankOne;
    double *prevRank = rankTwo;
    double *tempRank;

    int * *inlinksData = (int * *)malloc(npages*sizeof(int *));
    int *numInlinks    = malloc(npages*sizeof(int));
    int *numOutlinks   = malloc(npages*sizeof(int));
    int inlinks_count  = 0;


// populate arrays for storing page data
    node *tempNode = list;
    node *tempLinkNode;
    for (int i = 0; i < npages; i++) {
        numOutlinks[i] = tempNode->page->noutlinks;

        if (tempNode->page->inlinks != NULL) {
            numInlinks[i]  = inlinks_count;
            inlinksData[i] = malloc(numInlinks[i] * sizeof(int));
            tempLinkNode   = tempNode->page->inlinks;
            for (int j = 0; j <= numInlinks[i]; j++) {
                inlinksData[i][j] = tempLinkNode->page->index;
                printf("{%s}\n", tempLinkNode->page->name);
                if (tempLinkNode->next == NULL) break;
                printf("~\n");
                tempLinkNode = tempLinkNode->next;
                printf("%d\n", numInlinks[i]);
                numInlinks[i]++;

            }

        } else {
            numInlinks[i] = 0;
        }
        tempNode      = tempNode->next;
    }

// initialise ranks to 1/N (of pages)
    for (int i = 0; i < npages; i++) {
        rankOne[i] = 1.0/npages;
        rankTwo[i] = 1.0/npages;
    }

    while (vectNorm > squaredEpsilon) {

        // calculate each page rank
        for (int i = 0; i < npages; i++) {
            vectNorm = 0;
            if (numInlinks[i] != 0) {
                for (int j = 0; j < numInlinks[i]; j++) {
                    vectNorm = vectNorm + prevRank[inlinksData[i][j]]/(numOutlinks[inlinksData[i][j]]);
                }
            }
            currRank[i] = calcDampener + dampener * vectNorm;
        }

        vectNorm = 0;

        // produce new vector norm
        for (int i = 0; i < npages; i++) {
            minusedVectCalc = (currRank[i] - prevRank[i]);
            vectNorm        = vectNorm + (minusedVectCalc * minusedVectCalc);
        }

        tempRank = currRank;
        currRank = prevRank;
        prevRank = tempRank;

    }

    tempNode = list;
    for (int i = 0; i < npages; i++) {
        printf("%s %.4lf\n", tempNode->page->name, prevRank[i]);
        tempNode = tempNode->next;
    }

    free(rankOne);
    free(rankTwo);
    for (int i = 0; i < npages; i++) {
        if (numInlinks[i] != 0) {
            free(inlinksData[i]);
        }
    }
    free(inlinksData);
    free(numInlinks);
    free(numOutlinks);
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
