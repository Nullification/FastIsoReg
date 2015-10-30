/* --- BiPAVA ---
 *  Bidirectional Implementation of Pool Adjacent Violators Algorithm
 *  By Wail Alkowaileet
 *  Center for Complex Engineering Systems at KACST and MIT
 *  w.alkowaileet(at)cces-kacst-mit.org
 *  FastIsoReg is a modification to the previous implentation to accelorate the process
 *  of the Isotonic Regression.
 *  Copyright (C) 2015  Wail Alkowaileet
 * ---
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2003	The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

#include <R.h>
#include "modreg.h"

typedef struct list list;

//LinkedList
struct list
{
    list *next;
    double avg;
    double sum;
    int counter;
    
    int blockStartIndex;
    int blockEndIndex;
    int obselete; //deprecated
};



//Is the element in current block range ?
int checkBoundires(list *node, int index)
{
    return (node->blockStartIndex <= index && node->blockEndIndex >= index);
}

//Get the value that is recent (whereas it's stored in the Array or the list)
double getRecentValue(double *x, list *block, int i)
{
    double retVal = x[i];
    
    if(block != NULL && checkBoundires(block, i))
        retVal = block->avg;
    
    return retVal;
}

//[x]->[y]->[z] After inserting [a] will be: [a]->[x]->[y]->[z]
list* insertFirst(list *head, int index)
{
    list *newHead =(list*)malloc(sizeof(list));
    newHead->next = head;
    
    newHead->avg = 0;
    newHead->sum = 0;
    newHead->counter = 0;
    
    newHead->blockEndIndex = newHead->blockStartIndex = index;
    return newHead;
}

//Merge two adjacent blocks into one block by:
//  1) Computing the new average.
//  2) combine the ranges.
list* merge(list *block1, list *block2)
{
    block1->blockEndIndex = block2->blockEndIndex;
    block1->counter += block2->counter;
    block1->sum += block2->sum;
    block1->avg = block1->sum/block1->counter;
    
    block1->next = block2->next;
    free(block2);
    
}


//
list* addBackwardValue(list* node, double value, int startIndex)
{
    node->sum += value;
    node->counter++;
    node->avg = node->sum/node->counter;
    
    
    node->blockStartIndex = startIndex;
    return node;
}

list* addForwardValue(list* node, double value, int endIndex)
{
    
    if(node->next != NULL && checkBoundires(node->next, endIndex))
    {
        merge(node, node->next);
    }
    else
    {
        node->sum += value;
        node->counter++;
        node->avg = node->sum/node->counter;
        
        node->blockEndIndex = endIndex;
    }
    
    return node;
}



void deleteList(list *head)
{
    list *node = head;
    list *temp = NULL;
    
    while(node != NULL)
    {
        temp = node;
        node = node->next;
        free(temp);
    }
}

double *getNewValues(double *x, int n, list *blocks, double *y)
{
    int i=0, j;
    list *node = blocks;
    while(i<n)
    {
        if(node != NULL && checkBoundires(node, i))
        {
            for(j=i;j<=node->blockEndIndex;j++)
                y[j] = node->avg;
            
            i = node->blockEndIndex+1;
            node = node->next;
        }
        else
        {
            y[i] = x[i];
            i++;
        }
    }
    
    
    deleteList(blocks);
    
    return y;
}

// Faster implementation of isoreg
double *PAV(double *x, int n, double *y)
{
    //Definitions;
    const double DOUBLE_MAX = 1e+200;
    const int false = 0;
    const int true = 1;
    
    list *blocks = NULL; // list of blocks
    int i, j, k, blockStartIndex;
    int doAvg = false;
    int forward = false;
    
    
    
    
    //start from the last element of the array
    i=n-1;

    
    
    while(i>=1)
    {
        if(getRecentValue(x, blocks, i) > getRecentValue(x, blocks, i-1) && !doAvg)
        {
            i--;
            continue;
        }
        
        //There is a violation        
        

        if(!doAvg)
        {
            blocks = insertFirst(blocks, i);
            //Do Backward updating
            blockStartIndex = i;
            doAvg = true; //Start updating the avgs.
            k=i-1; //Violater index
        }
        
        
        if(forward)
        {
            addForwardValue(blocks, x[i], i);
            forward = false;
            blockStartIndex = blocks->blockEndIndex;
        }
        else
            addBackwardValue(blocks, x[i], i);
        
        //update backward
        for(j=k;j>=0 && doAvg;j--)
        {
            if(blocks->avg < x[j] )
            {
                addBackwardValue(blocks, x[j], j);
            }
            else
                doAvg = false;
        }
        
        //Move to the first element that has NOT been visited yet
        i=j+1;
        
        if (blockStartIndex+1>=n || !(doAvg = getRecentValue(x, blocks->next, blockStartIndex+1) < blocks->avg))
        {
            continue;
        }
        
        //After updating backward, a violation appeared at the end of the block
        //Go there and fix it
        k=blocks->blockStartIndex-1; //This will hold the first index to look at after finishing the forward update
        i=(++blockStartIndex); // First element after the block range
        forward = true; // Run forward update
        
        
    }
    
    return getNewValues(x, n, blocks, y);
    
}



SEXP isoreg(SEXP y)
{
    int n = LENGTH(y), i, ip, known, n_ip;
    double tmp, slope;
    SEXP yc, yf, iKnots, ans;
    const char *anms[] = {"y", "yf", ""};
    
    double *yPAV;

    /* unneeded: y = coerceVector(y, REALSXP); */

    PROTECT(ans = mkNamed(VECSXP, anms));

    SET_VECTOR_ELT(ans, 0, y);
    SET_VECTOR_ELT(ans, 1, yf = allocVector(REALSXP, n));
    
    PAV(REAL(y), n, REAL(yf));
    
    
    
    
    UNPROTECT(1);
    return(ans);
}
