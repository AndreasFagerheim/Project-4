#ifndef ISING_H
#define ISING_H


#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

using namespace std;


void  **matrix(int, int, int);
void free_matrix(void **);
double ran2(long *idum);
#endif // LIB_H
