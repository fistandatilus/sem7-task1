#ifndef TABLE_PRINT_H
#define TABLE_PRINT_H

#include <cstdio>


void print_header(const char *title, char **headers, int width, FILE *fp=stdin);
void print_row(const double *data, int width, int height=1, FILE *fp=stdin);


#endif //TABLE_PRINT_H
