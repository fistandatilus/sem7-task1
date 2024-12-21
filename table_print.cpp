#include "table_print.h"

void print_header(const char *title, char **headers, int width, FILE *fp) {
    fprintf(fp, "\n\\begin{tabular}{ |");
    for (int i = 0; i < width; i++) fprintf(fp, "l|");
    fprintf(fp, "}\n\\hline\n\\multicolumn{%d}{|c|}{%s} \\\\\n\\hline\n", width, title);
    for (int i = 0; i < width-1; i++) fprintf(fp, "%s &", headers[i]);
    fprintf(fp, "%s ", headers[width-1]);
    fprintf(fp, "\\\\\n\\hline\n");
}


void print_row(const double *data, int width, int height, FILE *fp) {
    for(int row = 0; row < height; row++) {
        for (int collumn = 1; collumn < width; collumn++) {
            fprintf(fp, "& $%e$ ", data[(collumn-1)*height + row]);
        }
        fprintf(fp, "\\\\\n");
    }
    fprintf(fp, "\\\\\n\\hline\n");
}
