#ifndef __TEXTIO_H__
#define __TEXTIO_H__

#define TEXTIO_MAXBUFFER 64000
#define TEXTIO_COLUMNMAJOR +1
#define TEXTIO_ROWMAJOR -1
#define TEXTIO_COMMENT_CHAR '#'

/* 
  CSV file when sep = ','
  TSV file when sep = '\t' 
*/
int textio_write_double_matrix(const char *filename,
                               const double *a,
												  		 int m,
															 int n,
															 const char *formatspec,
															 char sep,
															 const char* head_comment,
															 const char* foot_comment)
{
	char default_formatspec[] = "%.16e";
	char str1[8];
	char str2[8];
	if (formatspec != NULL) {
		sprintf(str1, "%s%c", formatspec, sep);
		sprintf(str2, "%s\n", formatspec);
	} else {
		sprintf(str1, "%s%c", default_formatspec, sep);
		sprintf(str2, "%s\n", default_formatspec);
	}
	FILE *pfile = NULL;
	pfile = fopen(filename, "w");
	if (!pfile) return 0;
	if (head_comment != NULL) {
		fprintf(pfile, "%c %s\n", TEXTIO_COMMENT_CHAR, head_comment);
	}
	int e = 0;
	for (int ii = 0; ii < m; ii++) {
		for (int jj = 0; jj < (n - 1); jj++) {
			fprintf(pfile, str1, a[jj * m + ii]);
			e++;
		}
		fprintf(pfile, str2, a[(n - 1) * m + ii]);
		e++;
	}
	if (foot_comment != NULL) {
		fprintf(pfile, "%c %s\n", TEXTIO_COMMENT_CHAR, foot_comment);
	}
	fclose(pfile);
	return e;
}

/* Read a matrix from a textfile; only count elements if buf==NULL, otherwise store elements in buf */
int textio_read_table_file(const char *filename,
                           int *m,
													 int *n,
													 double *buf,
													 int maxbuf,
													 int colStride,
													 int rowStride,
													 char sep,
													 bool detectComments)
{
	static char linebuffer[TEXTIO_MAXBUFFER];
	char delims[2];
	sprintf(delims, "%c", sep);

	*m = 0; 
	*n = 0;
	FILE *fp = NULL;

	fp = fopen(filename, "r");
	if (!fp) return 0;

	int ret,idx;
	int totalcounter = 0;
	int rowcounter = 0;
	int columncounter = 0;
	int firstcolumncount = 0;
	double dmy;
	char *tmp,*pch;
	
	while (true) {
		tmp = fgets(linebuffer, TEXTIO_MAXBUFFER, fp);
		if (tmp == NULL) break;
		if (detectComments && *tmp == TEXTIO_COMMENT_CHAR) continue;
		rowcounter++;
		columncounter = 0;
		pch = strtok(tmp, delims);
		while (pch != NULL) {
			ret = sscanf(pch, "%lf", &dmy);
			if (ret != 1) break;
			columncounter++;
			if (buf != NULL) { /* calculate the index for storage based on colStride,rowStride ints */
				idx = colStride * (columncounter - 1) + rowStride * (rowcounter - 1);
				if (idx < maxbuf && idx >= 0) buf[idx] = dmy;
			}
			totalcounter++;
			pch = strtok(NULL, delims);
		}
		if (rowcounter == 1) {
			firstcolumncount = columncounter;
		}
		if (columncounter != firstcolumncount) {
			fclose(fp);
			return -rowcounter;
		}
	}
	*m = rowcounter;
	*n = firstcolumncount;
	fclose(fp);
	return totalcounter;
}

int textio_read_double_matrix(const char *filename,
                              double **dat,
															int *m,
															int *n,
															int data_order,
															char sep,
															bool detectComments) 
{
	int ret = -1, rows, cols;
	int ret2 = -1, rows2, cols2;
	*dat = NULL;
	ret = textio_read_table_file(filename, 
	                             &rows, 
															 &cols, 
															 NULL, 
															 -1, 
															 -1, 
															 -1, 
															 sep,
															 detectComments);
	if (ret <= 0) {
		return ret;	/* failure code */
	}
	if (ret != rows * cols) {
		return 0;	/* also a failure code */
	}
	*dat = (double*) malloc(sizeof(double) * ret);
	if (*dat == NULL)
	  return 0;
	if (data_order == TEXTIO_COLUMNMAJOR)
	{ /* colStride=rows, rowStride=1 */
		ret2 = textio_read_table_file(filename, 
		                              &rows2, 
																	&cols2, 
																	*dat, 
																	ret, 
																	rows, 
																	1, 
																	sep,
																	detectComments); 
	} 
	else if (data_order == TEXTIO_ROWMAJOR)
	{ /* colStride=1, rowStride=cols */
		ret2 = textio_read_table_file(filename, 
		                              &rows2, 
																	&cols2, 
																	*dat, 
																	ret, 
																	1, 
																	cols, 
																	sep,
																	detectComments);
	}
	if (ret2 != ret || ret2 != rows2 * cols2) 
	  return 0;
	*m = rows2;
	*n = cols2;
	return ret2;
}

#endif
