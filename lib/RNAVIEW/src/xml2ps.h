 /* #define PSPIONT 12    the letter size for postscript file */ 

extern void read_sequence(char *inpfile, char *resname, long *author_seq, long *nres);

extern void extract_sequence(FILE *inp, char *resname, long *nres);

extern void read_xy_coord(char *inpfile, double **xy, long *resid_idx, long *num_xy);



extern void read_pair_type(char *inpfile,char **pair_type,long **npair_idx,long *npair,

                    long *nhelix, long **helix_idx, long *helix_length,

                    long *nsing, long *sing_st, long *sing_end);

extern void get_xyz_coord(FILE *inp, double *x, double *y, double *z);

extern void read_O3prime_P_xyz(char *inpfile, double **o3_prime_xyz, double **p_xyz, long *npo3);

extern void get_chain(long nres, double **a, double **b, long *nchain,long **chain_idx);

extern void link_chain(long nchain,long **chain_idx,  double **xy, long *broken);

extern void label_ps_resname(long num_res, char *resname,  double **xy, long *sugar_syn);

extern void label_5p_3p(FILE *psfile, long i, long **chain_idx, double **xy);

extern void write_5p_3p(long k1, long k2, double a, double **xy, char *labelP);

extern void label_seq_number(long nres, long nhelix, long **helix_idx,

                      long *helix_length,long nsing, long *sing_st,

                      long *sing_end, double **xy, long *author_seq);

extern void label_seq(long k1, long k2, double a, double **xy, long key, long *author_seq);

extern void draw_LW_diagram(long npair, char **pair_type, char *resname,

                     long **npair_idx, double **xy);

extern void extract_author_seq(FILE *inp, long *author_seq, long *nseq);



extern double h_width(long nhelix, long **helix_idx, long *helix_length, double **xy);







extern double slope(long k1, long k2,  double **xy);

extern void get_value(FILE *inp, char *value);

extern double dist(double *a,  double *b,  long n);

extern void usage(void);

extern void element_in_bracket(FILE *inp,char *item,long *size,char *lett,long *key);

extern void get_num_residue(long size_item, char *lett, long *nres);


extern void get_base_pair(long *num_pair, long **base_pair, char **edge_type,

                   long **pair_idx, long *num_helix, long **helix_st_end,

                   long *helix_len);

extern void get_ss_xy(double **xy, long *num_xy);

extern void make_a_line(char *str);

extern void search_item(long *num_in, char lett[], char *item1, char *item2);

extern void get_model_id(char *lett, char *identifer, char *item, long *model_id);

extern long num_of_pair(char *inpfile);

extern void get_residue_num(char *str, long *nres1, long *nres2, long *seq);

extern void read_bs_pair(char *inpfile, long *npair, char *edge_type, char *cis_tran,

                  char *resname, long *chain_id, long *seq, long **num_idx);





extern void LW_shapes(char *bseq, long k1, long k2, char *pair_type, double *x, double *y, double at, double ratio);



extern void nrerror(char error_text[]);





/*******************/







extern void generate_ps_file(long num_res, char *resname,  double **xy);

extern void ps_head(long *bbox);

extern double xml_dmax(double a, double b);

extern double xml_dmin(double a, double b);

extern void xml_max_dmatrix(double **d, long nr, long nc, double *maxdm);

extern void xml_min_dmatrix(double **d, long nr, long nc, double *mindm);

extern void xml_move_position(double **d, long nr, long nc, double *mpos);

extern long xml_round(double d);

extern void xml_xy4ps(long num, double **oxy, double ps_size, long n);



extern void get_position(FILE *inp, long *position);

extern void get_xy_position(FILE *inp, double *x, double *y);

/*******************/



extern void twopoints(double *xy0, double a, double d, double *xy1, double *xy2);



extern void line(FILE *psfile, double *xy1, double *xy2);

extern void square(FILE *psfile, long fill, double *xy1, double *xy2, double d, double a);

extern void circle(FILE *psfile, long fill, double *xy1, double *xy2, double r);

extern void triangle(FILE *psfile, long fill, double *xy1, double *xy3, double d,double a);

extern void color_line(FILE *psfile, char *color, double *xy1, double *xy2);

extern void dashline_red(FILE *psfile, double *xy1, double *xy2);

extern void dashline(FILE *psfile, double *xy1, double *xy2);

extern void dotline(FILE *psfile, double *xy1, double *xy2);



extern void del_extension(char *pdbfile, char *parfile);

