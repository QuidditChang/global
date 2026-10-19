/* Isolated replay instrumentation. Never changes the production fingerprint. */
#define REPLAY_DIAG_BLOCKS 96
static FILE *replay_diag_file;
static char replay_diag_base[2048];
static struct { const char *name; PicesSHA sha; size_t bytes, calls; } replay_diag_blocks[REPLAY_DIAG_BLOCKS];
static int replay_diag_count;
static PicesSHA replay_diag_zero_legacy;
static void replay_diag_begin(struct All_variables *E) {
 char path[2100];
 snprintf(replay_diag_base,sizeof(replay_diag_base),"%s.checkpoint_diag.%d",E->control.data_file,E->parallel.me);
 snprintf(path,sizeof(path),"%s.tsv",replay_diag_base);
 replay_diag_file=fopen(path,"w");
 if(!replay_diag_file)pices_fail(E,"cannot open replay checkpoint diagnostic");
 replay_diag_count=0;pices_sha_init(&replay_diag_zero_legacy);
 fprintf(replay_diag_file,"scope\trunning_state_only\nrank\t%d\nstep\t%d\n",E->parallel.me,E->monitor.solution_cycles);
}
static void replay_diag_scalar(const char *name,int index,double value) {
 if(replay_diag_file)fprintf(replay_diag_file,"scalar\t%s\t%d\t%.17g\t%a\n",name,index,value,value);
}
static void replay_diag_add(PicesSHA *s,const void *data,size_t bytes,const char *name) {
 int i;
 pices_sha_add(s,data,bytes); /* Original bytes/order, unconditionally. */
 if(!replay_diag_file)return;
 for(i=0;i<replay_diag_count;i++)if(!strcmp(name,replay_diag_blocks[i].name))break;
 if(i==replay_diag_count) {
  if(i==REPLAY_DIAG_BLOCKS){fprintf(stderr,"too many replay diagnostic blocks\n");MPI_Abort(MPI_COMM_WORLD,72);}
  replay_diag_count++;replay_diag_blocks[i].name=name;
  replay_diag_blocks[i].bytes=replay_diag_blocks[i].calls=0;pices_sha_init(&replay_diag_blocks[i].sha);
 }
 pices_sha_add(&replay_diag_blocks[i].sha,data,bytes);
 replay_diag_blocks[i].bytes+=bytes;replay_diag_blocks[i].calls++;
 if(!strcmp(name,"coupled")) {
  double copy[29];
  if(bytes!=sizeof(copy)){fprintf(stderr,"unexpected coupled diagnostic layout\n");MPI_Abort(MPI_COMM_WORLD,72);}
  memcpy(copy,data,bytes);copy[18]=copy[20]=copy[21]=copy[22]=copy[23]=0;
  pices_sha_add(&replay_diag_zero_legacy,copy,bytes);
 } else pices_sha_add(&replay_diag_zero_legacy,data,bytes);
}
static void replay_diag_forcing(const char *path,const char *sha) {
 if(replay_diag_file)fprintf(replay_diag_file,"forcing\t%s\t%s\n",path,sha);
}
static void replay_diag_finish(struct All_variables *E,const char *actual,const char *expected,int mismatch) {
 char path[2100],sha[65];FILE *f;size_t offset=0;int i;
 for(i=0;i<replay_diag_count;i++) {
  pices_sha_end(&replay_diag_blocks[i].sha,sha);
  fprintf(replay_diag_file,"block\t%s\t%zu\t%zu\t%s\n",replay_diag_blocks[i].name,replay_diag_blocks[i].calls,replay_diag_blocks[i].bytes,sha);
 }
 pices_sha_end(&replay_diag_zero_legacy,sha);
 fprintf(replay_diag_file,"candidate_zero_unused_fields_sha256\t%s\n",sha);
 while(actual[offset]&&expected[offset]&&actual[offset]==expected[offset])offset++;
 fprintf(replay_diag_file,"metadata_equal\t%d\nfirst_difference_offset\t%zu\n",!mismatch,offset);
 if(ferror(replay_diag_file)||fclose(replay_diag_file))pices_fail(E,"cannot close replay diagnostic");
 replay_diag_file=NULL;
 snprintf(path,sizeof(path),"%s.actual.json",replay_diag_base);f=fopen(path,"w");
 if(!f)pices_fail(E,"cannot open saved metadata diagnostic");
 if(fputs(actual,f)<0||fclose(f))pices_fail(E,"cannot write saved metadata diagnostic");
 snprintf(path,sizeof(path),"%s.expected.json",replay_diag_base);f=fopen(path,"w");
 if(!f)pices_fail(E,"cannot open expected metadata diagnostic");
 if(fputs(expected,f)<0||fclose(f))pices_fail(E,"cannot write expected metadata diagnostic");
}
