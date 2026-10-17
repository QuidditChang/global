/* Read-only population accounting. Bins 0..8 are exact; bin 9 means >=9.
 * Layout: histogram[10], previous-to-current transitions[100], then
 * deficient cells and additions required for thresholds 2,4,8 (six values). */
#ifndef PICES_COVERAGE_STATS_H
#define PICES_COVERAGE_STATS_H
#define PICES_COVERAGE_STATS_SIZE 116
static void pices_coverage_count(long long *s,int count,int previous,int valid)
{
    const int threshold[3]={2,4,8};
    int bin=count<9?count:9,i;
    s[bin]++;
    if(valid)s[10+(previous<9?previous:9)*10+bin]++;
    for(i=0;i<3;i++)if(count<threshold[i]) {
        s[110+2*i]++;
        s[111+2*i]+=threshold[i]-count;
    }
}
#endif
