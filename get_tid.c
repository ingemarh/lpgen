#include <time.h>
#include <stdio.h>
double get_tid (void);
void print_tid (double diff_time, int nth, char *decoder);
void print_tid_clutter (double diff_time, int ncth, int nlth, char *decoder);

double
get_tid ()
{
  struct timespec tp;
  clockid_t clk_id;

  clk_id = CLOCK_MONOTONIC;
  clock_gettime (clk_id, &tp);
  return (double) (tp.tv_sec + tp.tv_nsec / 1.e9);
}

void
print_tid (double diff_time, int nth, char *decoder)
{
  if (diff_time < 1.e-3)
    printf ("%s: %d threads: %d us", decoder, nth, (int) (diff_time * 1.e6));
  else if ((diff_time < 1.e-1) && (diff_time > 1.e-3))
    printf ("%s: %d threads: %d ms", decoder, nth, (int) (diff_time * 1.e3));
  else
    printf ("%s: %d threads: %.3lf s", decoder, nth, diff_time);
}

void
print_tid_clutter (double diff_time, int ncth, int nlth, char *decoder)
{
  if (diff_time < 1.e-3)
    printf ("%s: [%d %d] threads: %d us", decoder, ncth, nlth,
	    (int) (diff_time * 1.e6));
  else if ((diff_time < 1.e-1) && (diff_time > 1.e-3))
    printf ("%s: [%d %d] threads: %d ms", decoder, ncth, nlth,
	    (int) (diff_time * 1.e3));
  else
    printf ("%s: [%d %d] threads: %.3lf s", decoder, ncth, nlth, diff_time);
}
