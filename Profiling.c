/* Profiling functions .... return elapsed CPU time etc. 
   These functions seem the most likely to get broken by
   different architectures/operating systems &c */

#if defined(__osf__) || defined(__aix__) || defined(__sunos__) || defined(__sgi)
#include <sys/time.h>
#include <sys/resource.h>
#define RUSAGE_STYLE_TIME
#elif defined(_UNICOS) || defined(__hpux)
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#define TIMES_STYLE_TIME
#endif

/* ===============================================
   Function to return currently elapsed CPU usage
   =============================================== */

double CPU_time()

{
#if defined(RUSAGE_STYLE_TIME)
  struct rusage rusage;
  double time;

  getrusage(RUSAGE_SELF, &rusage);
  time = rusage.ru_utime.tv_sec + 1.0e-6 * rusage.ru_utime.tv_usec;
#elif defined(TIMES_STYLE_TIME)
  struct tms time_now;
  time_t utime;
  long sometime;

  double time;
  static double initial_time;
  static int visit = 0;

  if (visit == 0)
  {
    sometime = times(&time_now);
    initial_time = (double)time_now.tms_utime / (double)CLK_TCK;
    visit++;
  }

  sometime = times(&time_now);
  time = (double)time_now.tms_utime / (double)CLK_TCK - initial_time;
  return (time);

#else /* stupid, break nothing "timer" */
  static double time;

  time += 0.0001;
#endif

  return (time);
}
