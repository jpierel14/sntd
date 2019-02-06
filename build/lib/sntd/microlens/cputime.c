/*
 * Fake the convex call CPUTIME
 */
float
cputime_(time)
double *time;
{
   long clock();
   static int first = 1;

   if(first) {
      first = 0;
      (void)clock();
   }
   
   return(1e-6*clock() - *time);
}
