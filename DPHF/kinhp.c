void init()
{  printf("\033h\033J");
}

void disp(n, t, h, q, x, xmax, xmin) 
int n; 
double t,h,q; 
double x[],xmax[],xmin[];
{  int i,row,col;
   printf("\033\ht = %e",t);
   printf("     h = %e",h);
   printf("     q = %e",q);
   printf("\n");
   for (i= 0; i<n ;i++)
   {  printf("\n%2d:",i);
      printf("   %e",x[i]);
      printf("   %e",xmax[i]);
      printf("   %e",xmin[i]);
   }
}
