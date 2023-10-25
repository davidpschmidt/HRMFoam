#define GCC
#if defined (GCC) 
#define FORTNAME(A) A##_
#else 
#define FORTNAME(A) A##__
#endif
