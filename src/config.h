/************************************************************************** 
  Configuration file for kc. 
  
  (C) Copyright, 1992-1996 by 
  Kenneth Geisshirt (kneth@fatou.ruc.dk)          Keld Nielsen (kn@kiku.dk)
  Dept. of Life Sciences and Chemistry       Dept. of Theoretical Chemistry
  Roskilde University                              University of Copenhagen
  P.O. Box 260                                         Universitetsparken 5
  4000 Roskilde                                        2100 Copenhagen East
  Denmark                                                           Denmark

  All platform dependent details are placed in this file. The platform 
  is communicated through a macro (_PLATFORM_*_).

  See kc.tex for details. 

  Last updated: 4 January 1996 KG
**************************************************************************/

#ifndef _CONFIG_
#define _CONFIG_

#define VERSION "1.10"
#define STRING_LENGTH 25

#ifdef _PLATFORM_DOS_
#define MALLOCTYPE void
#endif

#ifdef _PLATFORM_GCC_
#define MALLOCTYPE void
#endif

#ifdef _PLATFORM_HPUX_
#define MALLOCTYPE char
#endif

#ifdef _PLATFORM_CONVEX_
#define MALLOCTYPE void
#define MAXFLOAT 1.0e38
#endif

#ifdef _PLATFORM_ULTRIX_
#define MALLOCTYPE char
#endif

#ifdef _PLATFORM_SGI_
#define MALLOCTYPE void
#endif

#ifdef _PLATFORM_AIX_
#define MALLOCTYPE void
#endif

#endif /* _CONFIG_ */
