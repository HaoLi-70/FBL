
#ifndef STRING2_h
#define STRING2_h

/******************************************************************************************/

#include <stdio.h>
#include <ctype.h>
#include <string.h>

/******************************************************************************************/

extern int Indx_Char(char *str, char c, int order);

extern void Trim1(char *str);

extern void Trim2(char *str);

extern void Trim(char *str, int indx);

extern void String_Copy(char *str1, char *str2, long len, int trim_flag);

/******************************************************************************************/

#endif /* STRING2_h */
