
#include "STRING2.h"

/**********************************************************************************************/

extern int Indx_Char(char *str, char c, int order){
    
    /******************************************************************************************
     Purpose:
     Find the position of the character in a string.
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     str, the input string.
     c, the charaster.
     order, the order of the charaster.
     return:
     return the position of the charaster.
     ******************************************************************************************/
    
    if (str == NULL || *str == '\0')
        return -1;
    if (order <= 0)
        return -2;
    
    char *p = str;
    int indx=0 ,num =0;
    while (*p != '\0' && num!=order ) {
        if (*p == c)
            num++;
        p++;
        indx++;
    }
    if(strlen(str) > indx){
        return indx;
    }else{
        return -3;
    }
}


extern void Trim1(char *str){
    
    /******************************************************************************************
     Purpose:
     Remove leading spaces in a string.
     Record of revisions:
     29 Nov. 2019
     Input parameters:
     str, the input string.
     leading and trailing spaces.
     Output parameters:
     str, the input string.
     ******************************************************************************************/
    
    if (str == NULL || *str == '\0')
        return;
    
    int len = 0;
    char *p = str;
    
    while (*p != '\0' && isspace(*p))
    {
        p++;
        len++;
    }
    
    memmove(str, p, strlen(str) - len + 1);
    
    return;
}


extern void Trim2(char *str){
    
    /******************************************************************************************
     Purpose:
     Remove trailing spaces in a string.
     Record of revisions:
     29 Nov. 2019
     Input parameters:
     str, the input string.
     leading and trailing spaces.
     Output parameters:
     str, the input string.
     ******************************************************************************************/
    
    if (str == NULL || *str == '\0')
        return;
    
    long len = strlen(str);
    char *p = str + len - 1;
    while (p >= str  && isspace(*p))
    {
        *p = '\0';
        p--;
    }
    
    return;
}


extern void Trim(char *str, int indx){
    
    /******************************************************************************************
     Purpose:
     Remove spaces in a string.
     Record of revisions:
     29 Nov. 2019
     Input parameters:
     str, the input string.
     indx, indx = 1, for the leading spaces; indx =2, for the trailing spaces; indx=3 for both
        leading and trailing spaces.
     Output parameters:
     str, the input string.
     ******************************************************************************************/
    
    switch (indx) {
        case 1:
            Trim1(str);
            break;
            
        case 2:
            Trim2(str);
            break;
            
        case 3:
            Trim1(str);
            Trim2(str);
            break;
            
        default:
            break;
    }
    
    return;
}


extern void String_Copy(char *str1, char *str2, long len, int trim_flag){
    
    /******************************************************************************************
     Purpose:
     Copy a string from str2 to str1.
     Record of revisions:
     29 Nov. 2019
     Input parameters:
     str2, the input string.
     len, the length.
     trim_flag, if the flag > 0, the leading and trailing spaces in str1 will be removed.
     Input parameters:
     str1, the output string.
     ******************************************************************************************/
    
    memmove(str1, str2, len);
    str1[len] = '\0';
    
    if (trim_flag > 0) {
        Trim(str1, 3);
    }
    
    return;
}

/**********************************************************************************************/

