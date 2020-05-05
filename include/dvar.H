#ifndef DVAR_H
#define DVAR_H


#ifdef __cplusplus
#include <map>
extern std::map<char *, size_t> dvar;
void register_dvar(char *name, size_t sz);
#endif

    
#endif
