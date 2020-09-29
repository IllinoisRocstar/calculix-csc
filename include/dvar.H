#ifndef DVAR_H
#define DVAR_H

#ifdef __cplusplus
#include <map>
#include <string>
extern std::map<std::string, size_t> dvar;
void register_dvar(std::string, size_t sz);
#endif

#endif
