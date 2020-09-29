#include "dvar.H"
#include <map>
#include <string>

// extern std::map<char *, size_t> dvar;

// void register_dvar(char *name, size_t sz) {
void register_dvar(std::string name, size_t sz) {
  /*std::cout << "Register " << std::string(name) << " size " << sz <<
   * std::endl;*/
  dvar[name] = sz;
}
