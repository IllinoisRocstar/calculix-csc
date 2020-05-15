#include <gtest.h>
#include <iostream>
#include <fstream>
#include <memory>
#include "clcxInterface.H"


std::string arg_fName1;
std::string arg_fName2;

// Aux functions
int compareFiles(std::string fn1, std::string fn2)
{
    std::fstream f1,f2;
    int pass = 1;
	f1.open(fn1, std::fstream::in);
	f2.open(fn2, std::fstream::in);
	char string1[256], string2[256];
	int j = 0;
	while(!f1.fail())
	{
		f1.getline(string1,256);
		f2.getline(string2,256);
		j++;
		if(strcmp(string1, string2) != 0)
		{
            std::cout << j << "-line strings are not equal" << "\n";
			std::cout << "   " << string1 << "\n";
			std::cout << "   " << string2 << "\n";
            pass = 0;
            break;
		}
	}
	return pass;
}

int testBeam(const std::string &testFN, const std::string &refFN) {
  auto _ci = std::make_shared<clcx_interface>();
  std::cout << "Test input file : " << testFN << std::endl;
  _ci->set_jobname(testFN);
  _ci->init();
  _ci->run();
  return compareFiles(testFN+".dat", refFN);
}

TEST(CSC, Beam) {
  EXPECT_EQ(1, testBeam(arg_fName1, arg_fName2));
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 3);
  arg_fName1 = argv[1];
  arg_fName2 = argv[2];

  // run tests
  int res = RUN_ALL_TESTS();

  return res;
}
