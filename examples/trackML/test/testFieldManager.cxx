#include <FieldManager.h>

int main(int argc, char const *argv[])
{
  /* code */
  float parFieldVector[3] = {0, 0, 2}; // Tesla
  
  trackml::FieldManager fm;
  fm.UseConstantMagField(parFieldVector, "tesla");
  fm.SetEpsilonRK(0.0003);
  fm.CreateFieldAndSolver();
  return 0;
}