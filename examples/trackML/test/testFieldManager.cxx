#include <Particle.h>
#include <Track.h>
#include <FieldManager.h>

int main(int argc, char const *argv[])
{
  using ThreeVector       = vecgeom::Vector3D<double>;
  using Track = vectorflow::Track;
  float parFieldVector[3] = {0, 0, 2}; // Tesla
  
  trackml::FieldManager fm;
  fm.UseConstantMagField(parFieldVector, "tesla");
  fm.SetEpsilonRK(0.0003);
  fm.CreateFieldAndSolver();

  geantphysics::Particle::DefineParticles();
  const char *pname = "e-";
  const geantphysics::Particle *particle = geantphysics::Particle::GetParticleByName(std::string(pname));
  if (!particle) {
    printf("AddPrimary: Error: particle %s was not defined\n", pname);
    return 1;
  }
  return 0;
}
