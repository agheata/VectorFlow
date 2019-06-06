#ifndef TRACKML_FIELD_MANAGER
#define TRACKML_FIELD_MANAGER

#include <vectorFlow/SystemOfUnits.h>
#include <vectorFlow/Typedefs.h>
#include <base/Vector3D.h>

#include <Geant/MagFieldEquation.h>
#include <Geant/GUFieldPropagator.h>
#include <Geant/ScalarMagFieldEquation.h>
#include <Geant/UniformMagField.h>
#include <Geant/StepperFactory.h>
#include <Geant/CashKarp.h>
#include <Geant/ScalarIntegrationDriver.h>
#include <Geant/SimpleIntegrationDriver.h>

namespace trackml {
/// Manager for the 2 Tesla uniform magnetic field from TrackML
class FieldManager {

public:
  void SetUseRK(bool flag) { fUseRK = flag; }
  void SetEpsilonRK(double val) { fEpsilonRK = val; }
  void SetMinimumStep(double v) { fMinimumStepInField = v; }

  // Inquiry methods
  bool IsFieldUniform() { return fUseUniformField; }
  double GetEpsilonRK() { return fEpsilonRK; }
  double GetMinimumStep() { return fMinimumStepInField; }
  bool IsFieldCreated() { return fCreatedField; }
  VVectorField *GetField() { return fpField; }
  static void GetFieldValue(const vecgeom::Vector3D<double> &Position,
                            vecgeom::Vector3D<double> &MagFieldValue,
                            double &bmag);
private:
  double fEpsilonRK = 3.0e-5;
  double fMinimumStepInField = 1.0e-4; // GV units = cm;
  vecgeom::Vector3D<float> fMagFieldValue;
  bool fUseUniformField = true;
  bool fZeroField = true;
  bool fUseRK = true;

  bool fCreatedField;
  UniformMagField *fUniformField = nullptr;
  VVectorField *fpField = nullptr;
  GUFieldPropagator *fPropagator = nullptr; ///< This is the field propagator

public:
  FieldManager() {}
  ~FieldManager(){};

  template <class Field_t>
  void CreateSolverForField(Field_t *ptrField)
  {
    fPropagator = CreatePropagator<Field_t>(*ptrField, fEpsilonRK, fMinimumStepInField);
    fCreatedField = true;
  }

  void UseConstantMagField(float fieldVal[3], const char *Units = 0)
  {
    const char *methodName = "UserFieldConstruction::UseConstantMagField";
    bool defaultUsed       = false;
    double unit            = 1;

    if (Units == 0 || strcmp(Units, "kilogauss") == 0) {
      unit        = geant::units::kilogauss;
      defaultUsed = (Units == 0);
    } else if ((strcmp(Units, "gauss") == 0) || (strcmp(Units, "Gauss") == 0)) {
      unit = geant::units::gauss;
    } else if ((strcmp(Units, "tesla") == 0) || (strcmp(Units, "Tesla") == 0)) {
      unit = geant::units::tesla;
    } else {
      unit        = geant::units::kilogauss;
      defaultUsed = (Units == 0);
    }
    if (defaultUsed) printf("%s - WARNING: No units provided - using kilogauss as default unit", methodName);
    fMagFieldValue = vecgeom::Vector3D<float>(fieldVal[0] * unit, fieldVal[1] * unit, fieldVal[2] * unit);

    fUseUniformField = true;
    fZeroField       = (fMagFieldValue.Mag2() == 0.0);
  }

  /** @brief Create the global magnetic field and classes to integrate it. Register field. */
  /** @description  Must call the templated CreateSolverForField method.                   */
  void CreateFieldAndSolver(bool uniform = true)
  {
    printf("CreateFieldAndSolver:  Use uniform= %d  Value= %f %f %f - kilogauss.  Zero-Flag= %d",
           fUseUniformField, fMagFieldValue[0] / geant::units::kilogauss,
           fMagFieldValue[1] / geant::units::kilogauss, fMagFieldValue[2] / geant::units::kilogauss, fZeroField);

    if (uniform) {
      fUseUniformField = true;
      fUniformField = new UniformMagField(fMagFieldValue);
      bool isUniform      = true;
      fpField = fUniformField;

      CreateSolverForField<UniformMagField>(fUniformField);

      if (fZeroField) {
        printf(" Zero Magnetic Field configured.");
      }
    } else {
      // magnetic field map here
      printf("CreateFieldAndSolver: Field map not supported yet\n");
    }
  }

  template <typename Field_t>
  inline GUFieldPropagator *CreatePropagator(Field_t &gvField, double relTolerance, double minStep)
  {
    // Create scalar driver
    constexpr unsigned int kNposmom = 6; // Position 3-vec + Momentum 3-vec
    using Equation_t = ScalarMagFieldEquation<Field_t, kNposmom>;
    auto gvEquation  = new Equation_t(&gvField);
    auto aStepper = StepperFactory::CreateStepper<Equation_t>(gvEquation);
    auto scalarDriver = new ScalarIntegrationDriver(minStep, aStepper, kNposmom, false);

    // Create flexible driver
    using FlexEquation_t = MagFieldEquation<Field_t>;
    auto gvFlexEquation  = new FlexEquation_t(&gvField);
    using StepperType = CashKarp<FlexEquation_t, kNposmom>;
    auto myStepper    = new StepperType(gvFlexEquation);
    using DriverType  = SimpleIntegrationDriver<StepperType, kNposmom>;
    auto flexibleDriver = new DriverType(minStep, myStepper, kNposmom, false);

    GUFieldPropagator *fieldPropagator = new GUFieldPropagator(scalarDriver, relTolerance, flexibleDriver);
    return fieldPropagator;
  }

};
} // namespace trackml
#endif
