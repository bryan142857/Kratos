#if !defined(DEM_D_LINEAR_HIGHSTIFFNESS_CL_H_INCLUDED)
#define DEM_D_LINEAR_HIGHSTIFFNESS_CL_H_INCLUDED

#include <string>
#include <iostream>
#include "DEM_D_Linear_viscous_Coulomb_CL.h"

namespace Kratos {

    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_D_Linear_HighStiffness : public DEM_D_Linear_viscous_Coulomb {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Linear_HighStiffness);

        DEM_D_Linear_HighStiffness() {}

        ~DEM_D_Linear_HighStiffness() {}

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        std::unique_ptr<DEMDiscontinuumConstitutiveLaw> CloneUnique() override;

        void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) override;

        void InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta = 0.0) override;

    };

} /* namespace Kratos.*/
#endif /* DEM_D_LINEAR_HIGHSTIFFNESS_CL_H_INCLUDED  defined */
