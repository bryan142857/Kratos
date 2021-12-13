// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Vahid Galavi
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_retention/brooks_corey_law.h"

namespace Kratos
{
//-------------------------------------------------------------------------------------------------
BrooksCoreyLaw::BrooksCoreyLaw()
    : RetentionLaw()
{
}

//-------------------------------------------------------------------------------------------------
BrooksCoreyLaw::BrooksCoreyLaw(const BrooksCoreyLaw& rOther)
    : RetentionLaw(rOther)
{
}

//-------------------------------------------------------------------------------------------------
RetentionLaw::Pointer BrooksCoreyLaw::Clone() const
{
    return Kratos::make_shared<BrooksCoreyLaw>(*this);
}

//-------------------------------------------------------------------------------------------------
BrooksCoreyLaw::~BrooksCoreyLaw()
{
}

//-------------------------------------------------------------------------------------------------
double BrooksCoreyLaw::
    CalculateSaturation(Parameters &rParameters)
{
    KRATOS_TRY;

    const double &p = rParameters.GetFluidPressure();
    const Properties &rMaterialProperties = rParameters.GetMaterialProperties();

    if (p > 0.0)
    {
        const double &satMax = rMaterialProperties[SATURATED_SATURATION];
        const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
        const double &pa     = rMaterialProperties[BROOKS_COREY_AIR_ENTRY_PRESSURE];
        const double &Landa     = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX];
        

        double sat = satMin + (satMax - satMin) *  pow(pa/p, Landa);
        return sat;
    }
    else
    {
        return rMaterialProperties[SATURATED_SATURATION];
    }

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double BrooksCoreyLaw::
    CalculateEffectiveSaturation(Parameters &rParameters)
{
    KRATOS_TRY;

    const double sat = CalculateSaturation(rParameters);

    const auto &rMaterialProperties = rParameters.GetMaterialProperties();
    const double &satMax = rMaterialProperties[SATURATED_SATURATION];
    const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];

    double effectiveSat = (sat - satMin) / (satMax - satMin);

    return effectiveSat;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double BrooksCoreyLaw::
    CalculateDerivativeOfSaturation(Parameters &rParameters)
{
    KRATOS_TRY;
    const double &p = rParameters.GetFluidPressure();

    if (p > 0.0)
    {
        const auto &rMaterialProperties = rParameters.GetMaterialProperties();
        const double &satMax = rMaterialProperties[SATURATED_SATURATION];
        const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
        const double &pa     = rMaterialProperties[BROOKS_COREY_AIR_ENTRY_PRESSURE];
        const double &Landa  = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX];
        

        double dSdp = (satMax - satMin) * (-Landa) * pow(pa,Landa)*pow(p, -Landa-1.0);
    
        return dSdp;
    }
    else
    {
        return 0.0;
    }

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double BrooksCoreyLaw::
    CalculateRelativePermeability(Parameters &rParameters)
{
    KRATOS_TRY;

    const double effSat = CalculateEffectiveSaturation(rParameters);

    const auto &rMaterialProperties = rParameters.GetMaterialProperties();
    const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
    const double &Landa  = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX];

    double relPerm = pow(effSat, ((2+3*Landa)/Landa)); 

    const double &minRelPerm = rMaterialProperties[MINIMUM_RELATIVE_PERMEABILITY];

    relPerm = std::max(relPerm, minRelPerm);

    return relPerm;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double BrooksCoreyLaw::
    CalculateBishopCoefficient(Parameters &rParameters)
{
    KRATOS_TRY;
     const double &p = rParameters.GetFluidPressure();
     
if (p > 0.0)
    {
        const auto &rMaterialProperties = rParameters.GetMaterialProperties();
        const double &Porosity = rMaterialProperties[POROSITY];
        const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
        const double &pa     = rMaterialProperties[BROOKS_COREY_AIR_ENTRY_PRESSURE];
        const double &Landa  = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX];
        const double &Beta  = rMaterialProperties[BROOKS_COREY_FITTING_PARAMETER];// which is considered between 0.4 to 0.7

        double BishopCo = pow(pa/p, Beta)+pow(pa/p, 1+Beta)*Porosity*(Landa/(Landa-1))*(1-satMin)*(1-pow(pa/p, Landa-1)); 
   
        return BishopCo;        

        
    }
    else
    {
        return 1.0;
    }


    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double& BrooksCoreyLaw::CalculateValue(RetentionLaw::Parameters& rParameterValues,
                                        const Variable<double>& rThisVariable,
                                        double& rValue)
{
    if (rThisVariable == DEGREE_OF_SATURATION)
    {
        rValue = this->CalculateSaturation(rParameterValues);
        return rValue;
    }
    else if (rThisVariable == EFFECTIVE_SATURATION)
    {
        rValue = this->CalculateEffectiveSaturation(rParameterValues);
        return rValue;
    }
    else if (rThisVariable == BISHOP_COEFICIENT)
    {
        rValue = this->CalculateBishopCoefficient(rParameterValues);
        return rValue;
    }
    else if (rThisVariable == DERIVATIVE_OF_SATURATION)
    {
        rValue = this->CalculateDerivativeOfSaturation(rParameterValues);
        return rValue;
    }
    else if (rThisVariable == RELATIVE_PERMEABILITY)
    {
        rValue = this->CalculateRelativePermeability(rParameterValues);
        return rValue;
    }

    return( rValue );
}

//------------------------- RETENSION LAW GENERAL FEATURES ----------------------------------------
//-------------------------------------------------------------------------------------------------
void BrooksCoreyLaw::
    InitializeMaterial(const Properties& rMaterialProperties,
                       const GeometryType& rElementGeometry,
                       const Vector& rShapeFunctionsValues)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void BrooksCoreyLaw::
    Initialize(Parameters &rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void BrooksCoreyLaw::
    InitializeSolutionStep(Parameters &rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void BrooksCoreyLaw::
    Finalize(Parameters &rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void BrooksCoreyLaw::
    FinalizeSolutionStep(Parameters &rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
int BrooksCoreyLaw::Check(const Properties& rMaterialProperties,
                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(!rMaterialProperties.Has(SATURATED_SATURATION))
                    << "SATURATED_SATURATION is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] < 0.0)
                    << "SATURATED_SATURATION cannot be less than 0 " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(RESIDUAL_SATURATION))
                    << "RESIDUAL_SATURATION is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[RESIDUAL_SATURATION] < 0.0)
                    << "RESIDUAL_SATURATION cannot be less than 0 " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(SATURATED_SATURATION))
                    << "SATURATED_SATURATION is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] > 1.0)
                    << "SATURATED_SATURATION cannot be greater than 1.0 " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(RESIDUAL_SATURATION))
                    << "RESIDUAL_SATURATION is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[RESIDUAL_SATURATION] > 1.0)
                    << "RESIDUAL_SATURATION cannot be greater than 1.0 " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(SATURATED_SATURATION))
                    << "SATURATED_SATURATION is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] < rMaterialProperties[RESIDUAL_SATURATION])
                    << "RESIDUAL_SATURATION cannot be greater than SATURATED_SATURATION " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(BROOKS_COREY_AIR_ENTRY_PRESSURE))
                    << "BROOKS_COREY_AIR_ENTRY_PRESSURE is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(!(rMaterialProperties[BROOKS_COREY_AIR_ENTRY_PRESSURE] > 0.0))
                    << "BROOKS_COREY_AIR_ENTRY_PRESSURE must be greater than 0 " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(MINIMUM_RELATIVE_PERMEABILITY))
                    << "MINIMUM_RELATIVE_PERMEABILITY is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(!(rMaterialProperties[MINIMUM_RELATIVE_PERMEABILITY] > 0.0))
                    << "MINIMUM_RELATIVE_PERMEABILITY must be greater than 0 " << std::endl;

    return 0;
}

} // Namespace Kratos