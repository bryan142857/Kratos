// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Sergio Jimenez
//                   Alejandro Cornejo
//  Collaborator:    Lucia Barbu
//

#if !defined(KRATOS_UNIFIED_FATIGUE_RULE_OF_MIXTURES_H_INCLUDED)
#define KRATOS_UNIFIED_FATIGUE_RULE_OF_MIXTURES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "custom_constitutive/generic_small_strain_isotropic_damage.h"
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/linear_plane_strain.h"


namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
/**
 * @class SerialParallelRuleOfMixturesLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This CL implements the serial-parallel rule of mixtures developed by F.Rastellini
 * @details
 * @author Alejandro Cornejo
 */
template <class TConstLawIntegratorType>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) UnifiedFatigueRuleOfMixturesLaw
    : public std::conditional<TConstLawIntegratorType::VoigtSize == 6, ElasticIsotropic3D, LinearPlaneStrain >::type
{
  public:
    ///@name Type Definitions
    ///@{

    /// The node definition
    typedef Node<3> NodeType;

    /// The geometry definition
    typedef Geometry<NodeType> GeometryType;


    /// Definition of the machine precision tolerance
    static constexpr double machine_tolerance = std::numeric_limits<double>::epsilon();
    static constexpr double threshold_tolerance = 1.0e-4;

    /// The define the working dimension size, already defined in the integrator
    static constexpr SizeType Dimension = TConstLawIntegratorType::Dimension;

    /// The define the Voigt size, already defined in the  integrator
    static constexpr SizeType VoigtSize = TConstLawIntegratorType::VoigtSize;

    typedef typename std::conditional<VoigtSize == 6, ElasticIsotropic3D, LinearPlaneStrain >::type BaseType;
    /// Counted pointer of SerialParallelRuleOfMixturesLaw
    KRATOS_CLASS_POINTER_DEFINITION(UnifiedFatigueRuleOfMixturesLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    UnifiedFatigueRuleOfMixturesLaw()
    {
    }

    /**
    * Constructor.
    */
    UnifiedFatigueRuleOfMixturesLaw(double HCFVolParticipation)
        : mHCFVolumetricParticipation(HCFVolParticipation)
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<UnifiedFatigueRuleOfMixturesLaw<TConstLawIntegratorType>>(*this);
    }

    // Copy constructor
    UnifiedFatigueRuleOfMixturesLaw(UnifiedFatigueRuleOfMixturesLaw const& rOther)
        :  BaseType(rOther),
            mpHCFConstitutiveLaw(rOther.mpHCFConstitutiveLaw),
            mpULCFConstitutiveLaw(rOther.mpULCFConstitutiveLaw),
            mHCFVolumetricParticipation(rOther.mHCFVolumetricParticipation)
    {
    }

    /**
    * Destructor.
    */
    ~UnifiedFatigueRuleOfMixturesLaw() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new constitutive law pointer
     * @param NewParameters The configuration parameters of the new constitutive law
     * @return a Pointer to the new constitutive law
     */
    ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override;

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 3;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() const override
    {
        return 6;
    };

    /**
     * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * to be called at the end of each solution step
     * (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void FinalizeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<bool>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (integer)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<int>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Matrix)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Matrix>& rThisVariable) override;

    /**
     * @brief Sets the value of a specified variable (bool)
     * @param rThisVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<bool>& rThisVariable,
        const bool& Value,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Sets the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<int>& rThisVariable,
        const int& rValue,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<double> &rThisVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (Vector)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Vector> &rThisVariable,
        const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) override;

    /**
     * @brief Returns the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    int& GetValue(const Variable<int>& rThisVariable, int& rValue) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;

    /**
     * @brief Returns the value of a specified variable (Matrix)
     * @param rThisVariable the variable to be returned
     * @return rValue output: the value of the specified variable
     */
    Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue) override;

    /**
     * @brief Calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue) override;

    /**
     * @brief Calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue) override;

    /**
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue) override;

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues) override;

    /**
     * @brief Initialize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * This method computes the stresses of the matrix/fiber according to its own CL
     * @param rValues the needed parameters for the CL calculation
     * @param rHCFStrainVector the strain vector of the HCF model
     * @param rULCFStrainVector the strain vector of the ULCF model
     * @param rHCFStressVector  the stress vector of the matrix
     * @param rULCFStressVector  the stress vector of the fiber
     */
    void IntegrateStressesOfHCFAndULCFModels(
        ConstitutiveLaw::Parameters& rValues,
        Vector rHCFStrainVector,
        Vector rULCFStrainVector,
        Vector& rHCFStressVector,
        Vector& rULCFStressVector);

    /**
     * This method computes the stresses of the matrix/fiber according to its own CL
     * @param values_HCF the needed parameters for the CL calculation
     */
    void CalculateMaterialResponseHCFModel(
        ConstitutiveLaw::Parameters& values_HCF);

    /**
     * This method computes the stresses of the matrix/fiber according to its own CL
     * @param values_ULCF the needed parameters for the CL calculation
     */
    void CalculateMaterialResponseULCFModel(
        ConstitutiveLaw::Parameters& values_ULCF);

    /**
     * This method computes the stresses of the matrix/fiber according to its own CL
     * @param values_HCF the needed parameters for the CL calculation
     */
    void FinalizeMaterialResponseHCFModel(
        double& rHighCycleFatiguePredictiveUniaxialStress,
        Vector& rHighCycleFatigueStressVector,
        ConstitutiveLaw::Parameters& values_HCF);

    /**
     * This method computes the stresses of the matrix/fiber according to its own CL
     * @param values_ULCF the needed parameters for the CL calculation
     */
    void FinalizeMaterialResponseULCFModel(
        double& rUltraLowCycleFatiguePredictiveUniaxialStress,
        Vector& rUltraLowCycleFatigueStressVector,
        ConstitutiveLaw::Parameters& values_ULCF);

    /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues);

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return true;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return true;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

  protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    double& GetHCFThreshold() { return mHCFThreshold; }
    double& GetULCFThreshold() { return mULCFThreshold; }
    double& GetDamage() { return mDamage; }
    double& GetPlasticDissipation() { return mPlasticDissipation; }
    Vector& GetPlasticStrain() { return mPlasticStrain; }

    void SetHCFThreshold(const double toHCFThreshold) { mHCFThreshold = toHCFThreshold; }
    void SetULCFThreshold(const double toULCFThreshold) { mULCFThreshold = toULCFThreshold; }
    void SetDamage(const double toDamage) { mDamage = toDamage; }
    void SetPlasticDissipation(const double toPlasticDissipation) { mPlasticDissipation = toPlasticDissipation; }
    void SetPlasticStrain(const array_1d<double, VoigtSize>& rPlasticStrain) { mPlasticStrain = rPlasticStrain; }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ConstitutiveLaw::Pointer mpHCFConstitutiveLaw;
    ConstitutiveLaw::Pointer mpULCFConstitutiveLaw;
    double mHCFVolumetricParticipation;
    double mFatigueReductionFactor = 1.0;
    Vector mPreviousStresses = ZeroVector(2); // [S_t-2, S_t-1]
    double mMaxStress = 0.0;
    double mMinStress = 0.0;
    bool mMaxDetected = false; // Maximum's indicator in the current cycle
    bool mMinDetected = false; // Minimum's indicator in the current cycle
    unsigned int mNumberOfCyclesGlobal = 1; // Total number of cycles in the whole analysis
    unsigned int mNumberOfCyclesLocal = 1; // Equivalent number of cycles for the current cyclic load
    bool mNewCycleIndicator = false; // New cycle identifier required for the advancing process.
    double mPreviousMaxStress = 0.0;
    double mPreviousMinStress = 0.0;
    double mFatigueReductionParameter = 0.0; // B0
    double mFatigueLimit = 0.0; // Endurance limit of the fatigue model.
    double mReversionFactorRelativeError = 0.0; // Relative error of the R = Smin / Smax between cycles inducing recalculation of Nlocal and advanciing process.
    double mMaxStressRelativeError = 0.0; // Relative error of Smax between cycles inducing recalculation of Nlocal and advanciing process.
    double mCyclesToFailure = 0.0; // Nf. Required for the advanciing process.
    double mWohlerStress = 1.0; // Normalised Wohler stress required for building the life prediction curves (SN curves)
    double mHCFThreshold = 0.0;
    double mULCFThreshold = 0.0;
    double mPlasticDissipation = 0.0;
    double mDamage = 0.0;
    Vector mPlasticStrain = ZeroVector(VoigtSize);
    double mPreviousCycleTime = 0.0; // Instanced variable used in the advanciing process for the conversion between time and number of cycles.
    double mPeriod = 0.0; // Instanced variable used in the advanciing process for the conversion between time and number of cycles.

    //Variable used while updating the volumetric participation.
    double mReferenceVolumetricParticipation = 0.0; //Reference volumetric participation when a new load block is detected
    double mReferenceDamage = 0.0; //Reference level when a new load block is detected
    double mReferencePlasticDissipation = 0.0; //Reference equivalent plastic dissipation when a new load block is detected
    double mReferenceFatigueReductionFactor = 1.0; //Reference equivalent plastic dissipation when a new load block is detected



    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("HCFConstitutiveLaw", mpHCFConstitutiveLaw);
        rSerializer.save("ULCFConstitutiveLaw", mpULCFConstitutiveLaw);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("HCFConstitutiveLaw", mpHCFConstitutiveLaw);
        rSerializer.load("ULCFConstitutiveLaw", mpULCFConstitutiveLaw);
    }

    ///@}

}; // Class SerialParallelRuleOfMixturesLaw

} // namespace Kratos

#endif