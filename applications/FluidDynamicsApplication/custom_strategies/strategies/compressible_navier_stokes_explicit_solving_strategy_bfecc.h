//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla, Edurd Gómez
//
//

#if !defined(KRATOS_EXPLICIT_SOLVING_STRATEGY_BFECC)
#define KRATOS_EXPLICIT_SOLVING_STRATEGY_BFECC

/* System includes */
#include <numeric>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "factories/factory.h"
#include "explicit_solving_strategy_bfecc.h"

namespace Kratos
{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class CompressibleNavierStokesExplicitSolvingStrategyBFECC : public ExplicitSolvingStrategyBFECC<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    // The base solving strategy class definition
    typedef SolvingStrategy<TSparseSpace, TDenseSpace> SolvingStrategyType;

    // The base class definition
    typedef ExplicitSolvingStrategyBFECC<TSparseSpace, TDenseSpace> BaseType;

    /// The definition of the current class
    typedef CompressibleNavierStokesExplicitSolvingStrategyBFECC<TSparseSpace, TDenseSpace, TButcherTableau> ClassType;

    // The explicit builder and solver definition
    typedef typename BaseType::ExplicitBuilderType ExplicitBuilderType;

    /// The DOF type
    typedef typename BaseType::DofType DofType;

    /// The local vector definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;
    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(CompressibleNavierStokesExplicitSolvingStrategyBFECC);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (empty)
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyBFECC()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyBFECC(
        ModelPart &rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param pExplicitBuilder The pointer to the explicit builder and solver
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyBFECC(
        ModelPart &rModelPart,
        typename ExplicitBuilderType::Pointer pExplicitBuilder,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, pExplicitBuilder, MoveMeshFlag, RebuildLevel)
    {
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyBFECC(
        ModelPart &rModelPart,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, MoveMeshFlag, RebuildLevel)
    {
    }

    /**
     * @brief Create method
     * @param rModelPart The model part to be computed
     * @param ThisParameters The configuration parameters
     */
    typename SolvingStrategyType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    /** Copy constructor.
     */
    CompressibleNavierStokesExplicitSolvingStrategyBFECC(const CompressibleNavierStokesExplicitSolvingStrategyBFECC &Other) = delete;

    /** Destructor.
     */
    ~CompressibleNavierStokesExplicitSolvingStrategyBFECC() override = default;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "compressible_navier_stokes_explicit_solving_strategy_bfecc"
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        std::stringstream s;
        s << "compressible_navier_stokes_explicit_solving_strategy_bfecc"
        return s.str();
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    int Check() override
    {
        int err_code = BaseType::Check();

        // Check that the process info already contains the OSS activation variable
        const auto& r_process_info = BaseType::GetModelPart().GetProcessInfo();
        KRATOS_ERROR_IF_NOT(r_process_info.Has(OSS_SWITCH)) << "OSS_SWITCH variable has not been set in the ProcessInfo container. Please set it in the strategy \'Initialize\'." << std::endl;

        // Shock capturing process check
        if (mpShockCapturingProcess) {
            mpShockCapturingProcess->Check();
        }

        return err_code;
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        KRATOS_TRY

        Parameters default_parameters = Parameters(R"(
        {
            "name" : "compressible_navier_stokes_explicit_explicit_solving_strategy_bfecc",
            "rebuild_level" : 0,
            "move_mesh_flag": false,
            "calculate_non_conservative_magnitudes" : true,
            "shock_capturing_settings" : { }
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "compressible_navier_stokes_explicit_solving_strategy_runge_kutta_4";
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        // Base class assign settings call
        BaseType::AssignSettings(ThisParameters);
        mCalculateNonConservativeMagnitudes = ThisParameters["calculate_non_conservative_magnitudes"].GetBool();

        SetUpShockCapturing(ThisParameters["shock_capturing_settings"]);

        if (mpShockCapturingProcess && !mCalculateNonConservativeMagnitudes) {
            KRATOS_WARNING("CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4")
                << "\'shock_capturing\' requires \'calculate_non_conservative_magnitudes\' to be active. Activating it." << std::endl;
            mCalculateNonConservativeMagnitudes = true;
        }
    }

    void SetUpShockCapturing(Parameters ShockCapturingParameters)
    {
        KRATOS_TRY

        auto defaults = Parameters(R"(
            {
                "type" : "physics_based",
                "Parameters" : {
                    "model_part_name" : ""
                }
            }
        )");
        defaults["Parameters"]["model_part_name"].SetString(BaseType::GetModelPart().Name());

        ShockCapturingParameters.ValidateAndAssignDefaults(defaults);
        ShockCapturingParameters.RecursivelyAddMissingParameters(defaults);

        using ShockCapturingFactoryType = Process::UniquePointer (*)(ModelPart&, Parameters);

        const std::map<const std::string, ShockCapturingFactoryType> shock_capturing_factory_map
        {
            {"none"         , [](ModelPart& m, Parameters p) -> Process::UniquePointer {return nullptr;}},
            {"physics_based", [](ModelPart& m, Parameters p) -> Process::UniquePointer {return Kratos::make_unique<ShockCapturingPhysicsBasedProcess>(m, p);}},
            {"entropy_based", [](ModelPart& m, Parameters p) -> Process::UniquePointer {return Kratos::make_unique<ShockCapturingEntropyViscosityProcess>(m, p);}}
        };

        const auto sc_type = ShockCapturingParameters["type"].GetString();

        const auto it = shock_capturing_factory_map.find(sc_type);
        if(it != shock_capturing_factory_map.end())
        {
            mpShockCapturingProcess = it->second(BaseType::GetModelPart(), ShockCapturingParameters["Parameters"]);
        }
        else
        {
            std::stringstream msg;
            msg << "Provided shock capturing type \""<< sc_type <<"\" does not match any of the available ones.\n";
            msg << "Please choose one from the following list:\n";
            for(const auto& keyvaluepair: shock_capturing_factory_map)
            {
                msg <<" - " << keyvaluepair.first << "\n";
            }
            msg << std::endl;

            KRATOS_ERROR << msg.str();
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Initialization of member variables and prior operations
     * In this method we call the base strategy initialize and initialize the time derivatives
     * This is required to prevent OpenMP errors as the time derivatives are stored in the non-historical database
     */
    void Initialize() override
    {
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();

        // Call the base finalize substep method
        BaseType::Initialize();

        // Initialize the postprocess non-historical variables
        block_for_each(r_model_part.Nodes(), [](Node<3>& rNode){
            rNode.SetValue(MACH, 0.0);
            rNode.SetValue(SOUND_VELOCITY, 0.0);
        });

        // If required, initialize the OSS projection variables
        if (r_process_info[OSS_SWITCH]) {
            block_for_each(r_model_part.Nodes(), [](Node<3>& rNode){
                rNode.SetValue(DENSITY_PROJECTION, 0.0);
                rNode.SetValue(TOTAL_ENERGY_PROJECTION, 0.0);
                rNode.SetValue(MOMENTUM_PROJECTION, ZeroVector(3));
            });
        }

        // If required, initialize the physics-based shock capturing variables
        if (mpShockCapturingProcess) {
            mpShockCapturingProcess->ExecuteInitialize();
        }
    }

    void InitializeSolutionStep() override
    {
        BaseType::InitializeSolutionStep();

        if (mCalculateNonConservativeMagnitudes) {
            CalculateNonConservativeMagnitudes();
        }

        if (mpShockCapturingProcess) {
            mpShockCapturingProcess->ExecuteInitializeSolutionStep();
        }
    }

    /**
     * @brief Finalize the Runge-Kutta step
     * In this method we calculate the final linearised time derivatives after the final update
     * These will be the time derivatives employed in the first RK4 sub step of the next time step
     */
    void FinalizeSolutionStep() override
    {
        // Call the base RK4 finalize substep method
        BaseType::FinalizeSolutionStep();

        // Apply the momentum slip condition
        if (mApplySlipCondition) {
            ApplySlipCondition();
        }

        // Postprocess the non-conservative magnitudes
        // This needs to be done before the shock capturing as it is based on these
        if (mCalculateNonConservativeMagnitudes) {
            CalculateNonConservativeMagnitudes();
        }

        // Perform the shock capturing detection and artificial values calculation
        // This needs to be done at the end of the step in order to include the future shock
        // capturing magnitudes in the next automatic dt calculation
        if (mpShockCapturingProcess) {
            mpShockCapturingProcess->ExecuteFinalizeSolutionStep();
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream ss;
        ss << "CompressibleNavierStokesExplicitSolvingStrategyBFECC with tableau " << mButcherTableau.Name();
        return ss.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

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

    /**
     * @brief Initialize the BFECC intermediate substep
     * This method is intended to implement all the operations required before each BFECC intermediate substep
     */
    virtual void InitializeBFECCIntermediateSubStep(SizeType SubstepIndex) override
    {
        if(SubstepIndex == 1)
        {
            StashDiffusiveConstants();
        }
    };

    /**
     * @brief Initialize the BFECC last substep
     * This method is intended to implement all the operations required before each BFECC last substep
     */
    virtual void InitializeBFECCLastSubStep(SizeType SubstepIndex) override
    {
        PopDiffusiveConstants();
    };


    void SubstepForward(const LocalSystemVectorType& rFixedDofsValues)
    {
        return PerformSubStep(1, 0.0, rFixedDofsValues, 1.0, false);
    }

    void SubstepBackward(const LocalSystemVectorType& rFixedDofsValues)
    {
        return PerformSubStep(2, 1.0, rFixedDofsValues,-1.0, false);
    }

    void FinalSubstep(const LocalSystemVectorType& rFixedDofsValues)
    {
        return PerformSubStep(3, 0.0, rFixedDofsValues, 1.0, true);
    }

    /**
     * @brief Performs an intermediate RK4 step
     * This functions performs all the operations required in an intermediate RK4 sub step
     * @param SubStepIndex The sub step index
     * @param TimeIntegrationTheta The point in time-step to evaulate in
     * @param rFixedDofsValues The vector containing the step n+1 values of the fixed DOFs
     * @param TimeDirection Whether going forawrds or backwards
     * @param LastSubstep Whether this is the last substep or not
     */
    virtual void PerformSubStep(
        const IndexType SubstepIndex,
        const double TimeIntegrationTheta,
        const LocalSystemVectorType& rFixedDofsValues,
        const double TimeDirection,
        bool LastSubstep)
    {
        KRATOS_TRY

        // Get the required data from the explicit builder and solver
        auto& explicit_bs = BaseType::GetExplicitBuilder();
        auto& r_dof_set = explicit_bs.GetDofSet();
        const auto& r_lumped_mass_vector = explicit_bs.GetLumpedMassMatrixVector();

        // Get model part and information
        const double dt = BaseType::GetDeltaTime();
        KRATOS_ERROR_IF(dt < 1.0e-12) << "ProcessInfo DELTA_TIME is close to zero." << std::endl;
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();

        // Perform the intermidate sub step update
        r_process_info.GetValue(TIME_INTEGRATION_THETA) = TimeIntegrationTheta;

        if(LastSubstep) {
            InitializeBFECCLastSubStep(SubstepIndex);
        } else {
            InitializeBFECCIntermediateSubStep(SubstepIndex);
        }

        explicit_bs.BuildRHS(r_model_part);

        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;

                // Save current value in the corresponding vector
                const double residual = it_dof->GetSolutionStepReactionValue();

                // Do the DOF update
                double& r_u = it_dof->GetSolutionStepValue(0);
                double& r_u_prev_step = it_dof->GetSolutionStepValue(0);

                if (!it_dof->IsFixed()) {
                    const double mass = r_lumped_mass_vector[i_dof];
                    r_u += TimeDirection * (dt / mass) * residual;
                } else {
                    r_u = integration_theta*rFixedDofsValues[i_dof] +(1-integration_theta)*r_u_prev_step;
                }
            }
        );

        if(LastSubstep) {
            FinalizeBFECCLastSubStep(SubstepIndex);
        } else {
            FinalizeBFECCIntermediateSubStep(SubstepIndex);
        }

        KRATOS_CATCH("SubstepIndex = " + std::to_string(SubStepIndex));
    }

    void CorrectErrorAfterForwardsAndBackwards()
    {
        auto& explicit_bs = BaseType::GetExplicitBuilder();
        auto& r_dof_set = explicit_bs.GetDofSet();

        block_for_each(r_dof_set, [&](Dof<double>& r_dof)
        {
            const double old_value = r_dof.GetSolutionStepValue(1);
            double& solution_step_value = r_dof.GetSolutionStepValue(0);
            const double error = (solution_step_value - old_value)/2.0;
            solution_step_value = old_value - error;
        }
        );
    }

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

    bool mApplySlipCondition = true;
    bool mCalculateNonConservativeMagnitudes = true;

    Process::UniquePointer mpShockCapturingProcess = nullptr;


    struct Stash {
        double conductivity = 0.0;
        double dynamic_viscosity = 0.0;
    } mDiffusionStash;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void CalculateOrthogonalSubScalesProjection()
    {
        // Get model part data
        auto& r_model_part = BaseType::GetModelPart();
        const int n_nodes = r_model_part.NumberOfNodes();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        const unsigned int block_size = r_process_info[DOMAIN_SIZE] + 2;

        // Get the required data from the explicit builder and solver
        // The lumped mass vector will be used to get the NODAL_AREA for the residuals projection
        const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Initialize the projection values
        block_for_each(r_model_part.Nodes(), [](Node<3>& rNode){
            rNode.GetValue(DENSITY_PROJECTION) = 0.0;
            rNode.GetValue(MOMENTUM_PROJECTION) = ZeroVector(3);
            rNode.GetValue(TOTAL_ENERGY_PROJECTION) = 0.0;
        });

        // Calculate the residuals projection
        std::tuple<double, double, array_1d<double,3>> oss_proj_tls;
        block_for_each(r_model_part.Elements(), oss_proj_tls, [&](Element& rElement, std::tuple<double, double, array_1d<double,3>>& rOssProjTLS){
            double& rho_proj = std::get<0>(rOssProjTLS);
            double& tot_ener_proj = std::get<1>(rOssProjTLS);
            array_1d<double,3>& mom_proj = std::get<2>(rOssProjTLS);
            rElement.Calculate(DENSITY_PROJECTION, rho_proj, r_process_info);
            rElement.Calculate(MOMENTUM_PROJECTION, mom_proj, r_process_info);
            rElement.Calculate(TOTAL_ENERGY_PROJECTION, tot_ener_proj, r_process_info);
        });

        // Do thhe nodal weighting
        // Note that to avoid calculating the NODAL_AREA we took from the first DOF of the lumped mass vector
        IndexPartition<>(n_nodes).for_each([&](const ModelPart::SizeType iNode){
            auto it_node = r_model_part.NodesBegin() + iNode;
            const double nodal_area = r_lumped_mass_vector[iNode * block_size];
            it_node->GetValue(DENSITY_PROJECTION) /= nodal_area;
            it_node->GetValue(MOMENTUM_PROJECTION) /= nodal_area;
            it_node->GetValue(TOTAL_ENERGY_PROJECTION) /= nodal_area;
        });

    }

    /**
     * @brief Calculates the non-conservative magnitudes
     * This function calculates the non-conservative magnitudes from the obtained conservative ones
     * The set of non-conservative magnitudes (VELOCITY, PRESSURE and TEMPERATURE) are saved in the
     * historical database. The speed of sound (SOUND_VELOCITY) and Mach number (MACH) are also
     * computed and stored in the non-historical database.
     */
    void CalculateNonConservativeMagnitudes()
    {
        // Get fluid properties from first element in mesh
        // Note that in here we assume that the fluid is single-phase
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_prop = r_model_part.ElementsBegin()->GetProperties();
        const double c_v = r_prop.GetValue(SPECIFIC_HEAT); // Specific heat at constant volume
        const double gamma = r_prop.GetValue(HEAT_CAPACITY_RATIO); // Heat capacity ratio
        const double R = (gamma - 1.0) * c_v; // Ideal gas constant

        // Loop the nodes to calculate the non-conservative magnitudes
        array_1d<double,3> aux_vel;
        block_for_each(r_model_part.Nodes(), aux_vel,[&] (Node<3> &rNode, array_1d<double,3>&rVelocity) {
            const auto& r_mom = rNode.FastGetSolutionStepValue(MOMENTUM);
            const double& r_rho = rNode.FastGetSolutionStepValue(DENSITY);
            const double& r_tot_ener = rNode.FastGetSolutionStepValue(TOTAL_ENERGY);
            rVelocity = r_mom / r_rho;
            const double temp = (r_tot_ener / r_rho - 0.5 * inner_prod(rVelocity, rVelocity)) / c_v;
            const double sound_velocity = std::sqrt(gamma * R * temp);
            rNode.FastGetSolutionStepValue(VELOCITY) = rVelocity;
            rNode.FastGetSolutionStepValue(TEMPERATURE) = temp;
            rNode.FastGetSolutionStepValue(PRESSURE) = r_rho * R * temp;
            rNode.GetValue(SOUND_VELOCITY) = sound_velocity;
            rNode.GetValue(MACH) = norm_2(rVelocity) / sound_velocity;
        });
    }

    /**
     * @brief Appy the slip condition
     * This method substracts the normal projection of the momentum for all the nodes flagged as SLIP
     * The correction is computed as m_slip = m - (m·n) x n. It is intended to be called after each RK substep
     */
    void ApplySlipCondition()
    {
        auto &r_model_part = BaseType::GetModelPart();

        block_for_each(r_model_part.Nodes(), [](Node<3>& rNode){
            if (rNode.Is(SLIP)) {
                auto unit_normal = rNode.FastGetSolutionStepValue(NORMAL);
                unit_normal /= norm_2(unit_normal);
                auto& r_mom = rNode.FastGetSolutionStepValue(MOMENTUM);
                const double r_mom_n = inner_prod(r_mom, unit_normal);
                noalias(r_mom) -= r_mom_n * unit_normal;
            }
        });
    }

    void StashDiffusiveConstants()
    {
        auto& model_part = BaseType::GetModelPart();
        if(model_part.NumberOfElements() == 0) return;

        auto& properties = model_part.ElementsBegin()->GetProperties();

        auto& r_conductivity = properties.GetValue(CONDUCTIVITY);
        auto& r_dynamic_viscosity = properties.GetValue(DYNAMIC_VISCOSITY);

        mDiffusionStash.conductivity = r_conductivity;
        mDiffusionStash.dynamic_viscosity = r_dynamic_viscosity;

        r_conductivity = 0;
        r_dynamic_viscosity = 0;
    }

    void PopDiffusiveConstants()
    {
        auto& model_part = BaseType::GetModelPart();
        if(model_part.NumberOfElements() == 0) return;

        auto& properties = model_part.ElementsBegin()->GetProperties();

        properties.SetValue(CONDUCTIVITY, mDiffusionStash.conductivity);
        properties.SetValue(DYNAMIC_VISCOSITY, mDiffusionStash.dynamic_viscosity);

        mDiffusionStash.conductivity = 0;
        mDiffusionStash.dynamic_viscosity = 0;
    }


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
};

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_SOLVING_STRATEGY_BFECC  defined */