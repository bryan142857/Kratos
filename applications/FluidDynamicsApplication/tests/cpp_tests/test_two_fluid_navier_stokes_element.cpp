//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez, Ruben Zorrilla, Uxue Chasco
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "includes/global_pointer_variables.h"
#include "custom_elements/two_fluid_navier_stokes.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"

#include "processes/find_nodal_neighbours_process.h"
#include "utilities/normal_calculation_utils.h"

namespace Kratos {
    namespace Testing {

        typedef ModelPart::IndexType									 IndexType;
        typedef ModelPart::NodeIterator					          NodeIteratorType;

        /** Checks the TwoFluidNavierStokes2D3N element.
         * Checks the LHS and RHS computation
         */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes2D3N, FluidDynamicsApplicationFastSuite)
        {

            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);
            

            // Process info creation
            double delta_time = 0.1;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0/(2.0*delta_time);
            bdf_coefs[1] = -2.0/delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR,0.0);


            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
            modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);

            // Define the nodal values
            Matrix vel_original(3,2);
            vel_original(0,0) = 0.0; vel_original(0,1) = 0.1;
            vel_original(1,0) = 0.1; vel_original(1,1) = 0.2;
            vel_original(2,0) = 0.2; vel_original(2,1) = 0.3;

            // Set the nodal DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
            }

            for(unsigned int i=0; i<3; i++){
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
                for(unsigned int k=0; k<2; k++){
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) =  0.5;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(9);
            Matrix LHS = ZeroMatrix(9,9);

            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)

            KRATOS_CHECK_NEAR(RHS(0), -32.05045961 , 1e-7);
            KRATOS_CHECK_NEAR(RHS(1), 30.58587996 , 1e-7);
            KRATOS_CHECK_NEAR(RHS(2), -0.201073848 , 1e-7);
            KRATOS_CHECK_NEAR(RHS(3), 103.3405559 , 1e-7);
            KRATOS_CHECK_NEAR(RHS(4), 100.0327221 , 1e-7);
            KRATOS_CHECK_NEAR(RHS(5), 0.06718133686 , 1e-7);
            KRATOS_CHECK_NEAR(RHS(6), 86.98930477 , 1e-7);
            KRATOS_CHECK_NEAR(RHS(7), 150.7991513 , 1e-7);
            KRATOS_CHECK_NEAR(RHS(8), -0.01610748885 , 1e-7);
        }

        // /** Checks the TwoFluidNavierStokes3D4N element
        //  * Checks the LHS and RHS for a cut element
        //  */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesCut3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);

            // Process info creation
            double delta_time = 0.1;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0/(2.0*delta_time);
            bdf_coefs[1] = -2.0/delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);

            // Define the nodal values
            Matrix vel_original(4,3);
            vel_original(0,0) = 0.0; vel_original(0,1) = 0.1; vel_original(0,2) = 0.2;
            vel_original(1,0) = 0.1; vel_original(1,1) = 0.2; vel_original(1,2) = 0.3;
            vel_original(2,0) = 0.2; vel_original(2,1) = 0.3; vel_original(2,2) = 0.4;
            vel_original(3,0) = 0.3; vel_original(3,1) = 0.4; vel_original(3,2) = 0.5;

            // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
            }

            for(unsigned int i=0; i<4; i++){
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
                for(unsigned int k=0; k<3; k++){
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) =  1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) =  1.0;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16,16);

            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);
            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            KRATOS_CHECK_NEAR(RHS(0), -99.67824766, 1e-7);
            KRATOS_CHECK_NEAR(RHS(1), 20.15785546, 1e-7);
            KRATOS_CHECK_NEAR(RHS(2), -231.8210039, 1e-7);
            KRATOS_CHECK_NEAR(RHS(3), -0.05448083831, 1e-7);
            KRATOS_CHECK_NEAR(RHS(4), 61.76382984, 1e-7);
            KRATOS_CHECK_NEAR(RHS(5), 31.04354526, 1e-7);
            KRATOS_CHECK_NEAR(RHS(6), -445.5102133, 1e-7);
            KRATOS_CHECK_NEAR(RHS(7), 0.1577814843, 1e-7);
            KRATOS_CHECK_NEAR(RHS(8), 78.28977281, 1e-7);
            KRATOS_CHECK_NEAR(RHS(9), 20.92404867, 1e-7);
            KRATOS_CHECK_NEAR(RHS(10), -435.4151357, 1e-7);
            KRATOS_CHECK_NEAR(RHS(11), 0.006595606811, 1e-7);
            KRATOS_CHECK_NEAR(RHS(12), 115.0227023, 1e-7);
            KRATOS_CHECK_NEAR(RHS(13), 41.93608137, 1e-7);
            KRATOS_CHECK_NEAR(RHS(14), -488.7122322, 1e-7);
            KRATOS_CHECK_NEAR(RHS(15), -0.2098962529, 1e-7);

        }

        // /** Checks the TwoFluidNavierStokes3D4N element
        //  * Checks the LHS and RHS for a cut element with surface tension
        //  */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesSurfaceTension3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);

            // Process info creation
            double delta_time = 0.1;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            modelPart.GetProcessInfo().SetValue(SURFACE_TENSION, true);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0/(2.0*delta_time);
            bdf_coefs[1] = -2.0/delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            pElemProp->SetValue(SURFACE_TENSION_COEFFICIENT, 1.0);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);

            // Define the nodal values
            Matrix vel_original(4,3);
            vel_original(0,0) = 0.0; vel_original(0,1) = 0.1; vel_original(0,2) = 0.2;
            vel_original(1,0) = 0.1; vel_original(1,1) = 0.2; vel_original(1,2) = 0.3;
            vel_original(2,0) = 0.2; vel_original(2,1) = 0.3; vel_original(2,2) = 0.4;
            vel_original(3,0) = 0.3; vel_original(3,1) = 0.4; vel_original(3,2) = 0.5;

            // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
            }

            for(unsigned int i=0; i<4; i++){
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
                for(unsigned int k=0; k<3; k++){
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }

            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) =  1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) =  1.0;

            pElement->GetGeometry()[0].SetValue(CURVATURE, 1.0);
            pElement->GetGeometry()[1].SetValue(CURVATURE, 1.0);
            pElement->GetGeometry()[2].SetValue(CURVATURE, 1.0);
            pElement->GetGeometry()[3].SetValue(CURVATURE, 1.0);

            Vector inner_pressure_grad = ZeroVector(3);
            Vector outer_pressure_grad = ZeroVector(3);

            for(unsigned int k=0; k<3; k++){
                inner_pressure_grad[k] = 0.5*k;
                outer_pressure_grad[k] = 1.5*k + 1.0;
            }

            pElement->GetGeometry()[0].SetValue(PRESSURE_GRADIENT, inner_pressure_grad);
            pElement->GetGeometry()[1].SetValue(PRESSURE_GRADIENT, outer_pressure_grad);
            pElement->GetGeometry()[2].SetValue(PRESSURE_GRADIENT, inner_pressure_grad);
            pElement->GetGeometry()[3].SetValue(PRESSURE_GRADIENT, outer_pressure_grad);

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Vector reference_RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16,16);

            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);
            
            reference_RHS[0] = -202.5894191;
            reference_RHS[1] = 2.719387838;
            reference_RHS[2]= -148.7936481;
            reference_RHS[3]= -0.07197384163;
            reference_RHS[4]= 116.6707765;
            reference_RHS[5]= 36.92585031;
            reference_RHS[6]= -489.8190839;
            reference_RHS[7]= 0.3700935082;
            reference_RHS[8]= 164.9620866;
            reference_RHS[9]= 33.24122903;
            reference_RHS[10]= -522.5500755;
            reference_RHS[11]= 0.02197384163;
            reference_RHS[12]= 203.4918275;
            reference_RHS[13]= 50.40200447;
            reference_RHS[14]= -566.6641308;
            reference_RHS[15]= -0.4200935082;

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            KRATOS_CHECK_VECTOR_NEAR(RHS, reference_RHS, 1e-7);

        }

        // /** Checks the TwoFluidNavierStokes3D4N element
        //  * Checks the LHS and RHS for a negative element (distance <= 0.0)
        //  */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesNegativeSide3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);

            // Process info creation
            double delta_time = 0.1;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0 / (2.0*delta_time);
            bdf_coefs[1] = -2.0 / delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3, 4 };
            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);

            // Define the nodal values
            Matrix vel_original(4, 3);
            vel_original(0, 0) = 0.0; vel_original(0, 1) = 0.1; vel_original(0, 2) = 0.2;
            vel_original(1, 0) = 0.1; vel_original(1, 1) = 0.2; vel_original(1, 2) = 0.3;
            vel_original(2, 0) = 0.2; vel_original(2, 1) = 0.3; vel_original(2, 2) = 0.4;
            vel_original(3, 0) = 0.3; vel_original(3, 1) = 0.4; vel_original(3, 2) = 0.5;

            // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node = modelPart.NodesBegin(); it_node < modelPart.NodesEnd(); ++it_node) {
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
            }

            for (unsigned int i = 0; i < 4; i++) {
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 3; k++) {
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = -1.0;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16, 16);

            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            KRATOS_CHECK_NEAR(RHS(0),4.375729989, 1e-7);
            KRATOS_CHECK_NEAR(RHS(1),2.636918779, 1e-7);
            KRATOS_CHECK_NEAR(RHS(2),58.47698649, 1e-7);
            KRATOS_CHECK_NEAR(RHS(3),0.6747571137, 1e-7);
            KRATOS_CHECK_NEAR(RHS(4),6.725748185, 1e-7);
            KRATOS_CHECK_NEAR(RHS(5),36.96298944, 1e-7);
            KRATOS_CHECK_NEAR(RHS(6),-599.6096953, 1e-7);
            KRATOS_CHECK_NEAR(RHS(7),-0.003324893955, 1e-7);
            KRATOS_CHECK_NEAR(RHS(8),21.62210894, 1e-7);
            KRATOS_CHECK_NEAR(RHS(9),33.30213947, 1e-7);
            KRATOS_CHECK_NEAR(RHS(10),-665.8383544, 1e-7);
            KRATOS_CHECK_NEAR(RHS(11),0.0220714121, 1e-7);
            KRATOS_CHECK_NEAR(RHS(12),26.01064318, 1e-7);
            KRATOS_CHECK_NEAR(RHS(13),50.44496624, 1e-7);
            KRATOS_CHECK_NEAR(RHS(14),-744.6519492, 1e-7);
            KRATOS_CHECK_NEAR(RHS(15), -0.79350363, 1e-7);

        }

        // /** Checks the TwoFluidNavierStokes3D4N element
        //  * Checks the LHS and RHS for a positive element (distance > 0.0)
        //  */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesPositiveSide3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);

            // Process info creation
            double delta_time = 0.1;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0 / (2.0*delta_time);
            bdf_coefs[1] = -2.0 / delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3, 4 };
            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);

            // Define the nodal values
            Matrix vel_original(4, 3);
            vel_original(0, 0) = 0.0; vel_original(0, 1) = 0.1; vel_original(0, 2) = 0.2;
            vel_original(1, 0) = 0.1; vel_original(1, 1) = 0.2; vel_original(1, 2) = 0.3;
            vel_original(2, 0) = 0.2; vel_original(2, 1) = 0.3; vel_original(2, 2) = 0.4;
            vel_original(3, 0) = 0.3; vel_original(3, 1) = 0.4; vel_original(3, 2) = 0.5;

            // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node = modelPart.NodesBegin(); it_node < modelPart.NodesEnd(); ++it_node) {
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
            }

            for (unsigned int i = 0; i < 4; i++) {
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 3; k++) {
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = 1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;
            pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = 1.0;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16, 16);

            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            KRATOS_CHECK_NEAR(RHS(0), 4.375729989, 1e-7);
            KRATOS_CHECK_NEAR(RHS(1), 2.636918779, 1e-7);
            KRATOS_CHECK_NEAR(RHS(2), 58.47698649, 1e-7);
            KRATOS_CHECK_NEAR(RHS(3), 0.6747571137, 1e-7);
            KRATOS_CHECK_NEAR(RHS(4), 6.725748185, 1e-7);
            KRATOS_CHECK_NEAR(RHS(5), 36.96298944, 1e-7);
            KRATOS_CHECK_NEAR(RHS(6), -599.6096953, 1e-7);
            KRATOS_CHECK_NEAR(RHS(7), -0.003324893955, 1e-7);
            KRATOS_CHECK_NEAR(RHS(8), 21.62210894, 1e-7);
            KRATOS_CHECK_NEAR(RHS(9), 33.30213947, 1e-7);
            KRATOS_CHECK_NEAR(RHS(10), -665.8383544, 1e-7);
            KRATOS_CHECK_NEAR(RHS(11), 0.0220714121, 1e-7);
            KRATOS_CHECK_NEAR(RHS(12), 26.01064318, 1e-7);
            KRATOS_CHECK_NEAR(RHS(13), 50.44496624, 1e-7);
            KRATOS_CHECK_NEAR(RHS(14), -744.6519492, 1e-7);
            KRATOS_CHECK_NEAR(RHS(15), -0.7935036318, 1e-7);

        }

        /** Checks the TwoFluidNavierStokes2D3N element in a hydrostatic case.
         *  Checks the computation of the RHS
         */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes2D3NHydrostatic, FluidDynamicsApplicationFastSuite)
        {

            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);

            // Process info creation
            double delta_time = 0.1;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);

            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0/(2.0*delta_time);
            bdf_coefs[1] = -2.0/delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-03);
            Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 2.0, 0.0, 0.0);  // 0 = node 1
            modelPart.CreateNewNode(2, 2.0, 2.0, 0.0);	// 1 = node 2
            modelPart.CreateNewNode(3, 0.0, 2.0, 0.0);	// 2 = node 3

            std::vector<ModelPart::IndexType> elemNodes1 {1, 2, 3};

            modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 1, elemNodes1, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);

            // Define the nodal values as 0 for hydrostatic case
            Matrix vel_original(3,2);
            vel_original(0,0) = 0.0; vel_original(0,1) = 0.0;
            vel_original(1,0) = 0.0; vel_original(1,1) = 0.0;
            vel_original(2,0) = 0.0; vel_original(2,1) = 0.0;

            // Setting equal nodal values for DENSITY, DYNAMIC_VISCOSITY, BODY_FORCE
            for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_X) = 0.0;
                it_node->FastGetSolutionStepValue(BODY_FORCE_Y) = -10.0;
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = 0.0;
            }

            // Setting the density (different for nodes since element cut by surface)
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DENSITY) = 2.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DENSITY) = 1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DENSITY) = 1.0;

            for(unsigned int i=0; i<3; i++){
                for(unsigned int k=0; k<2; k++){
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = vel_original(i,k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = vel_original(i,k);
                    // pElement->GetGeometry()[i].Fix(VELOCITY);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }

            pElement->GetGeometry()[0].Fix(VELOCITY_X);
            pElement->GetGeometry()[0].Fix(VELOCITY_Y);

            // Setting the density (different for nodes to define the position of the surface)
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;

            // Simon : Setting the pressure
            pElement->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) = 30.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) = 0.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE) = 0.0;

            pElement->GetGeometry()[0].Fix(PRESSURE);

            // Compute RHS and LHS
            Vector RHS = ZeroVector(9);
            Matrix LHS = ZeroMatrix(9,9);

            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            double det;
            MathUtils<double>::InvertMatrix(LHS, LHS, det);

            const Vector solVec = prod(LHS, RHS);

            // The remaining residuals in the velocities have the size of the boundary integrals over the enriched pressure.
            // If the "standard" pressure shape functions are used, the results do not hold.
            KRATOS_CHECK_NEAR(RHS(0), 0.0, 1e-7);		// U_x at node 1
            KRATOS_CHECK_NEAR(RHS(1), -17.5, 1e-7); 	// U_y at node 1
            KRATOS_CHECK_NEAR(RHS(2), 0.0, 1e-7);		// P   at node 1

            KRATOS_CHECK_NEAR(RHS(3), 7.5, 1e-7);		// U_x at node 2
            KRATOS_CHECK_NEAR(RHS(4), 0.0, 1e-7);		// U_y at node 2
            KRATOS_CHECK_NEAR(RHS(5), 0.0, 1e-7);		// P   at node 2

            KRATOS_CHECK_NEAR(RHS(6), -7.5, 1e-7);		// U_x at node 3
            KRATOS_CHECK_NEAR(RHS(7), -7.5, 1e-7);		// U_y at node 3
            KRATOS_CHECK_NEAR(RHS(8), 0.0, 1e-7);		// P   at node 3
        }




        /** Includes the TwoFluidNavierStokes2D3N element with the BEHR2004 slip boundary condition in a hydrostatic case.
         *  The BEHR2004 contribution is activated, if the flag SLIP is set for the NavierStokesWallCondition2D2N wall condition
         *  Checks the computation of the RHS for components in tangential direction
         */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes2D3NHydrostaticBehr, FluidDynamicsApplicationFastSuite){

            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);
            modelPart.AddNodalSolutionStepVariable(NORMAL);
            modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);

            // Process info creation
            double delta_time = 0.01;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.0);
            modelPart.GetProcessInfo().SetValue(EXTERNAL_PRESSURE, 0.0);
            modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0/(2.0*delta_time);
            bdf_coefs[1] = -2.0/delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 2.0, 0.0);
            modelPart.CreateNewNode(2, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 1.0, 0.0, 0.0);

            // Creation of elements
            std::vector<ModelPart::IndexType> elemNodes1 {1, 2, 3};
            std::vector<ModelPart::IndexType> elemNodes2 {2, 4, 3};
            modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 1, elemNodes1, pElemProp);
            modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 2, elemNodes2, pElemProp);
            Element::Pointer pElement1 = modelPart.pGetElement(1);
            Element::Pointer pElement2 = modelPart.pGetElement(2);

            std::vector<ModelPart::IndexType> condNodes1 {1, 2};
            std::vector<ModelPart::IndexType> condNodes2 {2, 4};
            modelPart.CreateNewCondition("NavierStokesWallCondition2D2N", 1, condNodes1, pElemProp);
            modelPart.CreateNewCondition("NavierStokesWallCondition2D2N", 2, condNodes2, pElemProp);
            Condition::Pointer pCondition1 = modelPart.pGetCondition(1);
            Condition::Pointer pCondition2 = modelPart.pGetCondition(2);

            pCondition1->SetFlags(SLIP);
            pCondition2->SetFlags(SLIP);

            // artificially assigning parents (regularly done by check_and_prepare_model_part_process)
            GlobalPointersVector<Element> wpParent1;
            wpParent1.push_back(GlobalPointer<Element>(modelPart.Elements()(1).get()));
            pCondition1->SetValue( NEIGHBOUR_ELEMENTS, wpParent1 );

            GlobalPointersVector<Element> wpParent2;
            wpParent2.push_back(GlobalPointer<Element>(modelPart.Elements()(2).get()));
            pCondition2->SetValue( NEIGHBOUR_ELEMENTS, wpParent2 );

            Vector elemRHS1 = ZeroVector(9);
            Vector elemRHS2 = ZeroVector(9);
            Matrix elemLHS = ZeroMatrix(9,9);

            Vector condRHS1 = ZeroVector(6);
            Vector condRHS2 = ZeroVector(6);
            Matrix condLHS = ZeroMatrix(6,6);

            for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
                it_node->FastGetSolutionStepValue(DENSITY) = 1000.0;
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Y) = -10.0;
                it_node->SetValue(SLIP_LENGTH, 1.0e10);
            }

            for(unsigned int i=0; i<3; i++){
                for(unsigned int k=0; k<2; k++){
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = 0.0;
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.0;
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.0;
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }

            pElement1->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    = 10000.0;
            pElement1->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    = 20000.0;
            pElement1->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)    = 20000.0;

            for(unsigned int i=0; i<3; i++){
                for(unsigned int k=0; k<2; k++){
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = 0.0;
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.0;
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.0;
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }

            pElement2->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    = 20000.0;
            pElement2->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    = 30000.0;
            pElement2->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)    = 20000.0;

            pElement1->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement1->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -2.0;
            pElement1->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -2.0;
            pElement2->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -2.0;
            pElement2->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -3.0;
            pElement2->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -2.0;

            // Initialization
            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement1->Initialize(r_process_info);
            pElement2->Initialize(r_process_info);

            pElement1->CalculateLocalSystem(elemLHS, elemRHS1, r_process_info);
            pElement2->CalculateLocalSystem(elemLHS, elemRHS2, r_process_info);

            FindNodalNeighboursProcess find_nodal_neighbours_process(modelPart);
            find_nodal_neighbours_process.Execute();

            NormalCalculationUtils find_nodal_normal_utility;
            find_nodal_normal_utility.CalculateOnSimplex(modelPart, 2);

            pCondition1->Initialize(r_process_info);
            pCondition2->Initialize(r_process_info);

            // Computing locel contributions
            pCondition1->CalculateLocalSystem(condLHS, condRHS1, r_process_info);
            pCondition2->CalculateLocalSystem(condLHS, condRHS2, r_process_info);

            // Assembly of the residual for node 2 (node between the conditions)
            Vector contriFromElem1 = ZeroVector(3);
            Vector contriFromElem2 = ZeroVector(3);
            Vector contriFromCond1 = ZeroVector(3);
            Vector contriFromCond2 = ZeroVector(3);

            for (unsigned i = 0; i < 3; i++){
                contriFromElem1[i] = elemRHS1[3+i];
                contriFromElem2[i] = elemRHS2[0+i];

                contriFromCond1[i] = condRHS1[3+i];
                contriFromCond2[i] = condRHS2[0+i];
            }

            Vector residualAtNodeTwo = contriFromElem1 + contriFromElem2 + contriFromCond1 + contriFromCond2;
            Vector normalAtNodeTwo = pElement1->GetGeometry()[1].FastGetSolutionStepValue(NORMAL);

            Vector tangentialComponent;
            array_1d<double,3> residualAtNodeTwoVector;
            array_1d<double,3> normalAtNodeTwoVector;

            tangentialComponent = MathUtils<double>::CrossProduct( normalAtNodeTwo, residualAtNodeTwo );
            tangentialComponent = - MathUtils<double>::CrossProduct( normalAtNodeTwo, tangentialComponent );

            KRATOS_CHECK_NEAR( tangentialComponent[0], 0.0, 1e-7);
            KRATOS_CHECK_NEAR( tangentialComponent[1], 0.0, 1e-7);
            KRATOS_CHECK_NEAR( tangentialComponent[2], 0.0, 1e-7);
        }


        /** Includes 3 elements of TwoFluidNavierStokes3D4N element
         *  and 3 conditions of type BEHR2004 slip boundary condition in a hydrostatic case.
         *  The BEHR2004 contribution is activated, if the flag SLIP is set for the NavierStokesWallCondition3D3N wall condition
         *  Checks the computation of the RHS for components in tangential direction
         */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes3D4NHydrostaticBehr, FluidDynamicsApplicationFastSuite){

            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);
            modelPart.AddNodalSolutionStepVariable(NORMAL);
            modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);

            // Process info creation
            double delta_time = 0.01;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.0);
            modelPart.GetProcessInfo().SetValue(EXTERNAL_PRESSURE, 0.0);
            modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0/(2.0*delta_time);
            bdf_coefs[1] = -2.0/delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 0.0001);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation (following MDPA file from GiD)
            // if y coordinate is altered, the pressure must be adapted (all other: random values)
            modelPart.CreateNewNode(1,  2.0,  0.0,  -1.0);		// y = 0
            modelPart.CreateNewNode(2, -1.0,  7.0,  2.0);		// y = 7
            modelPart.CreateNewNode(3,  0.0,  5.0,  5.0);		// y = 5  (will be used to check)
            modelPart.CreateNewNode(4, -3.0, 10.0,  3.0);		// y = 10
            modelPart.CreateNewNode(5,  6.0, 10.0,  1.0);		// y = 10

            // Creation of elements (following MDPA file from GiD)
            std::vector<ModelPart::IndexType> elemNodes1 {5, 4, 2, 3};     	// start at position 12
            std::vector<ModelPart::IndexType> elemNodes2 {3, 4, 2, 1};		// start at position 0
            std::vector<ModelPart::IndexType> elemNodes3 {3, 2, 5, 1};		// start at position 0

            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 1, elemNodes1, pElemProp);
            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 2, elemNodes2, pElemProp);
            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 3, elemNodes3, pElemProp);

            Element::Pointer pElement1 = modelPart.pGetElement(1);
            Element::Pointer pElement2 = modelPart.pGetElement(2);
            Element::Pointer pElement3 = modelPart.pGetElement(3);

            // Creation of conditions (following MDPA file from GiD)
            std::vector<ModelPart::IndexType> condNodes1 {1, 5, 3};			// start at position 8
            std::vector<ModelPart::IndexType> condNodes2 {4, 1, 3};			// start at position 8
            std::vector<ModelPart::IndexType> condNodes3 {5, 4, 3};			// start at position 8

            modelPart.CreateNewCondition("NavierStokesWallCondition3D3N", 1, condNodes1, pElemProp);
            modelPart.CreateNewCondition("NavierStokesWallCondition3D3N", 2, condNodes2, pElemProp);
            modelPart.CreateNewCondition("NavierStokesWallCondition3D3N", 3, condNodes3, pElemProp);

            Condition::Pointer pCondition1 = modelPart.pGetCondition(1);
            Condition::Pointer pCondition2 = modelPart.pGetCondition(2);
            Condition::Pointer pCondition3 = modelPart.pGetCondition(3);

            pCondition1->SetFlags(SLIP);
            pCondition2->SetFlags(SLIP);
            pCondition3->SetFlags(SLIP);

            // artificially assigning parents (regularly done by check_and_prepare_model_part_process)
            GlobalPointersVector<Element> wpParent1;
            wpParent1.push_back(GlobalPointer<Element>(modelPart.Elements()(3)));
            pCondition1->SetValue( NEIGHBOUR_ELEMENTS, wpParent1 );

            GlobalPointersVector<Element> wpParent2;
            wpParent2.push_back(GlobalPointer<Element>(modelPart.Elements()(2)));
            pCondition2->SetValue( NEIGHBOUR_ELEMENTS, wpParent2 );

            GlobalPointersVector<Element> wpParent3;
            wpParent3.push_back(GlobalPointer<Element>(modelPart.Elements()(1)));
            pCondition3->SetValue( NEIGHBOUR_ELEMENTS, wpParent3 );

            Vector elemRHS1 = ZeroVector(16);
            Vector elemRHS2 = ZeroVector(16);
            Vector elemRHS3 = ZeroVector(16);
            Matrix elemLHS = ZeroMatrix(16,16);

            Vector condRHS1 = ZeroVector(12);
            Vector condRHS2 = ZeroVector(12);
            Vector condRHS3 = ZeroVector(12);
            Matrix condLHS = ZeroMatrix(12,12);

            for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
                it_node->FastGetSolutionStepValue(DENSITY) = 1000.0;
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_X) = 0.0;
                it_node->FastGetSolutionStepValue(BODY_FORCE_Y) = -10.0;
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = 0.0;
                it_node->SetValue(SLIP_LENGTH, 1.0e10);
            }

            for(unsigned int i=0; i<4; i++){
                for (unsigned int k = 0; k < 3; k++){
                    // fixing the mesh
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;

                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;

                    pElement3->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement3->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement3->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }

            for(unsigned int timestep = 0; timestep < 3; timestep++){
                for ( unsigned int nnode = 0; nnode < 4; nnode++){
                    // setting fluid velocity to zero
                    pElement1->GetGeometry()[nnode].FastGetSolutionStepValue(VELOCITY, timestep)[0] = 0.0; 	// x
                    pElement1->GetGeometry()[nnode].FastGetSolutionStepValue(VELOCITY, timestep)[1] = 0.0;	// y
                    pElement1->GetGeometry()[nnode].FastGetSolutionStepValue(VELOCITY, timestep)[2] = 0.0;	// z

                    pElement2->GetGeometry()[nnode].FastGetSolutionStepValue(VELOCITY, timestep)[0] = 0.0; 	// x
                    pElement2->GetGeometry()[nnode].FastGetSolutionStepValue(VELOCITY, timestep)[1] = 0.0;	// y
                    pElement2->GetGeometry()[nnode].FastGetSolutionStepValue(VELOCITY, timestep)[2] = 0.0;	// z

                    pElement3->GetGeometry()[nnode].FastGetSolutionStepValue(VELOCITY, timestep)[0] = 0.0; 	// x
                    pElement3->GetGeometry()[nnode].FastGetSolutionStepValue(VELOCITY, timestep)[1] = 0.0;	// y
                    pElement3->GetGeometry()[nnode].FastGetSolutionStepValue(VELOCITY, timestep)[2] = 0.0;	// z
                }
            }


            // std::vector<ModelPart::IndexType> elemNodes1 {5, 4, 2, 3};     	// start at position 12
            // std::vector<ModelPart::IndexType> elemNodes2 {3, 4, 2, 1};		// start at position 0
            // std::vector<ModelPart::IndexType> elemNodes3 {3, 2, 5, 1};		// start at position 0

            pElement1->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    = 100000.0;
            pElement1->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    = 100000.0;
            pElement1->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)    = 130000.0;
            pElement1->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE)    = 15000.0;
            pElement1->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -10.0;
            pElement1->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -10.0;
            pElement1->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -13.0;
            pElement1->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = -15.0;

            pElement2->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    = 150000.0;
            pElement2->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    = 100000.0;
            pElement2->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)    = 130000.0;
            pElement2->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE)    = 200000.0;
            pElement2->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -15.0;
            pElement2->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -10.0;
            pElement2->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -13.0;
            pElement2->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = -20.0;

            pElement3->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    = 150000.0;
            pElement3->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    = 130000.0;
            pElement3->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)    = 100000.0;
            pElement3->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE)    = 200000.0;
            pElement3->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -15.0;
            pElement3->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -13.0;
            pElement3->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -10.0;
            pElement3->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = -20.0;

            // Assembly of the residual for node 4 (node between the 3 conditions)
            Vector contriFromElem1 = ZeroVector(3);
            Vector contriFromElem2 = ZeroVector(3);
            Vector contriFromElem3 = ZeroVector(3);

            Vector contriFromCond1 = ZeroVector(3);
            Vector contriFromCond2 = ZeroVector(3);
            Vector contriFromCond3 = ZeroVector(3);

            // Initialization
            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement1->Initialize(r_process_info);
            pElement2->Initialize(r_process_info);
            pElement3->Initialize(r_process_info);

            pElement1->InitializeSolutionStep(r_process_info);
            pElement2->InitializeSolutionStep(r_process_info);
            pElement3->InitializeSolutionStep(r_process_info);

            pElement1->CalculateLocalSystem(elemLHS, elemRHS1, r_process_info);
            pElement2->CalculateLocalSystem(elemLHS, elemRHS2, r_process_info);
            pElement3->CalculateLocalSystem(elemLHS, elemRHS3, r_process_info);

            FindNodalNeighboursProcess find_nodal_neighbours_process(modelPart);
            find_nodal_neighbours_process.Execute();

            NormalCalculationUtils find_nodal_normal_utility;
            find_nodal_normal_utility.CalculateOnSimplex(modelPart, 3);

            pCondition1->Initialize(r_process_info);
            pCondition2->Initialize(r_process_info);
            pCondition3->Initialize(r_process_info);

            // Computing local contributions of elemet
            pCondition1->CalculateLocalSystem(condLHS, condRHS1, r_process_info);
            pCondition2->CalculateLocalSystem(condLHS, condRHS2, r_process_info);
            pCondition3->CalculateLocalSystem(condLHS, condRHS3, r_process_info);

            for (unsigned int i = 0; i < 3; i++){
                contriFromElem1[i] = elemRHS1[12 + i];
                contriFromElem2[i] = elemRHS2[i];
                contriFromElem3[i] = elemRHS3[i];

                contriFromCond1[i] = condRHS1[8 + i];
                contriFromCond2[i] = condRHS2[8 + i];
                contriFromCond3[i] = condRHS3[8 + i];
            }

            Vector residualAtNodeTwo = 	contriFromElem1 + contriFromElem2 + contriFromElem3 +
                                        ( contriFromCond1 + contriFromCond2 + contriFromCond3 );

            Vector normalAtNodeTwo = pElement2->GetGeometry()[0].FastGetSolutionStepValue(NORMAL);

            double sumOfSquares = 0.0;
            for (int i = 0; i < 3; i++){
                sumOfSquares += normalAtNodeTwo[i] * normalAtNodeTwo[i];
            }

            normalAtNodeTwo /= sqrt(sumOfSquares);

            Vector tangentialComponent;
            Vector normalComponent;

            tangentialComponent = MathUtils<double>::CrossProduct( normalAtNodeTwo, residualAtNodeTwo );
            tangentialComponent = - MathUtils<double>::CrossProduct( normalAtNodeTwo, tangentialComponent );

            normalComponent = MathUtils<double>::Dot3( residualAtNodeTwo, normalAtNodeTwo ) * normalAtNodeTwo;

            KRATOS_CHECK_NEAR( tangentialComponent[0], 0.0, 1e-7);
            KRATOS_CHECK_NEAR( tangentialComponent[1], 0.0, 1e-7);
            KRATOS_CHECK_NEAR( tangentialComponent[2], 0.0, 1e-7);
        }



        /** Includes the TwoFluidNavierStokes2D3N element with the BEHR2004 slip boundary condition in a hydrostatic case.
         *  The BEHR2004 contribution is activated, if the flag SLIP is set for the NavierStokesWallCondition2D2N wall condition
         *  Checks the computation of the RHS for components in tangential direction
         */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes2D3NStressBehr, FluidDynamicsApplicationFastSuite){

            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);
            modelPart.AddNodalSolutionStepVariable(NORMAL);
            modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);

            // Process info creation
            double delta_time = 0.01;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.0);
            modelPart.GetProcessInfo().SetValue(EXTERNAL_PRESSURE, 0.0);
            modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0/(2.0*delta_time);
            bdf_coefs[1] = -2.0/delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1);
            Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 2.0, 0.0);
            modelPart.CreateNewNode(2, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 1.0, 0.0, 0.0);

            // Creation of elements
            std::vector<ModelPart::IndexType> elemNodes1 {1, 2, 3};
            std::vector<ModelPart::IndexType> elemNodes2 {2, 4, 3};
            modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 1, elemNodes1, pElemProp);
            modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 2, elemNodes2, pElemProp);
            Element::Pointer pElement1 = modelPart.pGetElement(1);
            Element::Pointer pElement2 = modelPart.pGetElement(2);


            std::vector<ModelPart::IndexType> condNodes1 {1, 2};
            std::vector<ModelPart::IndexType> condNodes2 {2, 4};
            modelPart.CreateNewCondition("NavierStokesWallCondition2D2N", 1, condNodes1, pElemProp);
            modelPart.CreateNewCondition("NavierStokesWallCondition2D2N", 2, condNodes2, pElemProp);
            Condition::Pointer pCondition1 = modelPart.pGetCondition(1);
            Condition::Pointer pCondition2 = modelPart.pGetCondition(2);

            pCondition1->SetFlags(SLIP);
            pCondition2->SetFlags(SLIP);

            // artificially assigning parents (regularly done by check_and_prepare_model_part_process)
            GlobalPointersVector<Element> wpParent1;
            wpParent1.push_back(GlobalPointer<Element>(modelPart.Elements()(1)));
            pCondition1->SetValue( NEIGHBOUR_ELEMENTS, wpParent1 );

            GlobalPointersVector<Element> wpParent2;
            wpParent2.push_back(GlobalPointer<Element>(modelPart.Elements()(2)));
            pCondition2->SetValue( NEIGHBOUR_ELEMENTS, wpParent2 );

            Vector elemRHS1 = ZeroVector(9);
            Vector elemRHS2 = ZeroVector(9);
            Matrix elemLHS = ZeroMatrix(9,9);

            Vector condRHS1 = ZeroVector(6);
            Vector condRHS2 = ZeroVector(6);
            Matrix condLHS = ZeroMatrix(6,6);

            for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
                it_node->FastGetSolutionStepValue(DENSITY) = 1000.0;
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);

                it_node->FastGetSolutionStepValue(BODY_FORCE_X) = 0.0;
                it_node->FastGetSolutionStepValue(BODY_FORCE_Y) = 0.0;
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = 0.0;
                it_node->SetValue(SLIP_LENGTH, 1.0e10);
            }

            for(unsigned int i=0; i<3; i++){
                for(unsigned int k=0; k<2; k++){
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = 0.0;
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.0;
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.0;
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }

            for(unsigned int i=0; i<3; i++){
                for(unsigned int k=0; k<2; k++){
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = 0.0;
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.0;
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.0;
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }


            for (unsigned int time = 0; time < 3; time++){
                    pElement1->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, time)[0] = 0.0;
                    pElement1->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, time)[1] = 0.0; // 1.0
                    pElement1->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, time)[2] = 0.0;

                    pElement2->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, time)[0] = 0.0;
                    pElement2->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, time)[1] = 0.0; // 1.0
                    pElement2->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, time)[2] = 0.0;
            }


            pElement1->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) = 0.0;
            pElement1->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) = 0.0;
            pElement1->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE) = 0.0;
            pElement2->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) = 0.0;
            pElement2->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) = 0.0;
            pElement2->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE) = 0.0;

            pElement1->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement1->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -2.0;
            pElement1->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -2.0;
            pElement2->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -2.0;
            pElement2->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -3.0;
            pElement2->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -2.0;

            // Initialization
            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement1->Initialize(r_process_info);
            pElement2->Initialize(r_process_info);

            pElement1->CalculateLocalSystem(elemLHS, elemRHS1, r_process_info);
            pElement2->CalculateLocalSystem(elemLHS, elemRHS2, r_process_info);

            FindNodalNeighboursProcess find_nodal_neighbours_process(modelPart);
            find_nodal_neighbours_process.Execute();

            NormalCalculationUtils find_nodal_normal_utility;
            find_nodal_normal_utility.CalculateOnSimplex(modelPart, 2);

            pCondition1->Initialize(r_process_info);
            pCondition2->Initialize(r_process_info);

            // Computing locel contributions
            pCondition1->CalculateLocalSystem(condLHS, condRHS1, r_process_info);
            pCondition2->CalculateLocalSystem(condLHS, condRHS2, r_process_info);

            // Assembly of the residual for node 2 (node between the conditions)
            Vector contriFromElem1 = ZeroVector(3);
            Vector contriFromElem2 = ZeroVector(3);
            Vector contriFromCond1 = ZeroVector(3);
            Vector contriFromCond2 = ZeroVector(3);

            for (unsigned i = 0; i < 2; i++){
                contriFromElem1[i] = elemRHS1[3+i];
                contriFromElem2[i] = elemRHS2[0+i];

                contriFromCond1[i] = condRHS1[3+i];
                contriFromCond2[i] = condRHS2[0+i];
            }

            Vector residualAtNodeTwo = contriFromElem1 + contriFromElem2 + contriFromCond1 + contriFromCond2;
            Vector normalAtNodeTwo = pElement1->GetGeometry()[1].FastGetSolutionStepValue(NORMAL);

            Vector tangentialComponent;
            array_1d<double,3> residualAtNodeTwoVector;
            array_1d<double,3> normalAtNodeTwoVector;

            tangentialComponent = MathUtils<double>::CrossProduct( normalAtNodeTwo, residualAtNodeTwo );
            tangentialComponent = - MathUtils<double>::CrossProduct( normalAtNodeTwo, tangentialComponent );

            KRATOS_CHECK_NEAR( tangentialComponent[0], 0.0, 1e-7);
            KRATOS_CHECK_NEAR( tangentialComponent[1], 0.0, 1e-7);
            KRATOS_CHECK_NEAR( tangentialComponent[2], 0.0, 1e-7);
        }


        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesDarcy3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);

            // Process info creation
            double delta_time = 0.1;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0 / (2.0*delta_time);
            bdf_coefs[1] = -2.0 / delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(LIN_DARCY_COEF, 1.0 / 4.339E-08);
            pElemProp->SetValue(NONLIN_DARCY_COEF, 1.0 / 5.086E-04);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3, 4 };
            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);

            // Define the nodal values
            Matrix vel_original(4, 3);
            vel_original(0, 0) = 0.0; vel_original(0, 1) = 0.1; vel_original(0, 2) = 0.2;
            vel_original(1, 0) = 0.1; vel_original(1, 1) = 0.2; vel_original(1, 2) = 0.3;
            vel_original(2, 0) = 0.2; vel_original(2, 1) = 0.3; vel_original(2, 2) = 0.4;
            vel_original(3, 0) = 0.3; vel_original(3, 1) = 0.4; vel_original(3, 2) = 0.5;

            // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node = modelPart.NodesBegin(); it_node < modelPart.NodesEnd(); ++it_node) {
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
            }

            for (unsigned int i = 0; i < 4; i++) {
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 3; k++) {
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = 1.0;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16, 16);

            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            KRATOS_CHECK_NEAR(RHS(0), -357.8243200052, 1e-7);
            KRATOS_CHECK_NEAR(RHS(1), 2616.9094454176, 1e-7);
            KRATOS_CHECK_NEAR(RHS(2), -5245.9530366262, 1e-7);
            KRATOS_CHECK_NEAR(RHS(3), -0.013507842388113, 1e-7);
            KRATOS_CHECK_NEAR(RHS(4), -14178.358643109, 1e-7);
            KRATOS_CHECK_NEAR(RHS(5), -9695.4830606294, 1e-7);
            KRATOS_CHECK_NEAR(RHS(6), -9397.4044901296, 1e-7);
            KRATOS_CHECK_NEAR(RHS(7), -0.016334635022297, 1e-7);
            KRATOS_CHECK_NEAR(RHS(8), -7384.0299128137, 1e-7);
            KRATOS_CHECK_NEAR(RHS(9), -11452.419170853, 1e-7);
            KRATOS_CHECK_NEAR(RHS(10), -12294.517818084, 1e-7);
            KRATOS_CHECK_NEAR(RHS(11), -0.03635048354971, 1e-7);
            KRATOS_CHECK_NEAR(RHS(12), -7158.5516510797, 1e-7);
            KRATOS_CHECK_NEAR(RHS(13), -9706.5109782994, 1e-7);
            KRATOS_CHECK_NEAR(RHS(14), -16434.188680389, 1e-7);
            KRATOS_CHECK_NEAR(RHS(15), -0.03380703903988, 1e-7);

        }

        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes2D3NNavierSlip, FluidDynamicsApplicationFastSuite){

            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);
            modelPart.AddNodalSolutionStepVariable(NORMAL);
            modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);

            // Process info creation
            double delta_time = 0.01;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.0);
            modelPart.GetProcessInfo().SetValue(EXTERNAL_PRESSURE, 0.0);
            modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0/(2.0*delta_time);
            bdf_coefs[1] = -2.0/delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0);
            Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(2, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 2.0, 0.0, 0.0);
            modelPart.CreateNewNode(4, 4.0, 0.0, 0.0);
            modelPart.CreateNewNode(5, 4.0, 1.0, 0.0);

            // Creation of elements
            std::vector<ModelPart::IndexType> elemNodes1 {1, 2, 3};
            std::vector<ModelPart::IndexType> elemNodes2 {1, 3, 5};
            std::vector<ModelPart::IndexType> elemNodes3 {5, 3, 4};
            modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 1, elemNodes1, pElemProp);
            modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 2, elemNodes2, pElemProp);
            modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 3, elemNodes3, pElemProp);
            Element::Pointer pElement1 = modelPart.pGetElement(1);
            Element::Pointer pElement2 = modelPart.pGetElement(2);
            Element::Pointer pElement3 = modelPart.pGetElement(3);


            std::vector<ModelPart::IndexType> condNodes1 {2, 3};
            std::vector<ModelPart::IndexType> condNodes2 {3, 4};
            modelPart.CreateNewCondition("NavierStokesWallCondition2D2N", 1, condNodes1, pElemProp);
            modelPart.CreateNewCondition("NavierStokesWallCondition2D2N", 2, condNodes2, pElemProp);
            Condition::Pointer pCondition1 = modelPart.pGetCondition(1);
            Condition::Pointer pCondition2 = modelPart.pGetCondition(2);

            pCondition1->SetFlags(SLIP);
            pCondition2->SetFlags(SLIP);

            // artificially assigning parents (regularly done by check_and_prepare_model_part_process)
            GlobalPointersVector<Element> wpParent1;
            wpParent1.push_back(GlobalPointer<Element>(modelPart.Elements()(1)));
            pCondition1->SetValue( NEIGHBOUR_ELEMENTS, wpParent1 );

            GlobalPointersVector<Element> wpParent2;
            wpParent2.push_back(GlobalPointer<Element>(modelPart.Elements()(3)));
            pCondition2->SetValue( NEIGHBOUR_ELEMENTS, wpParent2 );

            Vector elemRHS1 = ZeroVector(9);
            Vector elemRHS2 = ZeroVector(9);
            Vector elemRHS3 = ZeroVector(9);
            Matrix elemLHS = ZeroMatrix(9,9);

            Vector condRHS1Before = ZeroVector(6);
            Vector condRHS2Before = ZeroVector(6);

            Vector condRHS1After = ZeroVector(6);
            Vector condRHS2After = ZeroVector(6);

            Matrix condLHS = ZeroMatrix(6,6);

            for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
                it_node->FastGetSolutionStepValue(DENSITY) = 1.0;
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);

                it_node->FastGetSolutionStepValue(BODY_FORCE_X) = 0.0;
                it_node->FastGetSolutionStepValue(BODY_FORCE_Y) = 0.0;
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = 0.0;
                it_node->SetValue(SLIP_LENGTH, 1.0e10);
            }

            for(unsigned int i=0; i<3; i++){
                for(unsigned int k=0; k<3; k++){
                    for(unsigned int timestep=0; timestep<3; timestep++){
                        // remark: This flow field HAS a wall-normal component.
                        // This choice is made deliberately to see if the tangential projection works
                        pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, timestep)[k] = 1.0;
                        pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, timestep)[k] = 1.0;
                        pElement3->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, timestep)[k] = 1.0;

                        pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, timestep)[k] = 0.0;
                        pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, timestep)[k] = 0.0;
                        pElement3->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, timestep)[k] = 0.0;
                    }
                }
            }

            for(unsigned int i=0; i<3; i++){
                pElement1->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                pElement2->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                pElement3->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
            }

            pElement1->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement1->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -2.0;
            pElement1->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -2.0;
            pElement2->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -2.0;
            pElement2->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -3.0;
            pElement2->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -2.0;

            // Initialization
            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement1->Initialize(r_process_info);
            pElement2->Initialize(r_process_info);
            pElement3->Initialize(r_process_info);

            pElement1->CalculateLocalSystem(elemLHS, elemRHS1, r_process_info);
            pElement2->CalculateLocalSystem(elemLHS, elemRHS2, r_process_info);
            pElement3->CalculateLocalSystem(elemLHS, elemRHS3, r_process_info);

            FindNodalNeighboursProcess find_nodal_neighbours_process(modelPart);
            find_nodal_neighbours_process.Execute();

            NormalCalculationUtils find_nodal_normal_utility;
            find_nodal_normal_utility.CalculateOnSimplex(modelPart, 2);

            pCondition1->Initialize(r_process_info);
            pCondition2->Initialize(r_process_info);

            // before adding Navier Slip
            pCondition1->CalculateLocalSystem(condLHS, condRHS1Before, r_process_info);
            pCondition2->CalculateLocalSystem(condLHS, condRHS2Before, r_process_info);

            // adding a considerably small slip length at only one node (node 3)
            const double navier_slip_length = 0.0001;
            pElement1->GetGeometry()[2].SetValue( SLIP_LENGTH, navier_slip_length );

            pCondition1->CalculateLocalSystem(condLHS, condRHS1After, r_process_info);
            pCondition2->CalculateLocalSystem(condLHS, condRHS2After, r_process_info);

            // moniotring changes (after - before) in the residual
            // Here: no change expected
            const double res_x_node2 = ( condRHS1After[0] - condRHS1Before[0] );
            const double res_y_node2 = ( condRHS1After[1] - condRHS1Before[1] );
            const double res_p_node2 = ( condRHS1After[2] - condRHS1Before[2] );
            const double res_x_node4 = ( condRHS2After[3] - condRHS2Before[3] );
            const double res_y_node4 = ( condRHS2After[4] - condRHS2Before[4] );
            const double res_p_node4 = ( condRHS2After[5] - condRHS2Before[5] );
            const double res_y_node3 = ( ( condRHS1After[4] + condRHS2After[1] ) - ( condRHS1Before[4] + condRHS2Before[1] ) );
            const double res_p_node3 = ( ( condRHS1After[5] + condRHS2After[2] ) - ( condRHS1Before[5] + condRHS2Before[2] ) );
            // Here: change expected only at the x-component of the node with the altered slip length
            const double res_x_node3 = ( ( condRHS1After[3] + condRHS2After[0] ) - ( condRHS1Before[3] + condRHS2Before[0] ) );

            KRATOS_CHECK_NEAR( res_x_node2, 0.0, 1e-7);
            KRATOS_CHECK_NEAR( res_y_node2, 0.0, 1e-7);
            KRATOS_CHECK_NEAR( res_p_node2, 0.0, 1e-7);

            KRATOS_CHECK_NEAR( res_x_node3, - (0.5 * 4.0 * 1.0 * 1.0 * 1.0/navier_slip_length), 1e-7 );
            KRATOS_CHECK_NEAR( res_y_node3, 0.0, 1e-7);
            KRATOS_CHECK_NEAR( res_p_node3, 0.0, 1e-7);

            KRATOS_CHECK_NEAR( res_x_node4, 0.0, 1e-7);
            KRATOS_CHECK_NEAR( res_y_node4, 0.0, 1e-7);
            KRATOS_CHECK_NEAR( res_p_node4, 0.0, 1e-7);
        }


        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes3D4NNavierSlip, FluidDynamicsApplicationFastSuite){

            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);
            modelPart.AddNodalSolutionStepVariable(NORMAL);
            modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);

            // Process info creation
            double delta_time = 0.01;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.0);
            modelPart.GetProcessInfo().SetValue(EXTERNAL_PRESSURE, 0.0);
            modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0/(2.0*delta_time);
            bdf_coefs[1] = -2.0/delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR,0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation (following MDPA file from GiD)
            // Geometry configured such the normal at node 3 only with z-component
            modelPart.CreateNewNode(1,  0.0,  0.0,  0.0);		// y = 0
            modelPart.CreateNewNode(2,  0.0,  5.0,  0.0);		// y = 7
            modelPart.CreateNewNode(3,  0.0,  5.0,  5.0);		// y = 5  (will be used to check)
            modelPart.CreateNewNode(4, -5.0, 10.0,  0.0);		// y = 10
            modelPart.CreateNewNode(5,  5.0, 10.0,  0.0);		// y = 10

            // Creation of elements (following MDPA file from GiD)
            std::vector<ModelPart::IndexType> elemNodes1 {5, 4, 2, 3};     	// start at position 12
            std::vector<ModelPart::IndexType> elemNodes2 {3, 4, 2, 1};		// start at position 0
            std::vector<ModelPart::IndexType> elemNodes3 {3, 2, 5, 1};		// start at position 0

            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 1, elemNodes1, pElemProp);
            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 2, elemNodes2, pElemProp);
            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 3, elemNodes3, pElemProp);

            Element::Pointer pElement1 = modelPart.pGetElement(1);
            Element::Pointer pElement2 = modelPart.pGetElement(2);
            Element::Pointer pElement3 = modelPart.pGetElement(3);

            // Creation of conditions (following MDPA file from GiD)
            std::vector<ModelPart::IndexType> condNodes1 {1, 5, 3};			// start at position 8
            std::vector<ModelPart::IndexType> condNodes2 {4, 1, 3};			// start at position 8
            std::vector<ModelPart::IndexType> condNodes3 {5, 4, 3};			// start at position 8

            modelPart.CreateNewCondition("NavierStokesWallCondition3D3N", 1, condNodes1, pElemProp);
            modelPart.CreateNewCondition("NavierStokesWallCondition3D3N", 2, condNodes2, pElemProp);
            modelPart.CreateNewCondition("NavierStokesWallCondition3D3N", 3, condNodes3, pElemProp);

            Condition::Pointer pCondition1 = modelPart.pGetCondition(1);
            Condition::Pointer pCondition2 = modelPart.pGetCondition(2);
            Condition::Pointer pCondition3 = modelPart.pGetCondition(3);

            pCondition1->SetFlags(SLIP);
            pCondition2->SetFlags(SLIP);
            pCondition3->SetFlags(SLIP);

            // artificially assigning parents (regularly done by check_and_prepare_model_part_process)
            GlobalPointersVector<Element> wpParent1;
            wpParent1.push_back(GlobalPointer<Element>(modelPart.Elements()(3)));
            pCondition1->SetValue( NEIGHBOUR_ELEMENTS, wpParent1 );

            GlobalPointersVector<Element> wpParent2;
            wpParent2.push_back(GlobalPointer<Element>(modelPart.Elements()(2)));
            pCondition2->SetValue( NEIGHBOUR_ELEMENTS, wpParent2 );

            GlobalPointersVector<Element> wpParent3;
            wpParent3.push_back(GlobalPointer<Element>(modelPart.Elements()(1)));
            pCondition3->SetValue( NEIGHBOUR_ELEMENTS, wpParent3 );

            Vector elemRHS1 = ZeroVector(16);
            Vector elemRHS2 = ZeroVector(16);
            Vector elemRHS3 = ZeroVector(16);
            Matrix elemLHS = ZeroMatrix(16,16);

            Vector condRHS1 = ZeroVector(12);
            Vector condRHS2 = ZeroVector(12);
            Vector condRHS3 = ZeroVector(12);
            Matrix condLHS = ZeroMatrix(12,12);

            for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
                it_node->FastGetSolutionStepValue(DENSITY) = 1000.0;
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_X) = 0.0;
                it_node->FastGetSolutionStepValue(BODY_FORCE_Y) = 0.0;
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = 0.0;
                it_node->SetValue(SLIP_LENGTH, 1.0e10);
            }

            for(unsigned int i=0; i<4; i++){
                for(unsigned int k=0; k<3; k++){
                    for(unsigned int timestep=0; timestep<3; timestep++){
                        // remark: This flow field HAS a wall-normal component.
                        // This choice is made deliberately to see if the tangential projection works
                        pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, timestep)[k] = 1.0;
                        pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, timestep)[k] = 1.0;
                        pElement3->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, timestep)[k] = 1.0;

                        pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, timestep)[k] = 0.0;
                        pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, timestep)[k] = 0.0;
                        pElement3->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, timestep)[k] = 0.0;
                    }
                }
            }

            for(unsigned int i=0; i<4; i++){
                pElement1->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                pElement2->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                pElement3->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
            }

            pElement1->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -10.0;
            pElement1->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -10.0;
            pElement1->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -15.0;
            pElement1->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = -15.0;
            pElement2->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -15.0;
            pElement2->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -10.0;
            pElement2->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -15.0;
            pElement2->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = -20.0;
            pElement3->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -15.0;
            pElement3->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -15.0;
            pElement3->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -10.0;
            pElement3->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = -20.0;

            // Assembly of the residual for node 4 (node between the 3 conditions)
            Vector contriFromElem1 = ZeroVector(3);
            Vector contriFromElem2 = ZeroVector(3);
            Vector contriFromElem3 = ZeroVector(3);

            Vector contriFromCond1Before = ZeroVector(3);
            Vector contriFromCond2Before = ZeroVector(3);
            Vector contriFromCond3Before = ZeroVector(3);
            Vector contriFromCond1After = ZeroVector(3);
            Vector contriFromCond2After = ZeroVector(3);
            Vector contriFromCond3After = ZeroVector(3);

            // Initialization
            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement1->Initialize(r_process_info);
            pElement2->Initialize(r_process_info);
            pElement3->Initialize(r_process_info);

            pElement1->InitializeSolutionStep(r_process_info);
            pElement2->InitializeSolutionStep(r_process_info);
            pElement3->InitializeSolutionStep(r_process_info);

            pElement1->CalculateLocalSystem(elemLHS, elemRHS1, r_process_info);
            pElement2->CalculateLocalSystem(elemLHS, elemRHS2, r_process_info);
            pElement3->CalculateLocalSystem(elemLHS, elemRHS3, r_process_info);

            FindNodalNeighboursProcess find_nodal_neighbours_process(modelPart);
            find_nodal_neighbours_process.Execute();

            NormalCalculationUtils find_nodal_normal_utility;
            find_nodal_normal_utility.CalculateOnSimplex(modelPart, 3);

            pCondition1->Initialize(r_process_info);
            pCondition2->Initialize(r_process_info);
            pCondition3->Initialize(r_process_info);

            // Computing local contributions before setting a smaller slip length
            pCondition1->CalculateLocalSystem(condLHS, condRHS1, r_process_info);
            pCondition2->CalculateLocalSystem(condLHS, condRHS2, r_process_info);
            pCondition3->CalculateLocalSystem(condLHS, condRHS3, r_process_info);

            for (unsigned int i = 0; i < 3; i++){
                contriFromElem1[i] = elemRHS1[12 + i];
                contriFromElem2[i] = elemRHS2[i];
                contriFromElem3[i] = elemRHS3[i];

                contriFromCond1Before[i] = condRHS1[8 + i];
                contriFromCond2Before[i] = condRHS2[8 + i];
                contriFromCond3Before[i] = condRHS3[8 + i];
            }
            const Vector residualAtNodeTwoBefore = 	contriFromElem1 + contriFromElem2 + contriFromElem3 +
                                                    ( contriFromCond1Before + contriFromCond2Before + contriFromCond3Before );

            // change of the slip length
            const double navier_slip_length = 0.0001;
            pElement2->GetGeometry()[0].SetValue(SLIP_LENGTH, navier_slip_length);

            // Computing local contributions before setting a smaller slip length
            pCondition1->CalculateLocalSystem(condLHS, condRHS1, r_process_info);
            pCondition2->CalculateLocalSystem(condLHS, condRHS2, r_process_info);
            pCondition3->CalculateLocalSystem(condLHS, condRHS3, r_process_info);

            for (unsigned int i = 0; i < 3; i++){
                contriFromElem1[i] = elemRHS1[12 + i];
                contriFromElem2[i] = elemRHS2[i];
                contriFromElem3[i] = elemRHS3[i];

                contriFromCond1After[i] = condRHS1[8 + i];
                contriFromCond2After[i] = condRHS2[8 + i];
                contriFromCond3After[i] = condRHS3[8 + i];
            }
            const Vector residualAtNodeTwoAfter = 	contriFromElem1 + contriFromElem2 + contriFromElem3 +
                                                    ( contriFromCond1After + contriFromCond2After + contriFromCond3After );

            const Vector changesInResidual = residualAtNodeTwoAfter - residualAtNodeTwoBefore;
            const double integrationAreaForNode = ( 0.5*10.0*std::sqrt(50.0) + std::sqrt(50.0*50.0+2.0*25.0*25.0) ) / 3.0;

            KRATOS_CHECK_NEAR( changesInResidual[0], -(integrationAreaForNode * 1.0 * 1.0 * 1.0/navier_slip_length), 1e-7 );
            KRATOS_CHECK_NEAR( changesInResidual[1], -(integrationAreaForNode * 1.0 * 1.0 * 1.0/navier_slip_length), 1e-7 );
            KRATOS_CHECK_NEAR( changesInResidual[2], 0.0, 1e-7);
        }

        // Giving a value different to zero to source term in order to test it.

        // /** Checks the TwoFluidNavierStokes2D3N element with a source term in mass conservation equation
        //  * Checks the LHS and RHS for a cut element
        //  */

        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes2D3NError, FluidDynamicsApplicationFastSuite)
        {
             Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);
            

            // Process info creation
            double delta_time = 0.1;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0/(2.0*delta_time);
            bdf_coefs[1] = -2.0/delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR,10.0);


            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
            modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);


            // Define the nodal values
            Matrix vel_original(3, 2);
            vel_original(0, 0) = 0.0;
            vel_original(0, 1) = 0.1;
            vel_original(1, 0) = 0.1;
            vel_original(1, 1) = 0.2;
            vel_original(2, 0) = 0.2;
            vel_original(2, 1) = 0.3;

            // Set the nodal DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node = modelPart.NodesBegin(); it_node < modelPart.NodesEnd(); ++it_node)
            {
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
            }

            for (unsigned int i = 0; i < 3; i++)
            {
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 2; k++)
                {
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9 * vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75 * vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = 0.5;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(9);
            Matrix LHS = ZeroMatrix(9, 9);

            const auto &r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            Vector reference_RHS = ZeroVector(9);

            reference_RHS[0] = 5758.116396;
            reference_RHS[1] = 4790.604956;
            reference_RHS[2] = -28.59258428;
            reference_RHS[3] = -5809.097043;
            reference_RHS[4] = 2150.577447;
            reference_RHS[5] = -15.981999;
            reference_RHS[6] = 394.5059473;
            reference_RHS[7] =-3327.557419;
            reference_RHS[8] =-5.575416721;

            KRATOS_CHECK_VECTOR_NEAR(reference_RHS, RHS, 1e-2);
        }

        // /** Checks the TwoFluidNavierStokes3D4N element with a source term in mass conservation equation
        //  * Checks the LHS and RHS for a cut element
        //  */

        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesCut3D4NError, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);

            // Process info creation
            double delta_time = 0.1;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0/(2.0*delta_time);
            bdf_coefs[1] = -2.0/delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 10.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);
            // Define the nodal values
            Matrix vel_original(4, 3);
            vel_original(0, 0) = 0.0;
            vel_original(0, 1) = 0.1;
            vel_original(0, 2) = 0.2;
            vel_original(1, 0) = 0.1;
            vel_original(1, 1) = 0.2;
            vel_original(1, 2) = 0.3;
            vel_original(2, 0) = 0.2;
            vel_original(2, 1) = 0.3;
            vel_original(2, 2) = 0.4;
            vel_original(3, 0) = 0.3;
            vel_original(3, 1) = 0.4;
            vel_original(3, 2) = 0.5;

            // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node = modelPart.NodesBegin(); it_node < modelPart.NodesEnd(); ++it_node)
            {
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
            }

            for (unsigned int i = 0; i < 4; i++)
            {
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 3; k++)
                {
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9 * vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75 * vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = 1.0;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16, 16);

            const auto &r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            Vector reference_RHS = ZeroVector(16);

            reference_RHS[0] = 3968.830183;
            reference_RHS[1] = 2879.348152;
            reference_RHS[2]=3836.687427;
            reference_RHS[3] = -5.525948818;
            reference_RHS[4] = -3901.342382;
            reference_RHS[5] = 528.2576828;
            reference_RHS[6] = -658.3189566;
            reference_RHS[7] = -4.006423257;
            reference_RHS[8] = 155.8696564;
            reference_RHS[9] = -3271.153765;
            reference_RHS[10] = -357.8352521;
            reference_RHS[11] = -2.860193597;
            reference_RHS[12] = -66.48224427;
            reference_RHS[13]=757.5359384;
            reference_RHS[14] = -4420.514648;
            reference_RHS[15] = -4.374100994;

            KRATOS_CHECK_VECTOR_NEAR(reference_RHS, RHS, 1e-2);
        }

        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesCut3D4NMomentomCorrection, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(3);

            // Variables addition
            modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            modelPart.AddNodalSolutionStepVariable(DENSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            modelPart.AddNodalSolutionStepVariable(PRESSURE);
            modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            modelPart.AddNodalSolutionStepVariable(DISTANCE);

            // Process info creation
            double delta_time = 0.1;
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0/(2.0*delta_time);
            bdf_coefs[1] = -2.0/delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);
            modelPart.GetProcessInfo().SetValue(MOMENTUM_CORRECTION, true);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            auto p_cons_law = Kratos::make_shared<NewtonianTwoFluid3DLaw>();
            pElemProp->SetValue(CONSTITUTIVE_LAW, p_cons_law);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
            modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);
            // Define the nodal values
            Matrix vel_original(4, 3);
            vel_original(0, 0) = 0.0;
            vel_original(0, 1) = 0.1;
            vel_original(0, 2) = 0.2;
            vel_original(1, 0) = 0.1;
            vel_original(1, 1) = 0.2;
            vel_original(1, 2) = 0.3;
            vel_original(2, 0) = 0.2;
            vel_original(2, 1) = 0.3;
            vel_original(2, 2) = 0.4;
            vel_original(3, 0) = 0.3;
            vel_original(3, 1) = 0.4;
            vel_original(3, 2) = 0.5;

            // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node = modelPart.NodesBegin(); it_node < modelPart.NodesEnd(); ++it_node)
            {
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
                it_node->SetValue(DISTANCE_CORRECTION, 0.1);
            }

            for (unsigned int i = 0; i < 4; i++)
            {
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 3; k++)
                {
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9 * vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75 * vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = 1.0;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16, 16);

            const auto &r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            Vector reference_RHS = ZeroVector(16);

            reference_RHS[0] = -99.67824766;
            reference_RHS[1] = 20.15785546;
            reference_RHS[2]= -231.8210039;
            reference_RHS[3] = -0.05448083831;
            reference_RHS[4] = 61.76382984;
            reference_RHS[5] = 31.04354526;
            reference_RHS[6] = -445.5102133;
            reference_RHS[7] = 0.1577814843;
            reference_RHS[8] = 78.28977281;
            reference_RHS[9] = 20.92404867;
            reference_RHS[10] = -435.4151357;
            reference_RHS[11] = 0.006595606811;
            reference_RHS[12] = 115.0227023;
            reference_RHS[13]= 41.93608137;
            reference_RHS[14] = -488.7122322;
            reference_RHS[15] = -0.2098962529;

            KRATOS_CHECK_VECTOR_NEAR(reference_RHS, RHS, 1e-2);
        }

    } // namespace Testing
}  // namespace Kratos.
