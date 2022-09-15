//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//  Collaborators:    Vicente Mataix Ferrandiz
//                    Pablo Becker
//

// System includes

// External includes

// External includes
#include "utilities/math_utils.h"

namespace Kratos 
{

template<class TDataType>
TDataType MathUtils<TDataType>::Abs(const TDataType& rData)
{
    return rData > TDataType(0) ? rData : -rData;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TDataType>
TDataType MathUtils<TDataType>::Min(
    const TDataType& rValue1,
    const TDataType& rValue2
    )
{
    return rValue1 > rValue2 ? rValue2 : rValue1;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TDataType>
TDataType MathUtils<TDataType>::Max(
    const TDataType& rValue1,
    const TDataType& rValue2
    )
{
    return rValue1 > rValue2 ? rValue1 : rValue2;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TDataType>
void MathUtils<TDataType>::GeneralizedInvertMatrix(
    const MatrixType& rInputMatrix,
    MatrixType& rInvertedMatrix,
    TDataType& rInputMatrixDet,
    const TDataType Tolerance
    )
{
    const SizeType size_1 = rInputMatrix.size1();
    const SizeType size_2 = rInputMatrix.size2();

    if (size_1 == size_2) {
        InvertMatrix(rInputMatrix, rInvertedMatrix, rInputMatrixDet, Tolerance);
    } else if (size_1 < size_2) { // Right inverse
        if (rInvertedMatrix.size1() != size_2 || rInvertedMatrix.size2() != size_1) {
            rInvertedMatrix.resize(size_2, size_1, false);
        }
        const Matrix aux = prod(rInputMatrix, trans(rInputMatrix));
        Matrix auxInv;
        InvertMatrix(aux, auxInv, rInputMatrixDet, Tolerance);
        rInputMatrixDet = std::sqrt(rInputMatrixDet);
        noalias(rInvertedMatrix) = prod(trans(rInputMatrix), auxInv);
    } else { // Left inverse
        if (rInvertedMatrix.size1() != size_2 || rInvertedMatrix.size2() != size_1) {
            rInvertedMatrix.resize(size_2, size_1, false);
        }
        const Matrix aux = prod(trans(rInputMatrix), rInputMatrix);
        Matrix auxInv;
        InvertMatrix(aux, auxInv, rInputMatrixDet, Tolerance);
        rInputMatrixDet = std::sqrt(rInputMatrixDet);
        noalias(rInvertedMatrix) = prod(auxInv, trans(rInputMatrix));
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TDataType>
void MathUtils<TDataType>::Solve(
    MatrixType A,
    VectorType& rX,
    const VectorType& rB
    )
{
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
    AMatrix::LUFactorization<MatrixType, DenseVector<std::size_t> > lu_factorization(A);
    double determinant = lu_factorization.determinant();
    KRATOS_ERROR_IF(std::abs(determinant) <= ZeroTolerance) << "Matrix is singular: " << A << std::endl;
    rX = lu_factorization.solve(rB);
#else
    const SizeType size1 = A.size1();
    rX = rB;
    typedef permutation_matrix<SizeType> pmatrix;
    pmatrix pm(size1);
    int singular = lu_factorize(A,pm);
    KRATOS_DEBUG_ERROR_IF(singular == 1) << "Matrix is singular: " << A << std::endl;
    lu_substitute(A, pm, rX);
#endif // ifdef KRATOS_USE_AMATRIX
}

/***********************************************************************************/
/***********************************************************************************/

// Template declarations
template class MathUtils<double>;

} // namespace Kratos.