/*---------------------------------------------------------------------------*\

Taken from OpenFOAM 1.3 and modified to be compatible with OpenFoam 1.4.1
by David Schmidt


Class
    pressureTransmissiveFvPatchScalarField

Description
    applies a mixed boundary condition to pressure so that it can vary 
    about a fixed value "pInf".

SourceFiles
    pressureTransmissiveFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureTransmissiveFvPatchScalarFields_H
#define pressureTransmissiveFvPatchScalarFields_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class pressureTransmissiveFvPatch Declaration
\*---------------------------------------------------------------------------*/

class pressureTransmissiveFvPatchScalarField
:
    public mixedFvPatchScalarField
{

    // Private data

        //- Pressure of the free stream
        scalar pInf_;

        //- Relaxation length-scale
        scalar lInf_;

        //- Time index
        label curTimeIndex_;

        //- Old time patch field
        scalarField pb0_;

        //- Old-old time patch field
        scalarField pb00_;


public:

    //- Runtime type information
    TypeName("pressureTransmissive");


    // Constructors

        //- Construct from patch and internal field
        pressureTransmissiveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureTransmissiveFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pressureTransmissiveFvPatchScalarField
        //  onto a new patch
        pressureTransmissiveFvPatchScalarField
        (
            const pressureTransmissiveFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new pressureTransmissiveFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureTransmissiveFvPatchScalarField
        (
            const pressureTransmissiveFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(const DimensionedField<scalar, volMesh>& iF) const
        {
            return tmp<fvPatchScalarField>
            (
                new pressureTransmissiveFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the pressure at infinity
            scalar pInf() const
            {
                return pInf_;
            }

            //- Return reference to the pressure at infinity to allow adjustment
            scalar& pInf()
            {
                return pInf_;
            }

            //- Return the relaxation length-scale
            scalar lInf() const
            {
                return lInf_;
            }

            //- Return reference to the relaxation length-scale
            //  to allow adjustment
            scalar& lInf()
            {
                return lInf_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
