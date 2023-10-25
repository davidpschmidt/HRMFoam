/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "pressureTransmissiveFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicholsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "boundedBackwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pressureTransmissiveFvPatchScalarField::pressureTransmissiveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    pInf_(0.0),
    lInf_(0.0),
    curTimeIndex_(0),
    pb0_(p.size(), 0.0),
    pb00_(p.size(), 0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


pressureTransmissiveFvPatchScalarField::pressureTransmissiveFvPatchScalarField
(
    const pressureTransmissiveFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF ,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    pInf_(ptf.pInf_),
    lInf_(ptf.lInf_),
    curTimeIndex_(ptf.curTimeIndex_),
    pb0_(ptf.pb0_, mapper),
    pb00_(ptf.pb00_, mapper)
{}


pressureTransmissiveFvPatchScalarField::pressureTransmissiveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF ,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    pInf_(readScalar(dict.lookup("pInf"))),
    lInf_(readScalar(dict.lookup("lInf"))),
    curTimeIndex_(0),
    pb0_(p.size()),
    pb00_(p.size())
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(patchInternalField());
    }

    pb0_ = *this;
    pb00_ = *this;
    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (mag(pInf_) < SMALL)
    {
        FatalIOErrorIn
        (
            "pressureTransmissiveFvPatchScalarField::"
            "pressureTransmissiveFvPatchScalarField"
            "(const fvPatch&, const DimensionedField<scalar, volMesh>& iF, "
            "const dictionary&)",
            dict
        )   << "unphysical pInf_ specified (pInf_ = 0)" << endl
            << "    To switch off pressure relaxation specify pInf_ = -1"
            << exit(FatalIOError);
    }

    if (mag(lInf_) < SMALL)
    {
        FatalIOErrorIn
        (
            "pressureTransmissiveFvPatchScalarField::"
            "pressureTransmissiveFvPatchScalarField"
            "(const fvPatch&, const DimensionedField<scalar, volMesh>& iF, "
            "const dictionary&)",
            dict
        )   << "unphysical lInf_ specified (lInf_ = 0)" << endl
            << exit(FatalIOError);
    }
}


pressureTransmissiveFvPatchScalarField::pressureTransmissiveFvPatchScalarField
(
    const pressureTransmissiveFvPatchScalarField& ptpsf,
    const DimensionedField<scalar, volMesh>&  iF
)
:
    mixedFvPatchScalarField(ptpsf, iF),
    pInf_(ptpsf.pInf_),
    lInf_(ptpsf.lInf_),
    curTimeIndex_(0),
    pb0_(ptpsf.pb0_),
    pb00_(ptpsf.pb00_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void pressureTransmissiveFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    pb0_.autoMap(m);
    pb00_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void pressureTransmissiveFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const pressureTransmissiveFvPatchScalarField& ptptf =
        refCast<const pressureTransmissiveFvPatchScalarField>(ptf);

    pb0_.rmap(ptptf.pb0_, addr);
    pb00_.rmap(ptptf.pb00_, addr);
}


// Update the coefficients associated with the patch field
void pressureTransmissiveFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    word ddtScheme(patch().boundaryMesh().mesh().ddtScheme("p"));
    const Time& ldb = db().time();

    // Use time index to save oldTime patchField values
    if (curTimeIndex_ != ldb.timeIndex())
    {
        pb00_ = pb0_;
        pb0_ = *this;
        curTimeIndex_ = ldb.timeIndex();
    }

   // Lookup the velocity and compressibility of the patch
    const fvPatchScalarField& psip =
        patch().lookupPatchField<volScalarField, scalar>("thermo:psi");

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    scalar gamma = 1.3;

    // Calculate the speed of the pressure wave w
    // by summing the component of the velocity normal to the boundary
    // and the speed of sound (sqrt(gamma/psi))

//DPS:  Note that this is not suitable for incompressible flow
// because psi would be infinite

    scalarField w = (patch().nf() & Up) + sqrt(gamma/psip);
        /* This version is a bit unreliable
        ::Foam::max
        (
            (patch().nf() & Up)
          + (sign(patchInternalField() - *this))*sqrt(gamma/psip)
          , 0.0
        );
        */

    // Calculate the pressure wave coefficient alpha (See notes)
    scalarField alpha = w*ldb.deltaT().value()*patch().deltaCoeffs();

    // Calculate the pressure relaxation coefficient k (See notes)
    scalarField k = w*ldb.deltaT().value()/lInf_;


    // Non-reflecting outflow boundary
    // If pInf_ defined setup relaxation to that value.
    // If it is 0 or negative
    if (mag(pInf_) > SMALL)
    {
        if
        (
            ddtScheme == fv::EulerDdtScheme<scalar>::typeName
         || ddtScheme == fv::CrankNicholsonDdtScheme<scalar>::typeName
        )
        {
            refValue() = (pb0_ + k*pInf_)/(1.0 + k);
            valueFraction() = (1.0 + k)/(1.0 + alpha + k);
        }
        else if (ddtScheme == fv::backwardDdtScheme<scalar>::typeName)
        {
            refValue() = (2.0*pb0_ - 0.5*pb00_ + k*pInf_)/(1.5 + k);
            valueFraction() = (1.5 + k)/(1.5 + alpha + k);
        }
        else
        {
            FatalErrorIn
            (
                "pressureTransmissiveFvPatchScalarField::updateCoeffs()"
            )   << "    Unsupported temporal differencing scheme : "
                << ddtScheme
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            ddtScheme == fv::EulerDdtScheme<scalar>::typeName
         || ddtScheme == fv::CrankNicholsonDdtScheme<scalar>::typeName
        )
        {
            refValue() = pb0_;
            valueFraction() = 1.0/(1.0 + alpha);
        }
        else if
        (
            ddtScheme == fv::backwardDdtScheme<scalar>::typeName
         || ddtScheme == fv::boundedBackwardDdtScheme::typeName
        )
        {
            refValue() = (2.0*pb0_ - 0.5*pb00_)/1.5;
            valueFraction() = 1.5/(1.5 + alpha);
        }
        else
        {
            FatalErrorIn
            (
                "pressureTransmissiveFvPatchScalarField::updateCoeffs()"
            )   << "    Unsupported temporal differencing scheme : "
                << ddtScheme
                << exit(FatalError);
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


// Write
void pressureTransmissiveFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("pInf") << pInf_ << token::END_STATEMENT << endl;
    os.writeKeyword("lInf") << lInf_ << token::END_STATEMENT << endl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, pressureTransmissiveFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
