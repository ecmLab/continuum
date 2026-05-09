#ifndef ANISOTROPICDIFFUSIVITY_H
#define ANISOTROPICDIFFUSIVITY_H

#include "ADMaterial.h"

// AnisotropicDiffusivity class computes diffusivity based on the gradient of a coupled variable
class AnisotropicDiffusivity : public Material
{
public:
    static InputParameters validParams();
    // Constructor
    AnisotropicDiffusivity(const InputParameters & parameters);
  
protected:
    // Override to compute the property at each quadrature point
    virtual void computeQpProperties() override;
    const Real & _k0;         // Base diffusivity
    const Real & _w;          // Amplitude of the cosine term
    const Real & _lambda;     // Wavelength parameter
    const VariableGradient & _gradient; // Gradient of the coupled variable
    ADMaterialProperty<Real> & _k; // Property to store diffusivity
};

#endif // ANISOTROPICDIFFUSIVITY_H
