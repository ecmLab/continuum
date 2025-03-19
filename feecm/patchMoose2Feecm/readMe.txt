Comment 1:
all *.c files place are:      moose/modules/tensor_mechanics/src/materials/
one *.h file place is:        moose/modules/tensor_mechanics/include/materials/

Comment 2: Modify ADComputeIncrementalStrainBase.C
Add the following code to function ADComputeStrainBaseTempl:
   _deformation_gradient(                                                                               
                  this->template declareADProperty<RankTwoTensor>(_base_name + "deformation_gradient")),

Add the following code to function ADComputeIncrementalStrainBaseTempl:
   _deformation_gradient[_qp].setToIdentity();