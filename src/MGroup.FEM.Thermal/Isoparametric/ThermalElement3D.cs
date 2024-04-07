using System;
using System.Collections.Generic;
using MGroup.Constitutive.Thermal;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.MSolve.Numerics.Interpolation.GaussPointExtrapolation;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.DataStructures;

namespace MGroup.FEM.Thermal.Isoparametric
{
	public class ThermalElement3D : IThermalElementType
	{
		private readonly static IDofType[] nodalDOFTypes = new IDofType[] { ThermalDof.Temperature };
		private readonly IDofType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
		private readonly IThermalProperties material;

		public ThermalElement3D(IReadOnlyList<INode> nodes, IIsoparametricInterpolation3D interpolation,
		IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
		IGaussPointExtrapolation3D gaussPointExtrapolation, IThermalProperties material)
		{
			this.material = material;
			this.GaussPointExtrapolation = gaussPointExtrapolation;
			this.Nodes = nodes;
			this.Interpolation = interpolation;
			this.QuadratureForConsistentMass = quadratureForMass;
			this.QuadratureForStiffness = quadratureForStiffness;

			dofTypes = new IDofType[nodes.Count][];
			for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new IDofType[] { ThermalDof.Temperature };
		}
		public CellType CellType => Interpolation.CellType;
		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public IGaussPointExtrapolation3D GaussPointExtrapolation { get; }
		public IIsoparametricInterpolation3D Interpolation { get; }
		public IQuadrature3D QuadratureForConsistentMass { get; }
		public IQuadrature3D QuadratureForStiffness { get; }

		public bool ConstitutiveLawModified => throw new NotImplementedException();

		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		public IMatrix CapacityMatrix()
		{
			return BuildCapacityMatrix();
		}

		public Matrix BuildCapacityMatrix()
		{
			int numDofs = Nodes.Count;
			var capacity = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
				Matrix partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				capacity.AxpyIntoThis(partial, dA);
			}

			//WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
			capacity.Scale(material.Density * material.SpecialHeatCoeff);
			capacity.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;
			return capacity;
		}

		public IMatrix ConductivityMatrix()
		{
			return BuildConductivityMatrix();
		}
		public IMatrix PhysicsMatrix()
		{
			return ConductivityMatrix();
		}
		public Matrix BuildConductivityMatrix()
		{
			int numDofs = Nodes.Count;
			var conductivity = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				// Calculate the necessary quantities for the integration
				//Matrix constitutive = (Matrix)(materialsAtGaussPoints[gp].ConstitutiveMatrix); // ugly cast will be removed along with the legacy Matrix classes
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
				Matrix shapeGradientsCartesian =
					jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
				Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);

				// Contribution of this gauss point to the element stiffness matrix
				Matrix partialK = deformation.Transpose() * deformation;
				//Matrix partialΚ = deformation.Transpose() * (constitutive * deformation);
				//partialK.Scale(materialsAtGaussPoints[gaussPoint].ThermalConductivity);

				double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight; //TODO: this is used by all methods that integrate. I should cache it.
				conductivity.AxpyIntoThis(partialK, dA * material.ThermalConductivity);
				//conductivity.AxpyIntoThis(partialK, dA * 1);
			}

			conductivity.Scale(1);
			conductivity.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;
			return conductivity;
		}

		private Matrix BuildDeformationMatrix(Matrix shapeGradientsCartesian)
		{
			//TODO: isn't this just the transpose of [dNi/dxj]?
			var deformation = Matrix.CreateZero(3, Nodes.Count);
			for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
			{
				deformation[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
				deformation[1, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
				deformation[2, nodeIdx] = shapeGradientsCartesian[nodeIdx, 2];
			}
			return deformation;
		}

		public Matrix BuildShapeFunctionMatrix(double[] shapeFunctions) //TODO: reconsider this. As it is, it just returns the shape functions in a Matrix
		{
			//var array2D = new double[1, shapeFunctions.Length];
			//for (int i = 0; i < shapeFunctions.Length; ++i)
			//{
			//    array2D[0, i] = shapeFunctions[i];
			//}
			//return Matrix.CreateFromArray(array2D);
			return Matrix.CreateFromArray(shapeFunctions, 1, shapeFunctions.Length);
		}

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofTypes;

		public void ResetConstitutiveLawModified()
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateResponse(double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateResponseIntegral()
		{
			throw new NotImplementedException();
		}

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		//public double[] CalculateAccelerationResponseIntegral(IElementType element, IList<MassAccelerationLoad> loads)
		//{
		//	throw new NotImplementedException();
		//}

		public void SaveConstitutiveLawState(IHaveState externalState)
		{
			throw new NotImplementedException();
		}

		public void ClearConstitutiveLawState()
		{
			throw new NotImplementedException();
		}

		public void ClearConstitutiveLawStresses()
		{
			throw new NotImplementedException();
		}
	}
}
