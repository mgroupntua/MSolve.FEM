using System;
using System.Collections.Generic;
using MGroup.MSolve.DataStructures;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.MSolve.Numerics.Interpolation.GaussPointExtrapolation;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.FEM.ConvectionDiffusion.Isoparametric
{
	public class ConvectionDiffusionElement3D : IConvectionDiffusionElementType
	{
		private readonly static IDofType[] nodalDOFTypes = new IDofType[] { ConvectionDiffusionDof.UnknownVariable};
		private readonly IDofType[][] dofTypes;
		private readonly IConvectionDiffusionProperties material;

		public ConvectionDiffusionElement3D(IReadOnlyList<INode> nodes, IIsoparametricInterpolation3D interpolation,
		IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
		IGaussPointExtrapolation3D gaussPointExtrapolation, IConvectionDiffusionProperties material)
		{
			this.material = material;
			this.GaussPointExtrapolation = gaussPointExtrapolation;
			this.Nodes = nodes;
			this.Interpolation = interpolation;
			this.QuadratureForConsistentMass = quadratureForMass;
			this.QuadratureForStiffness = quadratureForStiffness;

			dofTypes = new IDofType[nodes.Count][];
			for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new IDofType[] { ConvectionDiffusionDof.UnknownVariable };
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

		public IMatrix DiffusionMatrix()
		{
			return BuildDiffusionMatrix();
		}

		public IMatrix ConvectionMatrix()
		{
			return BuildConvectionMatrix();
		}

		public IMatrix ProductionMatrix()
		{
			return BuildProductionMatrix();
		}

		public IMatrix PhysicsMatrix()
		{
			return DiffusionMatrix().Add(ConvectionMatrix()).Add(ProductionMatrix());
		}

		public double[] ProductionVector()
		{
			return BuildProductionVector();
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
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp], 1e-20);
				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				capacity.AxpyIntoThis(partial, dA);
			}
			
			capacity.ScaleIntoThis(material.CapacityCoeff);
			capacity.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;
			return capacity;
		}

		public Matrix BuildDiffusionMatrix()
		{
			int numDofs = Nodes.Count;
			var diffusion = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp], 1e-20);
				Matrix shapeGradientsCartesian =
					jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);

				var deformation = Matrix.CreateZero(3, Nodes.Count);
				for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
				{
					deformation[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
					deformation[1, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
					deformation[2, nodeIdx] = shapeGradientsCartesian[nodeIdx, 2];
				}

				Matrix partialK = deformation.Transpose() * deformation;

				double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
				diffusion.AxpyIntoThis(partialK, dA * material.DiffusionCoeff);
			}
			diffusion.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;

			return diffusion;
		}

		public Matrix BuildConvectionMatrix()
		{
			int numDofs = Nodes.Count;
			var convection = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp], 1e-20);
				
				Matrix shapeGradientsCartesian = jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);

				var deformationX = Matrix.CreateZero(1, Nodes.Count);
				var deformationY = Matrix.CreateZero(1, Nodes.Count);
				var deformationZ = Matrix.CreateZero(1, Nodes.Count);
				for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
				{
					deformationX[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
					deformationY[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
					deformationZ[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 2];
				}

				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
				Matrix partialConvectionMatrix = shapeFunctionMatrix.Transpose() * deformationX * material.ConvectionCoeff[0]
												+ shapeFunctionMatrix.Transpose() * deformationY * material.ConvectionCoeff[1]
												+ shapeFunctionMatrix.Transpose() * deformationZ * material.ConvectionCoeff[2];

				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;

				convection.AxpyIntoThis(partialConvectionMatrix, dA);
			}
			convection.MatrixSymmetry = material.ConvectionCoeff[0] == 0 && material.ConvectionCoeff[1] == 0 && material.ConvectionCoeff[2] == 0 ? LinearAlgebra.Providers.MatrixSymmetry.Symmetric : LinearAlgebra.Providers.MatrixSymmetry.NonSymmetric;


			return convection;
		}
		public Matrix BuildProductionMatrix()
		{
			int numDofs = Nodes.Count;
			var production = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
				Matrix partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp], 1e-20);
				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				production.AxpyIntoThis(partial, dA);
			}

			production.ScaleIntoThis(material.DependentSourceCoeff * (-1d));
			production.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;
			return production;
		}

		public double[] BuildProductionVector()
		{
			int numDofs = Nodes.Count;
			var productionVector = Matrix.CreateZero(numDofs, 1);
			var identity = Matrix.CreateIdentity(numDofs);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]).Transpose();
				var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp], 1e-20);
				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				productionVector.AxpyIntoThis(shapeFunctionMatrix, dA);
			}

			productionVector.ScaleIntoThis(material.IndependentSourceCoeff);

			double[,] vectorDouble = productionVector.CopyToArray2D();
			double[] result = new double[numDofs];

			for (int i = 0; i < numDofs; i++)
			{
				result[i] = vectorDouble[i, 0];
			}
			return result;
		}

		private Matrix BuildDeformationMatrix(Matrix shapeGradientsCartesian)
		{
			var deformation = Matrix.CreateZero(3, Nodes.Count);
			for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
			{
				deformation[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
				deformation[1, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
				deformation[2, nodeIdx] = shapeGradientsCartesian[nodeIdx, 2];
			}
			return deformation;
		}

		public Matrix BuildShapeFunctionMatrix(double[] shapeFunctions)
		{
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
