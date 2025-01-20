using System;
using System.Collections.Generic;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.MSolve.Numerics.Interpolation.GaussPointExtrapolation;
using MGroup.MSolve.Numerics.Interpolation.Inverse;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Discretization.Embedding;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Geometry.Coordinates;
using MGroup.LinearAlgebra.Implementations;


namespace MGroup.FEM.ConvectionDiffusion.Isoparametric
{
	public class ConvectionDiffusionElement2D : IConvectionDiffusionElementType, IEmbeddedHostElement
	{
		private readonly IDofType[][] dofTypes;
		private readonly IConvectionDiffusionProperties material;

		public ConvectionDiffusionElement2D(double thickness, IReadOnlyList<Node> nodes, IIsoparametricInterpolation2D interpolation,
			IQuadrature2D quadratureForStiffness, IQuadrature2D quadratureForConsistentMass,
			IGaussPointExtrapolation2D gaussPointExtrapolation,
			IConvectionDiffusionProperties material)
		{
			this.material = material;
			this.GaussPointExtrapolation = gaussPointExtrapolation;
			this.Nodes = nodes;
			this.Interpolation = interpolation;
			this.QuadratureForConsistentMass = quadratureForConsistentMass;
			this.QuadratureForStiffness = quadratureForStiffness;
			this.Thickness = thickness;

			dofTypes = new IDofType[nodes.Count][];
			for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new IDofType[] { ConvectionDiffusionDof.UnknownVariable};
		}

		public CellType CellType => Interpolation.CellType;

		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public int ID { get; set; }

		public int SubdomainID { get; set; }

		public IReadOnlyList<INode> Nodes { get; }

		public IGaussPointExtrapolation2D GaussPointExtrapolation { get; }

		public IIsoparametricInterpolation2D Interpolation { get; }

		public IQuadrature2D QuadratureForConsistentMass { get; }

		public IQuadrature2D QuadratureForStiffness { get; }

		public double Thickness { get; }

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
		public double[] ProductionVector()
		{
			return BuildProductionVector();
		}

		public IMatrix PhysicsMatrix()
		{
			return DiffusionMatrix().Add(ConvectionMatrix()).Add(ProductionMatrix());
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
				var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				capacity.AxpyIntoThis(partial, dA);
			}

			capacity.ScaleIntoThis(material.CapacityCoeff * Thickness);
			capacity.MatrixSymmetry = MatrixSymmetry.Symmetric;
			return capacity;
		}



		public Matrix BuildDiffusionMatrix()
		{
			int numDofs = Nodes.Count;
			var diffusion = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
				
				Matrix shapeGradientsCartesian = jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);

				var deformation = Matrix.CreateZero(2, Nodes.Count);
				for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
				{
					deformation[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
					deformation[1, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
				}

				Matrix partialK = deformation.Transpose() * deformation;

				double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
				diffusion.AxpyIntoThis(partialK, dA);
			}
			diffusion.ScaleIntoThis(material.DiffusionCoeff * Thickness);
			diffusion.MatrixSymmetry = MatrixSymmetry.Symmetric;

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
				var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
				
				Matrix shapeGradientsCartesian = jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);

				var deformationX = Matrix.CreateZero(1, Nodes.Count);
				var deformationY = Matrix.CreateZero(1, Nodes.Count);
				for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
				{
					deformationX[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
					deformationY[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
				}

				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);

				Matrix partialConvectionMatrix = shapeFunctionMatrix.Transpose() * deformationX * material.ConvectionCoeff[0]
												+ shapeFunctionMatrix.Transpose() * deformationY * material.ConvectionCoeff[1];

				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;

				convection.AxpyIntoThis(partialConvectionMatrix, dA * Thickness);
			}
			convection.MatrixSymmetry = material.ConvectionCoeff[0] == 0 && material.ConvectionCoeff[1] == 0 ?
				MatrixSymmetry.Symmetric :
				MatrixSymmetry.NonSymmetric;
			return convection;
		}

		public Matrix BuildProductionMatrix()
		{
			int numDofs = Nodes.Count;
			var production= Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
				Matrix partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
				var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				production.AxpyIntoThis(partial, dA);
			}

			production.ScaleIntoThis(material.DependentSourceCoeff * (-1d) * Thickness);
			production.MatrixSymmetry = MatrixSymmetry.Symmetric;
			return production;
		}

		public double[] BuildProductionVector()
		{
			int numDofs = Nodes.Count;
			var productionVector = Matrix.CreateZero(numDofs, 1);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]).Transpose();
				var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
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

		/// <summary>
		/// The shape function matrix is 1-by-n, where n = is the number of shape functions.
		/// </summary>
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

		private readonly IList<INode> embeddedNodes = new List<INode>();

		public EmbeddedNode BuildHostElementEmbeddedNode(IElementType element, INode node, IEmbeddedDOFInHostTransformationVector transformationVector)
		{
			IInverseInterpolation2D inverseInterpolation = Interpolation.CreateInverseMappingFor(Nodes);
			double[] naturalCoordinates = inverseInterpolation.TransformPointCartesianToNatural(new CartesianPoint(node.X, node.Y)).Coordinates;

			if (naturalCoordinates.Length == 0) return null;

			embeddedNodes.Add(node);
			var embeddedNode = new EmbeddedNode(node, element, transformationVector.GetDependentDOFTypes);
			for (int i = 0; i < naturalCoordinates.Length; i++)
				embeddedNode.Coordinates.Add(naturalCoordinates[i]);
			return embeddedNode;
		}

		public double[] GetShapeFunctionsForNode(IElementType element, EmbeddedNode node)
		{
			return Interpolation.EvaluateFunctionsAt(new NaturalPoint(node.Coordinates[0], node.Coordinates[1]));
		}
	}
}
