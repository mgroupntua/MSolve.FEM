using System;
using System.Collections.Generic;
using System.Diagnostics;
using MGroup.MSolve.DataStructures;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Embedding;
using MGroup.MSolve.Discretization.Entities;
using MGroup.LinearAlgebra.Implementations;

namespace MGroup.FEM.ConvectionDiffusion.Line
{
	/// <summary>
	/// ConvectionDiffusion element for unknown variable's transfer along a single direction. Can be used for 1D, 2D and 3D problems. Does not take into
	/// account geometric non-linearities.
	/// Authors: Orestis Papas, Christodoulou Theofilos.
	/// </summary>
	public class ConvectionDiffusionRod : IConvectionDiffusionElementType, IEmbeddedElement
	{
		private const int numNodes = 2;
		private const int numDofs = 2;
		private static readonly IDofType[][] dofTypes = {
			new IDofType[] { ConvectionDiffusionDof.UnknownVariable }, new IDofType[] { ConvectionDiffusionDof.UnknownVariable } };

		private readonly IConvectionDiffusionProperties material;

		public ConvectionDiffusionRod(IReadOnlyList<INode> nodes, double crossSectionArea, IConvectionDiffusionProperties material)
		{
			Debug.Assert(nodes.Count == 2, "Convection-Diffusion rod element must have exactly 2 nodes.");
			this.material = material;
			this.Nodes = nodes;
			this.CrossSectionArea = crossSectionArea;
			this.Length = nodes[0].CalculateEuclidianDistanceFrom(nodes[1]);
		}

		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public int ID { get; set; }

		public int SubdomainID { get; set; }

		public IReadOnlyList<INode> Nodes { get; }

		public CellType CellType { get; } = CellType.Line2;

		public double CrossSectionArea { get; }

		public double Length { get; }

		public bool ConstitutiveLawModified => throw new NotImplementedException();

		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		public IList<EmbeddedNode> EmbeddedNodes { get; } = new List<EmbeddedNode>();

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
			double coeff = material.CapacityCoeff * Length;
			double[,] firstTimeDerMatrix = { { coeff / 3d, coeff / 6d }, { coeff / 6d, coeff / 3d } };
			var matrix = Matrix.CreateFromArray(firstTimeDerMatrix);
			matrix.MatrixSymmetry = MatrixSymmetry.Symmetric;
			return matrix;
		}

		public Matrix BuildDiffusionMatrix()
		{
			double coeff = material.DiffusionCoeff * CrossSectionArea / Length;

			double[,] diffusionMatrix = { { coeff, -coeff }, { -coeff, coeff } };
			var matrix = Matrix.CreateFromArray(diffusionMatrix);
			matrix.MatrixSymmetry = MatrixSymmetry.Symmetric;
			return matrix;

		}

		public Matrix BuildConvectionMatrix()
		{
			double coeff = material.ConvectionCoeff[0] * CrossSectionArea;
			double[,] convectionMatrix = { { -coeff / 2d, coeff / 2d }, { -coeff / 2d, coeff / 2d } };
			var matrix = Matrix.CreateFromArray(convectionMatrix);
			matrix.MatrixSymmetry = material.ConvectionCoeff[0] == 0 ? MatrixSymmetry.Symmetric : MatrixSymmetry.NonSymmetric;
			return matrix;
		}

		public Matrix BuildProductionMatrix()
		{
			double coeff = material.DependentSourceCoeff * (-1d) * CrossSectionArea * Length;
			double[,] productionMatrix = { { coeff / 3d, coeff / 6d }, { coeff / 6d, coeff / 3d } };
			var matrix = Matrix.CreateFromArray(productionMatrix);
			matrix.MatrixSymmetry = MatrixSymmetry.Symmetric;
			return matrix;
		}

		public double[] BuildProductionVector()
		{
			double coeff = material.IndependentSourceCoeff * CrossSectionArea * Length;
			double[] productionVector= {coeff / 2d, coeff / 2d };
			return productionVector;
		}

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofTypes;

		public void ResetConstitutiveLawModified() => throw new NotImplementedException();

		public Tuple<double[], double[]> CalculateResponse( double[] localDisplacements) => throw new NotImplementedException();

		public double[] CalculateResponseIntegral() => throw new NotImplementedException();

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements) => throw new NotImplementedException();

		public void SaveConstitutiveLawState(IHaveState externalState) => throw new NotImplementedException();

		public void ClearConstitutiveLawState() => throw new NotImplementedException();

		public void ClearConstitutiveLawStresses() => throw new NotImplementedException();

		public Dictionary<IDofType, int> GetInternalNodalDOFs(IElementType element, INode node)
		{
			if (node.ID == this.Nodes[0].ID) return new Dictionary<IDofType, int> { { ConvectionDiffusionDof.UnknownVariable, 0 } };
			else if (node.ID == this.Nodes[1].ID) return new Dictionary<IDofType, int> { { ConvectionDiffusionDof.UnknownVariable, 1 } };
			else throw new ArgumentException($"GetInternalNodalDOFs: Node {node.ID} not found in element {element.ID}.");
		}

		public double[] GetLocalDOFValues(IElementType hostElement, double[] hostDOFValues)
		{
			return DofEnumerator.GetTransformedDisplacementsVector(hostDOFValues);
		}
	}
}
