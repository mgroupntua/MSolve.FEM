namespace MGroup.FEM.Thermal.Line
{
	using System;
	using System.Collections.Generic;
	using System.Diagnostics;

	using MGroup.Constitutive.Thermal;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.MSolve.DataStructures;
	using MGroup.MSolve.Discretization;
	using MGroup.MSolve.Discretization.Dofs;
	using MGroup.MSolve.Discretization.Embedding;
	using MGroup.MSolve.Discretization.Entities;

	/// <summary>
	/// Finite element for heat transfer along a single direction. Can be used for 1D, 2D and 3D problems. Does not take into 
	/// account geometric non-linearities.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class ThermalRod : IThermalElementType, IEmbeddedElement
	{
		private const int numNodes = 2;
		private const int numDofs = 2;
		private static readonly IDofType[][] dofTypes = {
			new IDofType[] { ThermalDof.Temperature }, new IDofType[] { ThermalDof.Temperature } };

		private readonly IThermalProperties material;

		public ThermalRod(IReadOnlyList<INode> nodes, double crossSectionArea, IThermalProperties material)
		{
			Debug.Assert(nodes.Count == 2, "Thermal rod element must have exactly 2 nodes.");
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

		public Matrix BuildCapacityMatrix()
		{
			double kdAL = material.SpecialHeatCoeff * material.Density * CrossSectionArea * Length;
			double[,] capacity = { { kdAL / 3.0, kdAL / 6.0 }, { kdAL / 6.0, kdAL / 3.0 } };
			var matrix = Matrix.CreateFromArray(capacity);
			matrix.MatrixSymmetry = MatrixSymmetry.Symmetric;
			return matrix;
		}

		public Matrix BuildConductivityMatrix()
		{

			double cAoverL = material.ThermalConductivity * CrossSectionArea / Length;
			double[,] conductivity = { { cAoverL, -cAoverL }, { -cAoverL, cAoverL } };
			var matrix = Matrix.CreateFromArray(conductivity);
			matrix.MatrixSymmetry = MatrixSymmetry.Symmetric;
			return matrix;
		}

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofTypes;

		public void ResetConstitutiveLawModified()
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateResponse( double[] localDisplacements)
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

		public IMatrix ConductivityMatrix()
		{
			return DofEnumerator.GetTransformedMatrix(BuildConductivityMatrix());
		}
		public IMatrix PhysicsMatrix()
		{
			return ConductivityMatrix();
		}

		public Dictionary<IDofType, int> GetInternalNodalDOFs(IElementType element, INode node)
		{
			if (node.ID == this.Nodes[0].ID) return new Dictionary<IDofType, int> { { ThermalDof.Temperature, 0 } };
			else if (node.ID == this.Nodes[1].ID) return new Dictionary<IDofType, int> { { ThermalDof.Temperature, 1 } };
			else throw new ArgumentException($"GetInternalNodalDOFs: Node {node.ID} not found in element {element.ID}.");
		}

		public double[] GetLocalDOFValues(IElementType hostElement, double[] hostDOFValues)
		{
			return DofEnumerator.GetTransformedDisplacementsVector(hostDOFValues);
		}
	}
}
